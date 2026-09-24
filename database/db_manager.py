"""
SQLite Database Persistence Layer for Incremental Sample Storage and Scalability.
Supports appending 1,000 to 10,000+ sample batches with duplicate checking.
"""

import sqlite3
import pandas as pd
from typing import Set, List, Tuple, Dict, Any
import os
import logging

logger = logging.getLogger(__name__)


class DatabaseManager:
    def __init__(self, db_path: str = "pharmacogenomics.db"):
        self.db_path = db_path
        self._init_db()

    def _get_connection(self):
        return sqlite3.connect(self.db_path)

    def _init_db(self):
        """Creates table schema if it does not exist."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS samples (
                    sample_id TEXT PRIMARY KEY,
                    gender TEXT,
                    date_of_birth TEXT,
                    native_place TEXT,
                    state TEXT,
                    region TEXT,
                    test_requested TEXT,
                    three_generations TEXT,
                    cyp2c19_2 TEXT,
                    cyp2c19_3 TEXT,
                    cyp2c19_17 TEXT,
                    diplotype TEXT,
                    phenotype TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            conn.commit()

    def get_existing_sample_ids(self) -> Set[str]:
        """Returns set of all sample IDs in database."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT sample_id FROM samples")
            rows = cursor.fetchall()
            return {row[0] for row in rows}

    def insert_batch(self, df: pd.DataFrame) -> Tuple[int, int]:
        """
        Inserts new samples into database.
        Returns tuple of (inserted_count, skipped_duplicate_count).
        """
        existing_ids = self.get_existing_sample_ids()
        
        inserted_count = 0
        skipped_count = 0
        
        records = []
        for _, row in df.iterrows():
            sid = str(row.get('Sample_ID', '')).strip()
            if not sid or sid in existing_ids:
                skipped_count += 1
                continue
                
            records.append((
                sid,
                str(row.get('Gender_Clean', row.get('Gender', ''))),
                str(row.get('Date_of_Birth', '')),
                str(row.get('Native_Place_Clean', row.get('Native_Place', ''))),
                str(row.get('State_Clean', row.get('State', ''))),
                str(row.get('Region', '')),
                str(row.get('Test_Requested', '')),
                str(row.get('Three_Generations_Residency', '')),
                str(row.get('CYP2C19*2', '')),
                str(row.get('CYP2C19*3', '')),
                str(row.get('CYP2C19*17', '')),
                str(row.get('Diplotype', '')),
                str(row.get('Phenotype', ''))
            ))
            existing_ids.add(sid)
            inserted_count += 1

        if records:
            with self._get_connection() as conn:
                cursor = conn.cursor()
                cursor.executemany("""
                    INSERT INTO samples (
                        sample_id, gender, date_of_birth, native_place, state, region,
                        test_requested, three_generations, cyp2c19_2, cyp2c19_3,
                        cyp2c19_17, diplotype, phenotype
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, records)
                conn.commit()

        return inserted_count, skipped_count

    def load_all_samples(self) -> pd.DataFrame:
        """Loads all historical samples into a Pandas DataFrame."""
        with self._get_connection() as conn:
            df = pd.read_sql_query("SELECT * FROM samples", conn)
            if not df.empty:
                # Rename columns back to standard dataframe schema
                df = df.rename(columns={
                    'sample_id': 'Sample_ID',
                    'gender': 'Gender_Clean',
                    'date_of_birth': 'Date_of_Birth',
                    'native_place': 'Native_Place',
                    'state': 'State',
                    'region': 'Region',
                    'test_requested': 'Test_Requested',
                    'three_generations': 'Three_Generations_Residency',
                    'cyp2c19_2': 'CYP2C19*2',
                    'cyp2c19_3': 'CYP2C19*3',
                    'cyp2c19_17': 'CYP2C19*17',
                    'diplotype': 'Diplotype',
                    'phenotype': 'Phenotype'
                })
            return df

    def delete_sample_by_id(self, sample_id: str):
        """Deletes a specific sample record by sample_id from database."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM samples WHERE sample_id = ?", (str(sample_id).strip(),))
            conn.commit()

    def clear_database(self):
        """Clears all records from database."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM samples")
            conn.commit()
            cursor.execute("VACUUM")
