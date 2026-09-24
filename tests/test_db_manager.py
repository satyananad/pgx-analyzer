"""
Unit tests for DatabaseManager persistence and permanent deletion operations.
"""

import os
import pytest
import pandas as pd
from database.db_manager import DatabaseManager


def test_database_manager_deletion():
    test_db = "test_pgx_temp.db"
    if os.path.exists(test_db):
        os.remove(test_db)
        
    db = DatabaseManager(db_path=test_db)
    df = pd.DataFrame([{
        'Sample_ID': 'SAMP_001',
        'Gender_Clean': 'Male',
        'Date_of_Birth': '1990-01-01',
        'Native_Place': 'Delhi',
        'State': 'Delhi',
        'Region': 'North India',
        'Test_Requested': 'CYP2C19',
        'Three_Generations_Residency': 'Yes',
        'CYP2C19*2': 'GA',
        'CYP2C19*3': 'GG',
        'CYP2C19*17': 'CC',
        'Diplotype': '*1/*2',
        'Phenotype': 'Intermediate Metabolizer'
    }, {
        'Sample_ID': 'SAMP_002',
        'Gender_Clean': 'Female',
        'Date_of_Birth': '1992-05-05',
        'Native_Place': 'Mumbai',
        'State': 'Maharashtra',
        'Region': 'West India',
        'Test_Requested': 'CYP2C19',
        'Three_Generations_Residency': 'Yes',
        'CYP2C19*2': 'GG',
        'CYP2C19*3': 'GG',
        'CYP2C19*17': 'CT',
        'Diplotype': '*1/*17',
        'Phenotype': 'Rapid Metabolizer'
    }])

    inserted, skipped = db.insert_batch(df)
    assert inserted == 2
    assert len(db.get_existing_sample_ids()) == 2

    # Delete single sample
    db.delete_sample_by_id('SAMP_001')
    remaining = db.get_existing_sample_ids()
    assert 'SAMP_001' not in remaining
    assert 'SAMP_002' in remaining

    # Clear whole database
    db.clear_database()
    assert len(db.get_existing_sample_ids()) == 0

    if os.path.exists(test_db):
        os.remove(test_db)
