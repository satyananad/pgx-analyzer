"""
Unit tests for Geographic Population Stratification (North, South, East, West India) and Database persistence.
"""

import pytest
import pandas as pd
from core.stratification import DemographicStratifier
from database.db_manager import DatabaseManager


def test_geographic_stratification():
    data = pd.DataFrame([
        {'Sample_ID': 'S1', 'State_Clean': 'Punjab', 'CYP2C19*2': 'GG', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CC', 'Gender_Clean': 'Male'},
        {'Sample_ID': 'S2', 'State_Clean': 'Andhra Pradesh', 'CYP2C19*2': 'GA', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CC', 'Gender_Clean': 'Female'},
        {'Sample_ID': 'S3', 'State_Clean': 'West Bengal', 'CYP2C19*2': 'AA', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CC', 'Gender_Clean': 'Male'},
        {'Sample_ID': 'S4', 'State_Clean': 'Gujarat', 'CYP2C19*2': 'GG', 'CYP2C19*3': 'GA', 'CYP2C19*17': 'CT', 'Gender_Clean': 'Female'}
    ])
    
    results = DemographicStratifier.stratify_and_analyze(data)
    
    assert 'regional' in results
    regional = results['regional']
    
    assert 'North India' in regional
    assert 'South India' in regional
    assert 'East India' in regional
    assert 'West India' in regional
    
    assert regional['North India']['sample_count'] == 1
    assert regional['South India']['sample_count'] == 1
    assert regional['East India']['sample_count'] == 1
    assert regional['West India']['sample_count'] == 1


def test_database_batch_insertion(tmp_path):
    db_file = tmp_path / "test_pgx.db"
    db = DatabaseManager(str(db_file))
    
    data = pd.DataFrame([
        {'Sample_ID': 'S101', 'Gender_Clean': 'Male', 'State_Clean': 'Punjab', 'Region': 'North India', 'CYP2C19*2': 'GG'},
        {'Sample_ID': 'S102', 'Gender_Clean': 'Female', 'State_Clean': 'Tamil Nadu', 'Region': 'South India', 'CYP2C19*2': 'GA'}
    ])
    
    inserted, skipped = db.insert_batch(data)
    assert inserted == 2
    assert skipped == 0
    
    # Re-inserting same data should skip duplicates
    inserted2, skipped2 = db.insert_batch(data)
    assert inserted2 == 0
    assert skipped2 == 2
