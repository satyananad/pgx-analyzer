"""
Unit tests for CYP2C19 Translator (core.cyp2c19).
"""

import pytest
import pandas as pd
from core.cyp2c19 import CYP2C19Translator


def test_diplotype_assignments():
    # Normal Metabolizer (*1/*1)
    row_nm = pd.Series({'CYP2C19*2': 'GG', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CC'})
    dip, pheno = CYP2C19Translator.assign_sample_diplotype(row_nm)
    assert dip == '*1/*1'
    assert pheno == 'Normal Metabolizer'

    # Poor Metabolizer (*2/*2)
    row_pm = pd.Series({'CYP2C19*2': 'AA', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CC'})
    dip, pheno = CYP2C19Translator.assign_sample_diplotype(row_pm)
    assert dip == '*2/*2'
    assert pheno == 'Poor Metabolizer'

    # Ultrarapid Metabolizer (*17/*17)
    row_um = pd.Series({'CYP2C19*2': 'GG', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'TT'})
    dip, pheno = CYP2C19Translator.assign_sample_diplotype(row_um)
    assert dip == '*17/*17'
    assert pheno == 'Ultrarapid Metabolizer'

    # Indeterminate (*2/*17 or *17/*2)
    row_im = pd.Series({'CYP2C19*2': 'GA', 'CYP2C19*3': 'GG', 'CYP2C19*17': 'CT'})
    dip, pheno = CYP2C19Translator.assign_sample_diplotype(row_im)
    assert dip == '*2/*17'
    assert pheno == 'Indeterminate'
