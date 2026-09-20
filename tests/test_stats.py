"""
Unit tests for Population Genetics Stats Engine (core.stats).
"""

import pytest
import pandas as pd
from core.stats import GeneticStatsEngine


def test_hwe_calculation():
    # Dataset in exact HWE: 100 GG, 200 GA, 100 AA (N=400, p=0.5, q=0.5)
    genotypes = pd.Series(['GG']*100 + ['GA']*200 + ['AA']*100)
    result = GeneticStatsEngine.analyze_snp(genotypes, 'CYP2C19*2')
    
    assert result['valid_samples'] == 400
    assert result['allele_freqs']['G'] == 0.5
    assert result['allele_freqs']['A'] == 0.5
    assert result['hwe']['chi2_stat'] == 0.0
    assert result['hwe']['p_value'] == 1.0
    assert result['hwe']['interpretation'] == 'In HWE'


def test_hwe_departure():
    # Extreme departure: 200 GG, 0 GA, 200 AA
    genotypes = pd.Series(['GG']*200 + ['AA']*200)
    result = GeneticStatsEngine.analyze_snp(genotypes, 'CYP2C19*2')
    
    assert result['hwe']['chi2_stat'] > 0
    assert result['hwe']['p_value'] < 0.05
    assert result['hwe']['interpretation'] == 'Departure from HWE'
