"""
Demographic & Geographic Stratification Engine for Population Pharmacogenomics.
Automatically groups datasets into Overall, Regional (North/South/East/West India),
State-Wise (Individual States), and Gender subgroups, computing full genetic and pharmacogenomic parameters for each.
"""

import pandas as pd
import numpy as np
from typing import Dict, Any, List
import re

from config import STATE_TO_REGION, CITY_TO_STATE, SNP_CONFIG
from core.stats import GeneticStatsEngine
from core.cyp2c19 import CYP2C19Translator


class DemographicStratifier:
    @staticmethod
    def assign_region(row: pd.Series) -> str:
        """
        Assigns geographic region (North India, South India, East India, West India)
        to a sample based on State, Native Place, or pre-classified text.
        """
        state_raw = str(row.get('State_Clean', row.get('State', ''))).strip()
        native_place_raw = str(row.get('Native_Place_Clean', row.get('Native_Place', ''))).strip().lower()
        
        # 1. Direct match on State
        if state_raw in STATE_TO_REGION:
            return STATE_TO_REGION[state_raw]
            
        # 2. Case & whitespace insensitive match
        state_title = state_raw.title()
        for k, v in STATE_TO_REGION.items():
            if k.lower() == state_raw.lower():
                return v

        # 3. Check for direct region string (e.g. "North India ", "South India ")
        state_clean = re.sub(r'\s+', ' ', state_raw).title()
        if 'North' in state_clean: return 'North India'
        if 'South' in state_clean: return 'South India'
        if 'East' in state_clean:  return 'East India'
        if 'West' in state_clean:  return 'West India'
        
        # 4. Predict State from Native Place (City/District gazetteer)
        for city, mapped_state in CITY_TO_STATE.items():
            if city in native_place_raw:
                return STATE_TO_REGION[mapped_state]
                
        return 'Unclassified / Other'

    @classmethod
    def stratify_and_analyze(cls, df: pd.DataFrame) -> Dict[str, Any]:
        """
        Performs end-to-end demographic stratification and statistical calculations.
        Returns detailed results dictionary for Overall, Regions, States, and Genders.
        """
        data = df.copy()
        
        # Ensure Region column is populated
        data['Region'] = data.apply(cls.assign_region, axis=1)
        
        # Apply CYP2C19 translation across full dataset
        data_processed, overall_pgx_summary = CYP2C19Translator.process_dataframe(data)
        
        results = {
            'overall': {},
            'regional': {},
            'state_wise': {},
            'gender': {},
            'processed_dataframe': data_processed
        }
        
        # Helper to compute full stats bundle for any dataframe slice
        def analyze_subgroup(sub_df: pd.DataFrame, group_name: str) -> Dict[str, Any]:
            n_sub = len(sub_df)
            missing_genotypes_sub = 0
            for snp_name in SNP_CONFIG.keys():
                if snp_name in sub_df.columns:
                    missing_genotypes_sub += int(sub_df[snp_name].isna().sum())
                    
            missing_pct_sub = (missing_genotypes_sub / (n_sub * len(SNP_CONFIG))) * 100.0 if n_sub > 0 else 0.0

            sub_results = {
                'sample_count': n_sub,
                'missing_count': missing_genotypes_sub,
                'missing_pct': round(missing_pct_sub, 2),
                'snps': {},
                'cyp2c19_summary': {}
            }
            
            # SNP Stats
            for snp_name in SNP_CONFIG.keys():
                if snp_name in sub_df.columns:
                    sub_results['snps'][snp_name] = GeneticStatsEngine.analyze_snp(sub_df[snp_name], snp_name)
                    
            # CYP2C19 Diplotype & Phenotype Stats
            _, pgx_sum = CYP2C19Translator.process_dataframe(sub_df)
            sub_results['cyp2c19_summary'] = pgx_sum
            
            return sub_results

        # 1. Overall Population Analysis
        results['overall'] = analyze_subgroup(data_processed, 'Overall Population')
        
        # 2. Regional Analysis (North, South, East, West India)
        standard_regions = ['North India', 'South India', 'East India', 'West India']
        for region in standard_regions:
            reg_df = data_processed[data_processed['Region'] == region]
            results['regional'][region] = analyze_subgroup(reg_df, region)
            
        # Catch any unclassified regions
        unclass_df = data_processed[~data_processed['Region'].isin(standard_regions)]
        if len(unclass_df) > 0:
            results['regional']['Unclassified / Other'] = analyze_subgroup(unclass_df, 'Unclassified')

        # 3. State-Wise Analysis (Individual Indian States)
        unique_states = data_processed['State_Clean'].unique()
        for state in unique_states:
            state_str = str(state).strip()
            if state_str and state_str.lower() != 'nan':
                st_df = data_processed[data_processed['State_Clean'] == state]
                if len(st_df) > 0:
                    results['state_wise'][state_str] = analyze_subgroup(st_df, state_str)

        # 4. Gender Analysis (Male, Female, Overall)
        genders = ['Male', 'Female']
        for gender in genders:
            gen_df = data_processed[data_processed['Gender_Clean'] == gender]
            results['gender'][gender] = analyze_subgroup(gen_df, gender)
            
        return results
