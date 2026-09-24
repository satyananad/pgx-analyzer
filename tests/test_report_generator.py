import pytest
import pandas as pd
from core.qc import DataQCEngine
from core.stratification import DemographicStratifier
from reports.report_generator import ReportGenerator

def test_report_generator():
    df = pd.DataFrame({
        'Sample ID': ['IND-001', 'IND-002', 'IND-003'],
        'Gender': ['Male', 'Female', 'Male'],
        'State': ['Tamil Nadu', 'Delhi', 'Maharashtra'],
        'CYP2C19*2': ['GG', 'GA', 'AA'],
        'CYP2C19*3': ['GG', 'GG', 'GG'],
        'CYP2C19*17': ['CC', 'CT', 'TT']
    })
    qc_engine = DataQCEngine(df)
    clean_df, qc_report = qc_engine.run_qc()
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    
    # 1. Test Class Methods
    excel_bytes = ReportGenerator.export_to_excel(full_results)
    assert len(excel_bytes) > 0
    
    pdf_bytes = ReportGenerator.export_to_pdf(full_results)
    assert len(pdf_bytes) > 0

    # 2. Test Instance Methods
    rep_gen = ReportGenerator(full_results, qc_report)
    assert len(rep_gen.generate_excel()) > 0
    assert len(rep_gen.generate_pdf()) > 0
    assert len(rep_gen.generate_csv()) > 0
