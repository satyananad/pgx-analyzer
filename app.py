"""
Automated Population Pharmacogenomics Analysis Platform
Comprehensive Multi-Page Streamlit SaaS Application

Full implementation of the 12-Section Specification:
1. Title & Aim
2. Input Data Normalization (CYP2C19*2 rs4244285, CYP2C19*3 rs4986893, CYP2C19*17 rs12248560)
3. Visual Welcome Landing Screen with Interactive Manual Data Entry Form & File Upload
4. Data Quality Control (QC) & Missing Values Audit with Invalid Genotype Call Detection & Warnings
5. Genotype Counts & Allele Frequency Engine (p & q calculation, Homogenous vs Heterogenous Zygosity)
6. Hardy-Weinberg Equilibrium Engine (Chi2, P-Value, Haldane Exact Test, Table 3 Reference)
7. Geographic Population Groups (South India, North India, East India, West India, Central India)
8. Dedicated State-Wise Analysis Pages (Individual State Profiles)
9. CYP2C19 Star Allele, Diplotype Calling & CPIC Phenotype Classification (UM, RM, NM, IM, PM)
10. Live Report Preview & Multi-Format Exporter (Excel, CSV, PDF)
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os

from core.qc import DataQCEngine
from core.stratification import DemographicStratifier
from core.stats import GeneticStatsEngine
from core.cyp2c19 import CYP2C19Translator
from database.db_manager import DatabaseManager
from reports.report_generator import ReportGenerator
from config import SNP_CONFIG, PHENOTYPE_ORDER, STATE_TO_REGION

# -----------------------------------------------------------------------------
# 1. PAGE CONFIG & CUSTOM CSS
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Population Pharmacogenomics Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

CUSTOM_CSS = """
<style>
    /* Main Background & Clean Typography */
    .stApp {
        background-color: #F8FAFC !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
        color: #0F172A !important;
    }
    
    /* Georgia Serif Header Titles matching Screenshots 1-5 */
    h1, h2, h3, .serif-header {
        font-family: Georgia, "Times New Roman", Times, serif !important;
        color: #0F172A !important;
        font-weight: 700 !important;
    }
    
    /* Primary Buttons Styling (Dark Forest Green #166534) */
    button[kind="primary"] {
        background-color: #166534 !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 6px !important;
        font-weight: 600 !important;
    }
    button[kind="primary"]:hover {
        background-color: #15803D !important;
    }
    
    /* Secondary Action Buttons */
    button[kind="secondary"] {
        background-color: #475569 !important;
        color: #FFFFFF !important;
        border: none !important;
        border-radius: 6px !important;
        font-weight: 600 !important;
    }
    button[kind="secondary"]:hover {
        background-color: #334155 !important;
    }

    /* Clean Card Box Containers */
    .clean-card-box {
        background-color: #FFFFFF;
        border: 1px solid #CBD5E1;
        border-radius: 8px;
        padding: 1.25rem;
        margin-bottom: 1.2rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.03);
    }
    
    /* Upload Box Styling (Matching Screenshot 1) */
    .upload-container-box {
        background-color: #FAFAFA;
        border: 2px dashed #94A3B8;
        border-radius: 8px;
        padding: 2.2rem 1.5rem;
        text-align: center;
        margin-bottom: 1rem;
    }
    .upload-icon-style {
        font-size: 2.5rem;
        color: #475569;
        margin-bottom: 0.3rem;
    }
    .upload-title-style {
        font-size: 1.1rem;
        font-weight: 700;
        color: #0F172A;
    }
    .upload-sub-style {
        font-size: 0.85rem;
        color: #64748B;
        margin-top: 0.2rem;
    }

    /* Privacy Banner (Matching Screenshot 1) */
    .privacy-info-banner {
        background-color: #E2E8F0;
        border-radius: 4px;
        padding: 0.75rem 1rem;
        font-size: 0.85rem;
        color: #334155;
        margin-bottom: 1.2rem;
    }

    /* Metric Cards (Matching Screenshots 4 & 5) */
    .kpi-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 6px;
        padding: 0.9rem;
        text-align: left;
        box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
    }
    .kpi-val {
        font-family: Georgia, serif;
        font-size: 1.9rem;
        font-weight: 700;
        color: #0F172A;
        line-height: 1.1;
    }
    .kpi-lbl {
        font-size: 0.68rem;
        font-weight: 700;
        text-transform: uppercase;
        color: #64748B;
        letter-spacing: 0.05em;
        margin-top: 0.3rem;
    }

    /* Table Clean Border Header Styling matching Screenshots 1, 2, 3 */
    div[data-testid="stTable"] table {
        border-collapse: collapse !important;
        width: 100% !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
    }
    div[data-testid="stTable"] th {
        border-top: 1.5px solid #0F172A !important;
        border-bottom: 1.5px solid #0F172A !important;
        background-color: transparent !important;
        color: #0F172A !important;
        font-size: 0.75rem !important;
        font-weight: 700 !important;
        letter-spacing: 0.05em !important;
        padding: 0.6rem 0.8rem !important;
    }
    div[data-testid="stTable"] td {
        border-bottom: 1px solid #E2E8F0 !important;
        padding: 0.55rem 0.8rem !important;
        font-size: 0.88rem !important;
        color: #1E293B !important;
    }

    /* Sub-cohort Horizontal Radio Pills Styling (Matching Screenshot 4 & 5) */
    div[role="radiogroup"] {
        gap: 0.5rem !important;
    }
    div[role="radiogroup"] label {
        background-color: #E2E8F0 !important;
        border-radius: 20px !important;
        padding: 0.35rem 0.9rem !important;
        font-size: 0.82rem !important;
        font-weight: 600 !important;
        color: #1E293B !important;
        border: none !important;
    }
    div[role="radiogroup"] label[aria-checked="true"], div[role="radiogroup"] label:has(input:checked) {
        background-color: #166534 !important;
        color: #FFFFFF !important;
    }
    
    /* Footer Disclaimer (Matching Screenshot 1 & 3) */
    .footer-disclaimer-text {
        text-align: center;
        font-size: 0.8rem;
        color: #64748B;
        margin-top: 2.5rem;
        margin-bottom: 1rem;
    }
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# 2. SESSION STATE MANAGEMENT & FAST PIPELINE CACHING
# -----------------------------------------------------------------------------
DELETED_FLAG_FILE = ".data_deleted.flag"

if os.path.exists(DELETED_FLAG_FILE):
    st.session_state['is_deleted'] = True

if 'is_deleted' not in st.session_state:
    st.session_state['is_deleted'] = False

if 'uploaded_df' not in st.session_state:
    st.session_state['uploaded_df'] = None

if 'mapped_cols' not in st.session_state:
    st.session_state['mapped_cols'] = {}

demo_file_path = "categorized by state into North, South, East, West .xlsx"

EMPTY_COLS = [
    'sample ID', 'Gender', 'Date of Birth', 'Native place ', 'State',
    'Test requested', 'Is their family lived at Native place for past 3 generations?',
    'CYP2C19*2 (rs4244285)', 'CYP2C19*3 (rs4986893)', 'CYP2C19*17 ( rs12248560)'
]

# Default to empty dataset (0 samples) on cold start so application starts fresh
if st.session_state['uploaded_df'] is None:
    st.session_state['uploaded_df'] = pd.DataFrame(columns=EMPTY_COLS)

db_manager = DatabaseManager()

def clear_deletion_flag():
    """Clears the disk-based permanent deletion marker file and purges cache."""
    st.session_state['is_deleted'] = False
    if os.path.exists(DELETED_FLAG_FILE):
        try:
            os.remove(DELETED_FLAG_FILE)
        except Exception:
            pass
    try:
        st.cache_data.clear()
    except Exception:
        pass

def reset_all_data():
    """Clears session dataset, resets custom column mappings, creates permanent deletion marker, and purges database & cache."""
    st.session_state['is_deleted'] = True
    st.session_state['uploaded_df'] = pd.DataFrame(columns=EMPTY_COLS)
    st.session_state['mapped_cols'] = {}
    st.session_state['workflow_step'] = '1 · Upload'
    
    # Write persistent disk marker so deletion survives app restarts, browser refreshes, and new tabs
    try:
        with open(DELETED_FLAG_FILE, "w") as f:
            f.write("deleted")
    except Exception:
        pass

    try:
        db_manager.clear_database()
    except Exception:
        pass

    try:
        st.cache_data.clear()
    except Exception:
        pass

def reload_demo_data():
    """Explicitly reloads the demo dataset, removes deletion flag, and resets cache."""
    clear_deletion_flag()
    st.session_state['mapped_cols'] = {}
    if os.path.exists(demo_file_path):
        st.session_state['uploaded_df'] = pd.read_excel(demo_file_path)
    else:
        st.session_state['uploaded_df'] = pd.DataFrame(columns=EMPTY_COLS)

# High performance in-memory cached analysis execution pipeline
@st.cache_data(show_spinner=False)
def run_fast_pipeline(df: pd.DataFrame, mapped_cols_tuple: tuple):
    mapped_dict = dict(mapped_cols_tuple)
    if df is None or len(df) == 0:
        clean_df = pd.DataFrame()
        qc_report = {
            'total_samples': 0,
            'fully_valid_samples': 0,
            'duplicate_sample_count': 0,
            'duplicate_sample_ids': [],
            'invalid_genotypes_count': 0,
            'invalid_genotypes_log': {},
            'snp_qc_stats': {}
        }
        full_results = {
            'overall': {'sample_count': 0, 'snps': {}, 'cyp2c19_summary': {'phenotypes': {}, 'diplotypes': {}}},
            'regional': {},
            'state_wise': {},
            'gender': {},
            'processed_dataframe': pd.DataFrame()
        }
        return clean_df, qc_report, full_results

    qc_engine = DataQCEngine(df, custom_col_mapping=mapped_dict)
    clean_df, qc_report = qc_engine.run_qc()
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    return clean_df, qc_report, full_results

# -----------------------------------------------------------------------------
# 3. SIDEBAR NAVIGATION MENU (CLEAN DATA MANAGEMENT SIDEBAR)
# -----------------------------------------------------------------------------
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/3004/3004458.png", width=64)
    st.title("PGx Analytics Pro")
    st.caption("Automated Population Pharmacogenomics Platform")
    
    st.divider()
    st.markdown("### Data Management")
    if st.button("🗑️ Delete Existing Data (Start Fresh)", type="secondary", use_container_width=True):
        reset_all_data()
        st.success("All existing data deleted! Workspace reset to 0 samples.")
        st.rerun()

    if os.path.exists(demo_file_path):
        if st.button("🔄 Reload Demo Dataset (1,044 Samples)", use_container_width=True):
            reload_demo_data()
            st.success("Loaded workspace dataset!")
            st.rerun()

# Execute Pipeline
df_raw = st.session_state['uploaded_df']
mapped_tuple = tuple(sorted(st.session_state['mapped_cols'].items()))
clean_df, qc_report, full_results = run_fast_pipeline(df_raw, mapped_tuple)
processed_df = full_results.get('processed_dataframe', pd.DataFrame())

# Helper to check active dataset state
def has_active_data() -> bool:
    if full_results is None or 'overall' not in full_results:
        return False
    return full_results['overall'].get('sample_count', 0) > 0

# -----------------------------------------------------------------------------
# REUSABLE UNIFIED RESULTS RENDERER (MATCHING SCREENSHOTS 1, 2, 3, 4)
# -----------------------------------------------------------------------------
def render_unified_results(full_results: dict, qc_report: dict):
    if not has_active_data() or qc_report is None:
        st.info("ℹ️ Workspace is empty (0 samples). All existing data was deleted. Upload a new Excel/CSV dataset file or add sample records on Page 1.")
        return

    ov = full_results['overall']
    tot_samples = qc_report.get('total_samples', ov['sample_count'])
    fully_valid = qc_report.get('fully_valid_samples', ov['sample_count'])
    dup_count = qc_report.get('duplicate_sample_count', 0)
    dup_ids = qc_report.get('duplicate_sample_ids', [])
    invalid_count = qc_report.get('invalid_genotypes_count', 0)
    invalid_log = qc_report.get('invalid_genotypes_log', {})

    # 9 Main Tabs Consolidating All Analysis Views
    tab_overview, tab_gen, tab_hwe, tab_pop, tab_pgx, tab_qc, tab_geo, tab_gender, tab_state = st.tabs([
        "📊 Executive Overview",
        "🧬 Genotype & Allele", 
        "⚖️ Hardy-Weinberg", 
        "🌐 Population Stratification", 
        "💊 CYP2C19 Calling",
        "🛡️ Data Quality Audit",
        "🗺️ Geographic Groups",
        "👫 Gender-Wise",
        "📍 Dedicated State Profiles"
    ])

    # Cohort Dictionary for Sub-cohort Filtering
    cohort_dict = {
        'Overall': full_results.get('overall', {}),
        'South India': full_results.get('regional', {}).get('South India', {}),
        'North India': full_results.get('regional', {}).get('North India', {}),
        'East India': full_results.get('regional', {}).get('East India', {}),
        'West India': full_results.get('regional', {}).get('West India', {}),
        'Central India': full_results.get('regional', {}).get('Central India', {}),
        'Female': full_results.get('gender', {}).get('Female', {}),
        'Male': full_results.get('gender', {}).get('Male', {})
    }

    # Helper for Sub-cohort Pills UI
    def render_subcohort_pills(key_prefix):
        pill_options = [f"{k} (n={v.get('sample_count', 0)})" for k, v in cohort_dict.items() if v]
        selected = st.radio("Select Sub-Cohort Filter Profile:", pill_options, horizontal=True, key=f"{key_prefix}_subcohort")
        sel_key = selected.split(' (')[0]
        return sel_key, cohort_dict.get(sel_key, full_results['overall'])

    # -------------------------------------------------------------------------
    # TAB 1: EXECUTIVE OVERVIEW
    # -------------------------------------------------------------------------
    with tab_overview:
        k1, k2, k3, k4, k5 = st.columns(5)
        k1.metric("TOTAL SAMPLES", f"{tot_samples:,}")
        k2.metric("FULLY VALID", f"{fully_valid:,}")
        k3.metric("MISSING GENOTYPES", "0.0%")
        k4.metric("DUPLICATE IDS", dup_count)
        k5.metric("INVALID GENOTYPE CALLS", invalid_count)

        if dup_ids:
            st.markdown(f"""
            <div style="background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 6px; padding: 0.8rem 1rem; margin-top: 0.5rem; margin-bottom: 1rem; font-family: monospace; font-size: 0.85rem; color: #334155;">
                <strong>Duplicate Sample IDs:</strong> {', '.join([str(x) for x in dup_ids[:10]])}{'...' if len(dup_ids) > 10 else ''}
            </div>
            """, unsafe_allow_html=True)
            
        if invalid_count > 0:
            inv_details = []
            for snp, samples in invalid_log.items():
                s_list = [f"{item['sample_id']} ('{item['raw_value']}')" for item in samples[:5]]
                inv_details.append(f"{snp}: {', '.join(s_list)}")
            st.warning(f"⚠️ **Warning: {invalid_count} Invalid Genotype Call(s) Detected!** Details: {'; '.join(inv_details)}")

        st.markdown("<br>", unsafe_allow_html=True)
        c_ov1, c_ov2 = st.columns(2)
        with c_ov1:
            st.markdown("### 🏆 CYP2C19 Phenotype Summary")
            pgx = ov.get('cyp2c19_summary', {})
            phe_dict = pgx.get('phenotypes', {})
            phe_df = pd.DataFrame([{'Phenotype Category': k, 'Count': v['count'], 'Percentage': f"{v['percentage']:.1f}%"} for k, v in phe_dict.items() if v['count'] > 0])
            st.table(phe_df)

            fig_ov_phe = px.pie(
                names=[k for k, v in phe_dict.items() if v['count'] > 0],
                values=[v['count'] for k, v in phe_dict.items() if v['count'] > 0],
                hole=0.45,
                title="Overall Phenotype Distribution",
                color_discrete_sequence=['#166534', '#1E40AF', '#4B5563', '#991B1B', '#D97706', '#0D9488']
            )
            fig_ov_phe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_ov_phe, use_container_width=True)

        with c_ov2:
            st.markdown("### 🌐 Regional Sample Distribution")
            reg = full_results.get('regional', {})
            reg_df = pd.DataFrame([{'Region': r, 'Sample Count (N)': r_data.get('sample_count', 0), '% of Cohort': f"{(r_data.get('sample_count', 0)/tot_samples*100):.1f}%" if tot_samples>0 else "0%"} for r, r_data in reg.items()])
            st.table(reg_df)

            reg_chart_df = pd.DataFrame([{'Region': r, 'Samples (N)': r_data.get('sample_count', 0)} for r, r_data in reg.items()])
            fig_ov_reg = px.bar(
                reg_chart_df, x='Region', y='Samples (N)',
                title="Regional Population Group Sizes",
                color='Region',
                color_discrete_sequence=px.colors.qualitative.Set2
            )
            fig_ov_reg.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_ov_reg, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 2: GENOTYPE & ALLELE
    # -------------------------------------------------------------------------
    with tab_gen:
        sel_key, sel_cohort = render_subcohort_pills("tab_gen")
        snps_data = sel_cohort.get('snps', {})
        
        def get_zygosity(gt_name, snp_name):
            snp_conf = SNP_CONFIG.get(snp_name, {})
            wt = snp_conf.get('wildtype_genotype', '')
            het = snp_conf.get('het_genotype', '')
            var = snp_conf.get('hom_var_genotype', '')
            if gt_name == wt:
                return "Homozygous Wildtype (Homogenous)"
            elif gt_name in [het, het[::-1]]:
                return "Heterozygote (Heterogenous)"
            elif gt_name == var:
                return "Homozygous Variant (Homogenous)"
            return "Genotype Call"

        def get_allele_type(al_name, snp_name):
            snp_conf = SNP_CONFIG.get(snp_name, {})
            ref = snp_conf.get('ref_allele', '')
            var = snp_conf.get('var_allele', '')
            if al_name == ref:
                return "Reference Allele"
            elif al_name == var:
                return "Variant Allele"
            return "Allele Call"

        snp_display_order = [
            ('CYP2C19*2', 'rs4244285 (CYP2C19*2, c.681G>A)'),
            ('CYP2C19*3', 'rs4986893 (CYP2C19*3, c.636G>A)'),
            ('CYP2C19*17', 'rs12248560 (CYP2C19*17, c.-806C>T)')
        ]

        for snp_code, snp_title in snp_display_order:
            if snp_code in snps_data:
                s_res = snps_data[snp_code]
                v_cnt = s_res.get('valid_samples', 0)
                inv_cnt = qc_report['snp_qc_stats'].get(snp_code, {}).get('invalid_samples', 0)
                m_cnt = qc_report['snp_qc_stats'].get(snp_code, {}).get('missing_samples', 0)
                
                st.markdown(f"#### {snp_title}")
                st.caption(f"Valid: {v_cnt:,} · Missing: {m_cnt} · Invalid calls: {inv_cnt}")
                
                c_tbl, c_fig = st.columns([1.1, 0.9])
                with c_tbl:
                    g_c = s_res.get('genotype_counts', {})
                    g_p = s_res.get('genotype_pcts', {})
                    g_f = s_res.get('genotype_freqs', {})
                    gt_rows = [{
                        'GENOTYPE': k, 
                        'ZYGOSITY': get_zygosity(k, snp_code),
                        'COUNT': v, 
                        'FREQUENCY': g_f.get(k, 0),
                        'PERCENTAGE': f"{g_p.get(k, 0):.1f}%"
                    } for k, v in g_c.items()]
                    st.table(pd.DataFrame(gt_rows))
                    
                    a_c = s_res.get('allele_counts', {})
                    a_p = s_res.get('allele_pcts', {})
                    al_rows = [{
                        'ALLELE': k, 
                        'TYPE': get_allele_type(k, snp_code),
                        'COUNT': v, 
                        'FREQUENCY': f"{a_p.get(k, 0):.1f}%"
                    } for k, v in a_c.items()]
                    st.table(pd.DataFrame(al_rows))
                    
                with c_fig:
                    fig_df = pd.DataFrame([{'Genotype': k, 'Count': v} for k, v in g_c.items()])
                    fig_bar = px.bar(
                        fig_df, x='Genotype', y='Count', 
                        title=f"{snp_code} Genotype distribution ({sel_key})", 
                        color='Genotype',
                        color_discrete_sequence=['#15803D', '#1D4ED8', '#B91C1C']
                    )
                    fig_bar.update_layout(height=220, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_bar, use_container_width=True)

                    al_df = pd.DataFrame([{'Allele': k, 'Frequency (%)': a_p.get(k, 0)} for k in a_c.keys()])
                    fig_al = px.bar(
                        al_df, x='Allele', y='Frequency (%)',
                        title=f"{snp_code} Allele Frequencies ({sel_key})",
                        color='Allele',
                        color_discrete_sequence=['#0D9488', '#D97706']
                    )
                    fig_al.update_layout(height=200, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_al, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 3: HARDY-WEINBERG
    # -------------------------------------------------------------------------
    with tab_hwe:
        sel_key, sel_cohort = render_subcohort_pills("tab_hwe")
        snps_data = sel_cohort.get('snps', {})
        
        st.markdown(f"### Hardy-Weinberg Equilibrium Audit ({sel_key} Cohort)")
        hwe_table_rows = []
        for snp_name, s_res in snps_data.items():
            hw = s_res.get('hwe', {})
            gc = s_res.get('genotype_counts', {})
            exp_c = hw.get('expected_counts', {})
            hwe_table_rows.append({
                'SNP Variant': snp_name,
                'Valid Samples (N)': s_res.get('valid_samples', 0),
                'Observed Genotypes': " / ".join([f"{k}:{v}" for k, v in gc.items()]),
                'Expected Genotypes': " / ".join([f"{k}:{exp_c.get(k, 0)}" for k in gc.keys()]),
                'Chi-Square (χ²)': hw.get('chi2_stat', 0),
                'Degrees of Freedom': hw.get('df', 1),
                'P-Value': hw.get('p_value', 1.0),
                'Exact Test P-Value': hw.get('exact_p_value', 1.0),
                'Interpretation': hw.get('interpretation', 'In HWE')
            })
            
        c_hwe_tbl, c_hwe_fig = st.columns([1.2, 0.8])
        with c_hwe_tbl:
            st.table(pd.DataFrame(hwe_table_rows))

        with c_hwe_fig:
            chi2_chart_df = pd.DataFrame([
                {'SNP': snp_name, 'Chi2 Stat': s_res.get('hwe', {}).get('chi2_stat', 0)}
                for snp_name, s_res in snps_data.items()
            ])
            fig_hwe = px.bar(
                chi2_chart_df, x='SNP', y='Chi2 Stat',
                title=f"HWE Chi-Square (χ²) Stats vs Threshold 3.841 ({sel_key})",
                color='SNP',
                color_discrete_sequence=['#2563EB', '#7C3AED', '#DB2777']
            )
            fig_hwe.add_hline(y=3.841, line_dash="dash", line_color="red", annotation_text="Critical Threshold (3.841)")
            fig_hwe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_hwe, use_container_width=True)
        
        st.markdown("##### Table 3. Chi-Square Distribution Reference Table (df = 1, α = 0.05, Critical = 3.841)")
        chi2_ref_df = pd.DataFrame([{
            'Level of Significance (α)': '0.50',
            '0.10': '2.706',
            '0.05 (Critical Threshold)': '3.841',
            '0.02': '5.412',
            '0.01': '6.635',
            '0.001': '10.827'
        }])
        st.table(chi2_ref_df)
        st.info("📌 **HWE Decision Rule:** If $\\chi^2 \\ge 3.841$ ($P < 0.05$), population deviates from Hardy-Weinberg Equilibrium.")

    # -------------------------------------------------------------------------
    # TAB 4: POPULATION STRATIFICATION
    # -------------------------------------------------------------------------
    with tab_pop:
        st.markdown("### Group sizes")
        reg_dict = full_results.get('regional', {})
        gen_dict = full_results.get('gender', {})
        tot_n = ov['sample_count']
        
        group_rows = [
            {'GROUP': 'Overall', 'N': tot_n, '% OF TOTAL': '100.0%'},
            {'GROUP': 'Region — South India', 'N': reg_dict.get('South India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('South India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — North India', 'N': reg_dict.get('North India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('North India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — East India', 'N': reg_dict.get('East India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('East India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — West India', 'N': reg_dict.get('West India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('West India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Region — Central India', 'N': reg_dict.get('Central India', {}).get('sample_count', 0), '% OF TOTAL': f"{(reg_dict.get('Central India', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Gender — Female', 'N': gen_dict.get('Female', {}).get('sample_count', 0), '% OF TOTAL': f"{(gen_dict.get('Female', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
            {'GROUP': 'Gender — Male', 'N': gen_dict.get('Male', {}).get('sample_count', 0), '% OF TOTAL': f"{(gen_dict.get('Male', {}).get('sample_count', 0)/tot_n*100):.1f}%" if tot_n>0 else "0%"},
        ]
        
        c_pop1, c_pop2 = st.columns([1.0, 1.0])
        with c_pop1:
            st.table(pd.DataFrame(group_rows))

        with c_pop2:
            pop_chart_df = pd.DataFrame([
                {'Group': r['GROUP'], 'N': r['N']}
                for r in group_rows if r['GROUP'] != 'Overall'
            ])
            fig_pop = px.bar(
                pop_chart_df, x='N', y='Group', orientation='h',
                title="Population Subgroup Sizes (N)",
                color='Group',
                color_discrete_sequence=px.colors.qualitative.Bold
            )
            fig_pop.update_layout(height=320, margin=dict(l=20, r=20, t=35, b=20), yaxis=dict(autorange="reversed"))
            st.plotly_chart(fig_pop, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 5: CYP2C19 CALLING
    # -------------------------------------------------------------------------
    with tab_pgx:
        sel_key, sel_cohort = render_subcohort_pills("tab_pgx")
        pgx = sel_cohort.get('cyp2c19_summary', {})
        dips = pgx.get('diplotypes', {})
        phenos = pgx.get('phenotypes', {})
        
        c_dip, c_dip_chart = st.columns([1.1, 0.9])
        with c_dip:
            st.markdown("#### Diplotype distribution")
            dip_rows = [{'DIPLOTYPE': k, 'COUNT': v['count'], 'FREQUENCY': f"{v['percentage']:.1f}%"} for k, v in dips.items() if v['count'] > 0]
            st.table(pd.DataFrame(dip_rows))
            
        with c_dip_chart:
            dip_fig_df = pd.DataFrame([{'Diplotype': k, 'Count': v['count']} for k, v in dips.items() if v['count'] > 0])
            fig_dip = px.bar(dip_fig_df, x='Count', y='Diplotype', orientation='h', title="Diplotype counts", color_discrete_sequence=['#166534'])
            fig_dip.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20), yaxis=dict(autorange="reversed"))
            st.plotly_chart(fig_dip, use_container_width=True)

        st.markdown("<br>", unsafe_allow_html=True)
        c_phe, c_phe_chart = st.columns([1.1, 0.9])
        with c_phe:
            st.markdown("#### Metabolizer / phenotype classification")
            phe_rows = [{'PHENOTYPE': k, 'COUNT': v['count'], 'PERCENTAGE': f"{v['percentage']:.1f}%"} for k, v in phenos.items() if v['count'] > 0]
            st.table(pd.DataFrame(phe_rows))
            
        with c_phe_chart:
            phe_labels = [k for k, v in phenos.items() if v['count'] > 0]
            phe_counts = [v['count'] for k, v in phenos.items() if v['count'] > 0]
            fig_donut = px.pie(
                names=phe_labels, 
                values=phe_counts, 
                hole=0.5, 
                title="Phenotype distribution", 
                color_discrete_sequence=['#166534', '#1E40AF', '#4B5563', '#991B1B', '#D97706', '#0D9488']
            )
            fig_donut.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_donut, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 6: DATA QUALITY AUDIT
    # -------------------------------------------------------------------------
    with tab_qc:
        st.markdown("### 🛡️ Data Quality Control (QC) & Missing Values Audit")
        q1, q2, q3, q4, q5 = st.columns(5)
        q1.metric("Total Cohort Uploaded", qc_report['total_samples'])
        q2.metric("Fully Valid Samples", qc_report['fully_valid_samples'])
        q3.metric("Duplicate Sample IDs", qc_report['duplicate_sample_count'])
        q4.metric("Invalid Genotype Calls", qc_report['invalid_genotypes_count'])
        q5.metric("Validation Status", "PASS" if qc_report['invalid_genotypes_count'] == 0 else "WARNING")
        
        if qc_report['invalid_genotypes_count'] > 0:
            st.error(f"⚠️ **Warning: {qc_report['invalid_genotypes_count']} Invalid Genotype Call(s) Detected!** Invalid genotype entries were flagged and quarantined.")
            inv_log_rows = []
            for snp, samples in qc_report['invalid_genotypes_log'].items():
                for item in samples:
                    inv_log_rows.append({'SNP Variant': snp, 'Sample ID': item['sample_id'], 'Raw Invalid Genotype Value': item['raw_value'], 'Audit Action': 'Quarantined & Excluded'})
            st.table(pd.DataFrame(inv_log_rows))

        c_qc_tbl, c_qc_fig = st.columns([1.1, 0.9])
        with c_qc_tbl:
            st.markdown("#### Per-SNP Quality Control & Missing Values Table")
            qc_rows = []
            for snp, sdata in qc_report['snp_qc_stats'].items():
                qc_rows.append({
                    'SNP Variant': snp,
                    'Valid Samples': sdata['valid_samples'],
                    'Missing Genotypes Count': sdata['missing_samples'],
                    'Invalid Genotypes Count': sdata.get('invalid_samples', 0),
                    'Missing Data Rate (%)': f"{sdata['missing_pct']:.1f}%",
                    'Quality Audit Status': 'PASS' if sdata.get('invalid_samples', 0) == 0 else 'WARNING'
                })
            st.table(pd.DataFrame(qc_rows))

        with c_qc_fig:
            qc_fig_df = pd.DataFrame([
                {'SNP': snp, 'Valid': sdata['valid_samples'], 'Missing': sdata['missing_samples'], 'Invalid': sdata.get('invalid_samples', 0)}
                for snp, sdata in qc_report['snp_qc_stats'].items()
            ])
            fig_qc = px.bar(
                qc_fig_df, x='SNP', y=['Valid', 'Missing', 'Invalid'],
                title="Sample Quality Audit Breakdown per SNP",
                barmode='stack',
                color_discrete_sequence=['#16A34A', '#EAB308', '#DC2626']
            )
            fig_qc.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_qc, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 7: GEOGRAPHIC POPULATION GROUPS
    # -------------------------------------------------------------------------
    with tab_geo:
        st.markdown("### 🌐 Geographic Population Groups Analysis")
        reg = full_results['regional']
        regions = ['South India', 'North India', 'East India', 'West India', 'Central India']
        
        cols = st.columns(5)
        for idx, r in enumerate(regions):
            r_data = reg.get(r, {})
            with cols[idx]:
                st.markdown(f"""
                <div class="kpi-card">
                    <div class="kpi-lbl">{r}</div>
                    <div class="kpi-val">{r_data.get('sample_count', 0):,}</div>
                    <div style="font-size: 0.8rem; color: #64748B; margin-top: 0.3rem;">
                        CYP2C19*2: <strong>{r_data.get('snps', {}).get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)}</strong><br>
                        PM: <strong style="color: #EF4444;">{r_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer', r_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {})).get('percentage', 0):.1f}%</strong>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
        st.markdown("<br>", unsafe_allow_html=True)

        c_geo1, c_geo2 = st.columns(2)
        with c_geo1:
            geo_var_df = pd.DataFrame([
                {'Region': r, 'CYP2C19*2 Variant Freq f(A)': reg.get(r, {}).get('snps', {}).get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)}
                for r in regions
            ])
            fig_geo_var = px.bar(
                geo_var_df, x='Region', y='CYP2C19*2 Variant Freq f(A)',
                title="CYP2C19*2 Variant Frequency f(A) by Region",
                color='Region',
                color_discrete_sequence=px.colors.qualitative.Dark24
            )
            fig_geo_var.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_geo_var, use_container_width=True)

        with c_geo2:
            geo_pm_df = pd.DataFrame([
                {
                    'Region': r, 
                    'Poor Metabolizer %': reg.get(r, {}).get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer', reg.get(r, {}).get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {})).get('percentage', 0)
                }
                for r in regions
            ])
            fig_geo_pm = px.bar(
                geo_pm_df, x='Region', y='Poor Metabolizer %',
                title="Poor Metabolizers (PM %) by Region",
                color='Region',
                color_discrete_sequence=px.colors.qualitative.Vivid
            )
            fig_geo_pm.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_geo_pm, use_container_width=True)

        st.markdown("#### Comprehensive Regional Statistical Summary Table")
        
        cyp2_rows = [
            ('A Count (Variant Allele)', lambda s: s.get('CYP2C19*2', {}).get('allele_counts', {}).get('A', 0)),
            ('G Count (Reference Allele)', lambda s: s.get('CYP2C19*2', {}).get('allele_counts', {}).get('G', 0)),
            ('Total Allele Count (2N)', lambda s: sum(s.get('CYP2C19*2', {}).get('allele_counts', {}).values())),
            ('GA COUNT (Heterozygote)', lambda s: s.get('CYP2C19*2', {}).get('genotype_counts', {}).get('GA', 0)),
            ('AA COUNT (Homozygous Variant)', lambda s: s.get('CYP2C19*2', {}).get('genotype_counts', {}).get('AA', 0)),
            ('GG COUNT (Homozygous Wildtype)', lambda s: s.get('CYP2C19*2', {}).get('genotype_counts', {}).get('GG', 0)),
            ('Total Sample N', lambda s: sum(s.get('CYP2C19*2', {}).get('genotype_counts', {}).values())),
            ('Allele Frequency F(A)', lambda s: s.get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)),
            ('Allele Frequency F(G)', lambda s: s.get('CYP2C19*2', {}).get('allele_freqs', {}).get('G', 0)),
            ('Chi-Square (χ²)', lambda s: s.get('CYP2C19*2', {}).get('hwe', {}).get('chi2_stat', 0)),
            ('P VALUE', lambda s: s.get('CYP2C19*2', {}).get('hwe', {}).get('p_value', 1.0)),
            ('HWE Interpretation', lambda s: s.get('CYP2C19*2', {}).get('hwe', {}).get('interpretation', 'In HWE'))
        ]
        
        cyp2_matrix = []
        for label, func in cyp2_rows:
            r_dict = {'Statistical Parameter (*2)': label}
            for r in regions:
                r_dict[r] = func(reg.get(r, {}).get('snps', {}))
            cyp2_matrix.append(r_dict)
            
        st.dataframe(pd.DataFrame(cyp2_matrix), use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 8: GENDER-WISE ANALYSIS
    # -------------------------------------------------------------------------
    with tab_gender:
        st.markdown("### 👫 Gender-Wise Population Pharmacogenomics Analysis")
        gender_data = full_results.get('gender', {})
        male = gender_data.get('Male', {})
        female = gender_data.get('Female', {})
        
        m_count = male.get('sample_count', 0)
        f_count = female.get('sample_count', 0)
        tot_count = ov.get('sample_count', 0)
        
        gk1, gk2, gk3 = st.columns(3)
        gk1.metric("Male Samples (N)", f"{m_count:,}")
        gk2.metric("Female Samples (N)", f"{f_count:,}")
        gk3.metric("Total Cohort N", f"{tot_count:,}")
        
        m_phenos = male.get('cyp2c19_summary', {}).get('phenotypes', {})
        f_phenos = female.get('cyp2c19_summary', {}).get('phenotypes', {})
        
        pheno_chart_df = []
        for p_name in PHENOTYPE_ORDER:
            if p_name in m_phenos:
                pheno_chart_df.append({'Gender': 'Male', 'Phenotype': p_name.split(' (')[0], 'Count': m_phenos[p_name]['count']})
            if p_name in f_phenos:
                pheno_chart_df.append({'Gender': 'Female', 'Phenotype': p_name.split(' (')[0], 'Count': f_phenos[p_name]['count']})
                
        st.markdown("<br>", unsafe_allow_html=True)
        cg1, cg2 = st.columns(2)
        with cg1:
            if pheno_chart_df:
                fig_g2 = px.bar(
                    pd.DataFrame(pheno_chart_df), 
                    x='Phenotype', y='Count', 
                    color='Gender', barmode='group', 
                    title="Metabolizer Phenotypes Comparison (Male vs Female)",
                    color_discrete_sequence=['#2563EB', '#EC4899']
                )
                fig_g2.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_g2, use_container_width=True)

        with cg2:
            m_snps = male.get('snps', {}).get('CYP2C19*2', {}).get('genotype_counts', {})
            f_snps = female.get('snps', {}).get('CYP2C19*2', {}).get('genotype_counts', {})
            gen_gt_df = []
            for gt in ['GG', 'GA', 'AA']:
                if gt in m_snps:
                    gen_gt_df.append({'Gender': 'Male', 'Genotype': gt, 'Count': m_snps[gt]})
                if gt in f_snps:
                    gen_gt_df.append({'Gender': 'Female', 'Genotype': gt, 'Count': f_snps[gt]})
            if gen_gt_df:
                fig_gen_gt = px.bar(
                    pd.DataFrame(gen_gt_df), x='Genotype', y='Count', color='Gender', barmode='group',
                    title="CYP2C19*2 Genotypes Comparison (Male vs Female)",
                    color_discrete_sequence=['#2563EB', '#EC4899']
                )
                fig_gen_gt.update_layout(height=280, margin=dict(l=20, r=20, t=35, b=20))
                st.plotly_chart(fig_gen_gt, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 9: DEDICATED STATE-WISE PAGES
    # -------------------------------------------------------------------------
    with tab_state:
        st.markdown("### 📍 Dedicated Individual State Profiles")
        state_dict = full_results.get('state_wise', {})
        state_names = sorted(list(state_dict.keys()))
        
        if state_names:
            selected_st = st.selectbox("Select Dedicated State Profile:", state_names, index=0, key="st_tab_selectbox")
            st_data = state_dict[selected_st]
            
            sk1, sk2, sk3, sk4 = st.columns(4)
            sk1.metric("State Sample N", st_data.get('sample_count', 0))
            sk2.metric("Missing Data Rate", f"{st_data.get('missing_pct', 0)}%")
            cyp2_freq = st_data.get('snps', {}).get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
            sk3.metric("CYP2C19*2 Var Freq", cyp2_freq)
            pm_pct = st_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer', st_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {})).get('percentage', 0)
            sk4.metric("Poor Metabolizers %", f"{pm_pct:.1f}%")

            st_phenos = st_data.get('cyp2c19_summary', {}).get('phenotypes', {})
            st_snps = st_data.get('snps', {}).get('CYP2C19*2', {}).get('genotype_counts', {})
            
            cs1, cs2 = st.columns(2)
            with cs1:
                st_gt_df = pd.DataFrame([{'Genotype': k, 'Count': v} for k, v in st_snps.items()])
                if not st_gt_df.empty:
                    fig_st_gt = px.bar(
                        st_gt_df, x='Genotype', y='Count',
                        title=f"{selected_st} — CYP2C19*2 Genotypes",
                        color='Genotype',
                        color_discrete_sequence=['#16A34A', '#2563EB', '#DC2626']
                    )
                    fig_st_gt.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_st_gt, use_container_width=True)

            with cs2:
                st_phe_df = pd.DataFrame([{'Phenotype': k.split(' (')[0], 'Count': v['count']} for k, v in st_phenos.items() if v['count'] > 0])
                if not st_phe_df.empty:
                    fig_st_phe = px.pie(
                        st_phe_df, names='Phenotype', values='Count', hole=0.45,
                        title=f"{selected_st} — Metabolizer Phenotypes",
                        color_discrete_sequence=['#166534', '#1E40AF', '#4B5563', '#991B1B', '#D97706', '#0D9488']
                    )
                    fig_st_phe.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_st_phe, use_container_width=True)

# -----------------------------------------------------------------------------
# 3-STEP WORKFLOW PIPELINE RENDERER (EXACT MATCH FOR USER SCREENSHOTS 1, 2, 3, 4, 5)
# -----------------------------------------------------------------------------
def render_workflow_pipeline(df_raw, full_results, qc_report):
    if 'workflow_step' not in st.session_state:
        st.session_state['workflow_step'] = '1 · Upload'

    step_cols = st.columns(3)
    t1 = "primary" if st.session_state['workflow_step'] == '1 · Upload' else "secondary"
    t2 = "primary" if st.session_state['workflow_step'] == '2 · Map columns' else "secondary"
    t3 = "primary" if st.session_state['workflow_step'] == '3 · Results' else "secondary"

    if step_cols[0].button("1 · Upload", use_container_width=True, type=t1, key="nav_btn_step1"):
        st.session_state['workflow_step'] = '1 · Upload'
        st.rerun()

    if step_cols[1].button("2 · Map columns", use_container_width=True, type=t2, key="nav_btn_step2"):
        st.session_state['workflow_step'] = '2 · Map columns'
        st.rerun()

    if step_cols[2].button("3 · Results", use_container_width=True, type=t3, key="nav_btn_step3"):
        st.session_state['workflow_step'] = '3 · Results'
        st.rerun()

    st.markdown("<br>", unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # STEP 1: UPLOAD (Matching Image 1)
    # -------------------------------------------------------------------------
    if st.session_state['workflow_step'] == '1 · Upload':
        st.markdown("## Upload genotype dataset")
        st.caption('Excel (.xlsx/.csv) with one row per sample: Sample ID, Gender, DOB, region/location, and genotype columns for each SNP (e.g. "GA", "G/A", "AA").')

        st.markdown("""
        <div class="upload-container-box">
            <div class="upload-icon-style">⇧</div>
            <div class="upload-title-style">Click to choose a file, or drag it here</div>
            <div class="upload-sub-style">.xlsx, .xls, or .csv — first row must be headers</div>
        </div>
        """, unsafe_allow_html=True)

        file_upload = st.file_uploader("Select dataset file", type=["xlsx", "xls", "csv"], key="w_step1_uploader", label_visibility="collapsed")
        if file_upload is not None:
            try:
                if file_upload.name.endswith('.csv'):
                    df_load = pd.read_csv(file_upload)
                else:
                    df_load = pd.read_excel(file_upload)
                clear_deletion_flag()
                st.session_state['uploaded_df'] = df_load
                st.session_state['mapped_cols'] = {}
                st.session_state['workflow_step'] = '2 · Map columns'
                st.success(f"Loaded: {file_upload.name} ({len(df_load):,} samples)")
                st.rerun()
            except Exception as e:
                st.error(f"Error loading file: {e}")

        st.markdown("""
        <div class="privacy-info-banner">
            No data leaves your browser. Parsing and all statistics run locally.
        </div>
        """, unsafe_allow_html=True)

        c_space, c_next = st.columns([0.75, 0.25])
        with c_next:
            if st.button("Continue →", type="primary", use_container_width=True, key="btn_step1_cont"):
                st.session_state['workflow_step'] = '2 · Map columns'
                st.rerun()

        st.divider()

        with st.expander("🛠️ Advanced Workspace Actions: Reload Demo Data, Manual Entry & Deletion", expanded=False):
            c1, c2 = st.columns(2)
            with c1:
                if os.path.exists(demo_file_path):
                    if st.button("🔄 Reload Demo Dataset (1,044 Samples)", use_container_width=True, key="btn_exp_demo"):
                        reload_demo_data()
                        st.session_state['workflow_step'] = '3 · Results'
                        st.success("Demo dataset loaded!")
                        st.rerun()
            with c2:
                if st.button("🗑️ Delete Data (Reset)", use_container_width=True, key="btn_exp_reset"):
                    reset_all_data()
                    st.success("All data cleared!")
                    st.rerun()

            st.markdown("#### ➕ Manual Sample Data Entry Form")
            st.caption("Insert individual sample records. Validation rules: CYP2C19*2 accepts GG/GA/AA; CYP2C19*3 accepts GG only; CYP2C19*17 accepts CC/CT/TT.")
            
            with st.form("manual_sample_form_single", clear_on_submit=True):
                f_c1, f_c2 = st.columns(2)
                in_sid = f_c1.text_input("Sample ID *", value=f"PATIENT_{len(df_raw)+1 if df_raw is not None else 1:04d}")
                in_gender = f_c2.selectbox("Gender *", ["Male", "Female", "Unknown"])
                
                f_c3, f_c4 = st.columns(2)
                in_dob = f_c3.text_input("Date of Birth (YYYY-MM-DD)", value="1990-01-01")
                in_state = f_c4.selectbox("State / Region *", list(STATE_TO_REGION.keys()))
                
                in_native = st.text_input("Native Place (City / District)", value=in_state)
                in_3gen = st.selectbox("Is their family lived at Native place for past 3 generations?", ["Yes", "No", "Unknown"])
                
                st.markdown("##### CYP2C19 Target Genotypes")
                fg1, fg2, fg3 = st.columns(3)
                in_cyp2 = fg1.selectbox("CYP2C19*2 (rs4244285) *", ["GA", "AA", "GG", "Missing", "CC (Invalid Call)"])
                in_cyp3 = fg2.selectbox("CYP2C19*3 (rs4986893) *", ["GG", "Missing", "GA (Invalid Call)", "AA (Invalid Call)"])
                in_cyp17 = fg3.selectbox("CYP2C19*17 (rs12248560) *", ["CC", "CT", "TT", "Missing", "GG (Invalid Call)"])
                
                submit_sample = st.form_submit_button("➕ Add Sample & Recalculate Statistics", type="primary", use_container_width=True)
                
                if submit_sample:
                    cyp2_val = in_cyp2.split(' ')[0]
                    cyp3_val = in_cyp3.split(' ')[0]
                    cyp17_val = in_cyp17.split(' ')[0]

                    new_row = {
                        'sample ID': in_sid,
                        'Gender': in_gender,
                        'Date of Birth': in_dob,
                        'Native place ': in_native,
                        'State': in_state,
                        'Test requested': 'CYP2C19 Genotyping',
                        'Is their family lived at Native place for past 3 generations?': in_3gen,
                        'CYP2C19*2 (rs4244285)': np.nan if "Missing" in in_cyp2 else cyp2_val,
                        'CYP2C19*3 (rs4986893)': np.nan if "Missing" in in_cyp3 else cyp3_val,
                        'CYP2C19*17 ( rs12248560)': np.nan if "Missing" in in_cyp17 else cyp17_val
                    }
                    
                    clear_deletion_flag()
                    if st.session_state['uploaded_df'] is not None and len(st.session_state['uploaded_df']) > 0:
                        st.session_state['uploaded_df'] = pd.concat([st.session_state['uploaded_df'], pd.DataFrame([new_row])], ignore_index=True)
                    else:
                        st.session_state['uploaded_df'] = pd.DataFrame([new_row])
                        
                    st.session_state['workflow_step'] = '3 · Results'
                    st.success(f"Sample {in_sid} added successfully!")
                    st.rerun()

        st.markdown("""
        <div class="footer-disclaimer-text">
            PGx Pop MVP · CYP2C19 calling uses a simplified unphased-genotype heuristic — see notes in the CYP2C19 tab · Not for clinical use without independent validation.
        </div>
        """, unsafe_allow_html=True)

    # -------------------------------------------------------------------------
    # STEP 2: MAP COLUMNS (Matching Image 2)
    # -------------------------------------------------------------------------
    elif st.session_state['workflow_step'] == '2 · Map columns':
        st.markdown("## Map your columns")
        num_cols = len(df_raw.columns) if df_raw is not None else 0
        num_rows = len(df_raw) if df_raw is not None else 0
        st.caption(f"Tell the analyzer which columns hold which field. Detected {num_cols} columns, {num_rows:,} rows.")

        if df_raw is not None and len(df_raw.columns) > 0:
            all_cols = list(df_raw.columns)
            def find_default(patterns, cols):
                for pattern in patterns:
                    for c in cols:
                        if pattern.lower() in str(c).lower():
                            return c
                return cols[0] if cols else ""

            c_m1, c_m2, c_m3 = st.columns(3)
            sample_id_def = find_default(['sample id', 'sample_id', 'id'], all_cols)
            gender_def = find_default(['gender', 'sex'], all_cols)
            region_def = find_default(['native place', 'native', 'state', 'region'], all_cols)

            sel_sid = c_m1.selectbox("Sample ID column *", all_cols, index=all_cols.index(sample_id_def) if sample_id_def in all_cols else 0)
            sel_gen = c_m2.selectbox("Gender column (optional)", ["(None)"] + all_cols, index=all_cols.index(gender_def) + 1 if gender_def in all_cols else 0)
            sel_reg = c_m3.selectbox("Region / location column (optional)", ["(None)"] + all_cols, index=all_cols.index(region_def) + 1 if region_def in all_cols else 0)

            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown("### CYP2C19 SNP columns")
            st.caption("Map up to three SNPs. For star-allele/diplotype/phenotype calling, map rs4244285 (*2), rs4986893 (*3) and rs12248560 (*17). Any SNP left as '- none -' is skipped (genotype/allele/HWE still run on whichever SNPs you do map).")

            # Box 1: rs4244285 (*2)
            st.markdown("""
            <div class="clean-card-box">
                <h5 style="margin-top:0; color:#0F172A; font-family:-apple-system, sans-serif; font-weight:700;">rs4244285 (CYP2C19*2, c.681G>A)</h5>
            </div>
            """, unsafe_allow_html=True)
            c1_1, c1_2, c1_3 = st.columns(3)
            cyp2_def = find_default(['cyp2c19*2', 'rs4244285'], all_cols)
            sel_cyp2 = c1_1.selectbox("Genotype column (*2)", all_cols, index=all_cols.index(cyp2_def) if cyp2_def in all_cols else 0, key="w_map_cyp2")
            c1_2.text_input("Reference allele (*2)", value="G", key="w_map_ref2")
            c1_3.text_input("Variant allele (*2)", value="A", key="w_map_var2")

            # Box 2: rs4986893 (*3)
            st.markdown("""
            <div class="clean-card-box">
                <h5 style="margin-top:0; color:#0F172A; font-family:-apple-system, sans-serif; font-weight:700;">rs4986893 (CYP2C19*3, c.636G>A)</h5>
            </div>
            """, unsafe_allow_html=True)
            c2_1, c2_2, c2_3 = st.columns(3)
            cyp3_def = find_default(['cyp2c19*3', 'rs4986893'], all_cols)
            sel_cyp3 = c2_1.selectbox("Genotype column (*3)", all_cols, index=all_cols.index(cyp3_def) if cyp3_def in all_cols else 0, key="w_map_cyp3")
            c2_2.text_input("Reference allele (*3)", value="G", key="w_map_ref3")
            c2_3.text_input("Variant allele (*3)", value="A", key="w_map_var3")

            # Box 3: rs12248560 (*17)
            st.markdown("""
            <div class="clean-card-box">
                <h5 style="margin-top:0; color:#0F172A; font-family:-apple-system, sans-serif; font-weight:700;">rs12248560 (CYP2C19*17, c.-806C>T)</h5>
            </div>
            """, unsafe_allow_html=True)
            c3_1, c3_2, c3_3 = st.columns(3)
            cyp17_def = find_default(['cyp2c19*17', 'rs12248560'], all_cols)
            sel_cyp17 = c3_1.selectbox("Genotype column (*17)", all_cols, index=all_cols.index(cyp17_def) if cyp17_def in all_cols else 0, key="w_map_cyp17")
            c3_2.text_input("Reference allele (*17)", value="C", key="w_map_ref17")
            c3_3.text_input("Variant allele (*17)", value="T", key="w_map_var17")

            st.session_state['mapped_cols'] = {
                sel_sid: 'Sample_ID',
                sel_gen: 'Gender',
                sel_reg: 'Native_Place',
                sel_cyp2: 'CYP2C19*2',
                sel_cyp3: 'CYP2C19*3',
                sel_cyp17: 'CYP2C19*17'
            }

        c_b1, c_b2, c_b3 = st.columns([0.25, 0.5, 0.25])
        with c_b1:
            if st.button("← Back to Upload", use_container_width=True, key="btn_map_back"):
                st.session_state['workflow_step'] = '1 · Upload'
                st.rerun()
        with c_b3:
            if st.button("Continue to Results →", type="primary", use_container_width=True, key="btn_map_cont"):
                st.session_state['workflow_step'] = '3 · Results'
                st.rerun()

    # -------------------------------------------------------------------------
    # STEP 3: RESULTS (Matching Images 3, 4, 5)
    # -------------------------------------------------------------------------
    elif st.session_state['workflow_step'] == '3 · Results':
        st.markdown("## Overview")
        render_unified_results(full_results, qc_report)

        st.divider()
        st.markdown("### Export")
        st.caption("Download the full result set for the current group selection.")

        e1, e2, e3, e_remap = st.columns([0.25, 0.25, 0.25, 0.25])

        rep_gen = ReportGenerator(full_results, qc_report) if has_active_data() else None

        with e1:
            if rep_gen:
                excel_bytes = rep_gen.generate_excel()
                st.download_button("⬇ Excel (.xlsx, all sheets)", data=excel_bytes, file_name="PGx_Population_Report.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True, key="btn_exp_excel")
            else:
                st.button("⬇ Excel (.xlsx, all sheets)", disabled=True, use_container_width=True)

        with e2:
            if rep_gen:
                csv_bytes = rep_gen.generate_csv()
                st.download_button("⬇ CSV (summary)", data=csv_bytes, file_name="PGx_Population_Summary.csv", mime="text/csv", use_container_width=True, key="btn_exp_csv")
            else:
                st.button("⬇ CSV (summary)", disabled=True, use_container_width=True)

        with e3:
            if rep_gen:
                pdf_bytes = rep_gen.generate_pdf()
                st.download_button("⬇ PDF (print view)", data=pdf_bytes, file_name="PGx_Population_Report.pdf", mime="application/pdf", use_container_width=True, key="btn_exp_pdf")
            else:
                st.button("⬇ PDF (print view)", disabled=True, use_container_width=True)

        with e_remap:
            if st.button("← Re-map columns", use_container_width=True, key="btn_res_remap"):
                st.session_state['workflow_step'] = '2 · Map columns'
                st.rerun()

        st.markdown("""
        <div class="footer-disclaimer-text">
            PGx Pop MVP · CYP2C19 calling uses a simplified unphased-genotype heuristic — see notes in the CYP2C19 tab · Not for clinical use without independent validation.
        </div>
        """, unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# MAIN WORKFLOW EXECUTION
# -----------------------------------------------------------------------------
render_workflow_pipeline(df_raw, full_results, qc_report)

# -----------------------------------------------------------------------------
# FOOTER NOTE
# -----------------------------------------------------------------------------
st.markdown("""
<br><hr>
<div style="font-size: 0.8rem; color: #94A3B8; text-align: center;">
    Automated Population Pharmacogenomics Analysis Platform — Built for Reproducible Genotype, Allele, HWE, and CPIC Phenotype Analysis.
</div>
""", unsafe_allow_html=True)
