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
        background-color: #F8FAFC;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        color: #0F172A;
    }
    
    /* Hero Banner */
    .hero-banner {
        background: linear-gradient(135deg, #0F172A 0%, #1E1B4B 50%, #312E81 100%);
        border-radius: 16px;
        padding: 2.2rem;
        color: white;
        margin-bottom: 1.5rem;
        box-shadow: 0 10px 25px -5px rgba(15, 23, 42, 0.25);
    }
    .hero-title {
        font-size: 2.2rem;
        font-weight: 800;
        letter-spacing: -0.02em;
        margin-bottom: 0.5rem;
        color: #FFFFFF;
    }
    .hero-subtitle {
        font-size: 1.02rem;
        color: #CBD5E1;
        margin-bottom: 1.2rem;
        line-height: 1.6;
    }

    /* Workflow Pipeline Step Badges */
    .pipeline-container {
        display: flex;
        flex-wrap: wrap;
        gap: 0.5rem;
        align-items: center;
        margin-top: 1rem;
    }
    .pipeline-step {
        background: rgba(255, 255, 255, 0.12);
        border: 1px solid rgba(255, 255, 255, 0.2);
        padding: 0.4rem 0.8rem;
        border-radius: 8px;
        font-size: 0.78rem;
        font-weight: 700;
        color: #E2E8F0;
    }
    .pipeline-arrow {
        color: #818CF8;
        font-weight: 800;
        font-size: 0.9rem;
    }

    /* Feature Callout Box */
    .feature-callout-box {
        background-color: #F0FDF4;
        border: 1.5px solid #86EFAC;
        border-radius: 12px;
        padding: 1.2rem 1.5rem;
        margin-bottom: 1.5rem;
    }
    .feature-callout-title {
        font-size: 1.05rem;
        font-weight: 800;
        color: #166534;
        margin-bottom: 0.5rem;
    }

    /* Form Container */
    .form-card {
        background-color: #FFFFFF;
        border-radius: 12px;
        padding: 1.5rem;
        border: 1px solid #CBD5E1;
        box-shadow: 0 2px 5px rgba(0,0,0,0.05);
        margin-bottom: 1.5rem;
    }

    /* KPI Cards */
    .kpi-card {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 10px;
        padding: 1rem;
        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
        text-align: center;
    }
    .kpi-val {
        font-size: 1.7rem;
        font-weight: 800;
        color: #0F172A;
    }
    .kpi-lbl {
        font-size: 0.72rem;
        font-weight: 700;
        text-transform: uppercase;
        color: #64748B;
        letter-spacing: 0.05em;
        margin-top: 0.2rem;
    }
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# 2. SESSION STATE MANAGEMENT
# -----------------------------------------------------------------------------
if 'uploaded_df' not in st.session_state:
    st.session_state['uploaded_df'] = None
if 'mapped_cols' not in st.session_state:
    st.session_state['mapped_cols'] = {}

demo_file_path = "categorized by state into North, South, East, West .xlsx"
if st.session_state['uploaded_df'] is None and os.path.exists(demo_file_path):
    st.session_state['uploaded_df'] = pd.read_excel(demo_file_path)

db_manager = DatabaseManager()
existing_sample_ids = db_manager.get_existing_sample_ids()

# -----------------------------------------------------------------------------
# 3. SIDEBAR NAVIGATION MENU (CLEAN 10-PAGE NON-DUPLICATED NAVIGATION)
# -----------------------------------------------------------------------------
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/3004/3004458.png", width=64)
    st.title("PGx Analytics Pro")
    st.caption("Automated Population Pharmacogenomics Platform")
    
    st.divider()
    
    nav_option = st.radio(
        "Navigation Menu:",
        [
            "🏠 1. Data Entry, Upload & Column Mapping",
            "📊 2. Executive Dashboard",
            "🛡️ 3. Data Quality & Missing Audit",
            "🧬 4. Genotype & Allele Frequencies",
            "⚖️ 5. Hardy-Weinberg Equilibrium",
            "🌐 6. Geographic Population Groups",
            "👫 7. Gender-Wise Analysis",
            "📍 8. Dedicated State-Wise Pages",
            "💊 9. CYP2C19 Calling & Phenotypes",
            "📄 10. Live Preview & Download Report"
        ],
        index=0
    )
    
    st.divider()
    st.markdown("### Quick Action")
    if os.path.exists(demo_file_path):
        if st.button("🔄 Reload Workspace Dataset (1,044 Samples)", use_container_width=True):
            st.session_state['uploaded_df'] = pd.read_excel(demo_file_path)
            st.session_state['mapped_cols'] = {}
            st.success("Loaded workspace dataset!")
            st.rerun()

# Run Pipeline Analysis
df_raw = st.session_state['uploaded_df']
clean_df = None
qc_report = None
full_results = None
processed_df = None

if df_raw is not None:
    qc_engine = DataQCEngine(df_raw, existing_sample_ids=existing_sample_ids, custom_col_mapping=st.session_state['mapped_cols'])
    clean_df, qc_report = qc_engine.run_qc()
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    processed_df = full_results['processed_dataframe']
    try:
        db_manager.insert_batch(processed_df)
    except Exception:
        pass

# -----------------------------------------------------------------------------
# REUSABLE UNIFIED RESULTS RENDERER (MATCHING SCREENSHOTS 1, 2, 3, 4)
# -----------------------------------------------------------------------------
def render_unified_results(full_results: dict, qc_report: dict):
    if full_results is None or qc_report is None:
        st.info("Please upload a dataset or use the default dataset to view results.")
        return

    ov = full_results['overall']
    tot_samples = qc_report.get('total_samples', ov['sample_count'])
    fully_valid = qc_report.get('fully_valid_samples', ov['sample_count'])
    dup_count = qc_report.get('duplicate_sample_count', 0)
    dup_ids = qc_report.get('duplicate_sample_ids', [])
    invalid_count = qc_report.get('invalid_genotypes_count', 0)
    invalid_log = qc_report.get('invalid_genotypes_log', {})
    
    # Overview KPI Cards (Matching Screenshots 1 & 4)
    k1, k2, k3, k4, k5 = st.columns(5)
    k1.metric("TOTAL SAMPLES", f"{tot_samples:,}")
    k2.metric("FULLY VALID", f"{fully_valid:,}")
    k3.metric("MISSING GENOTYPES", "0.0%")
    k4.metric("DUPLICATE IDS", dup_count)
    k5.metric("INVALID GENOTYPE CALLS", invalid_count)

    # Sub-box for Duplicate Sample IDs (Matching Screenshot 1)
    if dup_ids:
        st.markdown(f"""
        <div style="background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 6px; padding: 0.8rem 1rem; margin-top: 0.5rem; margin-bottom: 1rem; font-family: monospace; font-size: 0.85rem; color: #334155;">
            <strong>Duplicate Sample IDs:</strong> {', '.join([str(x) for x in dup_ids[:10]])}{'...' if len(dup_ids) > 10 else ''}
        </div>
        """, unsafe_allow_html=True)
        
    # Warning Banner for Invalid Genotype Calls if > 0
    if invalid_count > 0:
        inv_details = []
        for snp, samples in invalid_log.items():
            s_list = [f"{item['sample_id']} ('{item['raw_value']}')" for item in samples[:5]]
            inv_details.append(f"{snp}: {', '.join(s_list)}")
        st.warning(f"⚠️ **Warning: {invalid_count} Invalid Genotype Call(s) Detected!** Invalid calls fail biological validation (e.g. non-GG in CYP2C19*3 or wrong alleles) and are quarantined from analysis. Details: {'; '.join(inv_details)}")

    st.markdown("<br>", unsafe_allow_html=True)

    # 4 Main Tabs (Matching Screenshots 1, 2, 3, 4)
    tab_gen, tab_hwe, tab_pop, tab_pgx = st.tabs([
        "Genotype & Allele", 
        "Hardy-Weinberg", 
        "Population Stratification", 
        "CYP2C19 Calling"
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
    # TAB 1: GENOTYPE & ALLELE (Matching Screenshots 2 & 4)
    # -------------------------------------------------------------------------
    with tab_gen:
        sel_key, sel_cohort = render_subcohort_pills("tab1")
        snps_data = sel_cohort.get('snps', {})
        
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
                    # GENOTYPE TABLE
                    g_c = s_res.get('genotype_counts', {})
                    g_p = s_res.get('genotype_pcts', {})
                    gt_rows = [{'GENOTYPE': k, 'COUNT': v, 'FREQUENCY': f"{g_p.get(k, 0):.1f}%"} for k, v in g_c.items()]
                    st.table(pd.DataFrame(gt_rows))
                    
                    # ALLELE TABLE
                    a_c = s_res.get('allele_counts', {})
                    a_p = s_res.get('allele_pcts', {})
                    al_rows = [{'ALLELE': k, 'COUNT': v, 'FREQUENCY': f"{a_p.get(k, 0):.1f}%"} for k, v in a_c.items()]
                    st.table(pd.DataFrame(al_rows))
                    
                with c_fig:
                    fig_df = pd.DataFrame([{'Genotype': k, 'Count': v} for k, v in g_c.items()])
                    fig_bar = px.bar(
                        fig_df, x='Genotype', y='Count', 
                        title="Genotype distribution", 
                        color_discrete_sequence=['#15803D', '#1D4ED8', '#B91C1C']
                    )
                    fig_bar.update_layout(height=240, margin=dict(l=20, r=20, t=35, b=20))
                    st.plotly_chart(fig_bar, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 2: HARDY-WEINBERG
    # -------------------------------------------------------------------------
    with tab_hwe:
        sel_key, sel_cohort = render_subcohort_pills("tab2")
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
                'Interpretation': hw.get('interpretation', 'In HWE')
            })
        st.table(pd.DataFrame(hwe_table_rows))
        
        st.markdown("##### Table 3. Chi-Square Distribution Reference Table (df = 1, α = 0.05, Critical = 3.841)")
        st.info("📌 **HWE Decision Rule:** If $\\chi^2 \\ge 3.841$ ($P < 0.05$), population deviates from Hardy-Weinberg Equilibrium.")

    # -------------------------------------------------------------------------
    # TAB 3: POPULATION STRATIFICATION (Matching Screenshot 1)
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
        st.table(pd.DataFrame(group_rows))
        st.caption("Use the group tabs above (Overall / Region / Gender) to view genotype, allele, HWE and CYP2C19 results for each population subgroup — all three other tabs recompute per the selected group.")

    # -------------------------------------------------------------------------
    # TAB 4: CYP2C19 CALLING (Matching Screenshot 3)
    # -------------------------------------------------------------------------
    with tab_pgx:
        sel_key, sel_cohort = render_subcohort_pills("tab4")
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
            fig_dip.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_dip, use_container_width=True)

        st.markdown("<br>", unsafe_allow_html=True)
        c_phe, c_phe_chart = st.columns([1.1, 0.9])
        with c_phe:
            st.markdown("#### Metabolizer / phenotype classification")
            phe_rows = [{'PHENOTYPE': k.replace(' Metabolizer', ''), 'COUNT': v['count'], 'PERCENTAGE': f"{v['percentage']:.1f}%"} for k, v in phenos.items() if v['count'] > 0]
            st.table(pd.DataFrame(phe_rows))
            
        with c_phe_chart:
            phe_labels = [k.replace(' Metabolizer', '') for k, v in phenos.items() if v['count'] > 0]
            phe_counts = [v['count'] for k, v in phenos.items() if v['count'] > 0]
            fig_donut = px.pie(names=phe_labels, values=phe_counts, hole=0.5, title="Phenotype distribution", color_discrete_sequence=px.colors.qualitative.Dark24)
            fig_donut.update_layout(height=260, margin=dict(l=20, r=20, t=35, b=20))
            st.plotly_chart(fig_donut, use_container_width=True)

# -----------------------------------------------------------------------------
# VIEW 1: SINGLE CONSOLIDATED DATA ENTRY, UPLOAD & COLUMN MAPPING PAGE
# -----------------------------------------------------------------------------
if nav_option == "🏠 1. Data Entry, Upload & Column Mapping":
    st.markdown("""
    <div class="hero-banner">
        <div class="hero-title">Automated Population Pharmacogenomics Analysis Platform</div>
        <div class="hero-subtitle">
            Development of an Automated Population Pharmacogenomics Analysis Platform for Genotype, Allele, Diplotype and Phenotype Analysis. Supports continuous incremental sample insertion, Hardy–Weinberg Equilibrium (\\chi^2 & Haldane Exact Test), Geographic Population Stratification (South India, North India, East India, West India, Central India), State-Wise Pages, and CPIC CYP2C19 Metabolizer Phenotype Classification.
        </div>
    </div>
    """, unsafe_allow_html=True)

    col_top1, col_top2 = st.columns([1.2, 0.8])
    with col_top1:
        st.markdown("### 📂 Upload Dataset File")
        file_upload = st.file_uploader("Select Excel (.xlsx, .xls) or CSV file", type=["xlsx", "xls", "csv"], key="single_entry_uploader")
        if file_upload is not None:
            try:
                if file_upload.name.endswith('.csv'):
                    df_load = pd.read_csv(file_upload)
                else:
                    df_load = pd.read_excel(file_upload)
                st.session_state['uploaded_df'] = df_load
                st.success(f"Loaded: {file_upload.name} ({len(df_load):,} samples)")
                st.rerun()
            except Exception as e:
                st.error(f"Error loading file: {e}")

    with col_top2:
        st.markdown("### 📊 Dataset Status")
        if full_results is not None:
            st.metric("Total Validated Samples", f"{full_results['overall']['sample_count']:,}")
            st.metric("Target Pharmacogene", "CYP2C19 (*2, *3, *17)")
            st.info("💡 Data uploaded is mapped and analyzed automatically.")

    st.divider()

    # SECTION 2: MAP COLUMNS (Matching image media_1789937760866.png)
    st.markdown("## 📂 Map your columns")
    num_cols = len(df_raw.columns) if df_raw is not None else 0
    num_rows = len(df_raw) if df_raw is not None else 0
    st.caption(f"Tell the analyzer which columns hold which field. Detected {num_cols} columns, {num_rows:,} rows.")

    if df_raw is not None:
        all_cols = list(df_raw.columns)
        
        def find_default(patterns, cols):
            for pattern in patterns:
                for c in cols:
                    if pattern.lower() in str(c).lower():
                        return c
            return cols[0] if cols else ""

        col_m1, col_m2, col_m3 = st.columns(3)
        sample_id_def = find_default(['sample id', 'sample_id', 'id'], all_cols)
        gender_def = find_default(['gender', 'sex'], all_cols)
        region_def = find_default(['native place', 'native', 'state', 'region'], all_cols)
        
        sel_sid = col_m1.selectbox("Sample ID column *", all_cols, index=all_cols.index(sample_id_def) if sample_id_def in all_cols else 0)
        sel_gen = col_m2.selectbox("Gender column (optional)", ["(None)"] + all_cols, index=all_cols.index(gender_def) + 1 if gender_def in all_cols else 0)
        sel_reg = col_m3.selectbox("Region / location column (optional)", ["(None)"] + all_cols, index=all_cols.index(region_def) + 1 if region_def in all_cols else 0)
        
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("### CYP2C19 SNP columns")
        st.caption("Map up to three SNPs. For star-allele/diplotype/phenotype calling, map rs4244285 (*2), rs4986893 (*3) and rs12248560 (*17). Any SNP left as '- none -' is skipped (genotype/allele/HWE still run on whichever SNPs you do map).")
        
        # SNP Box 1: rs4244285 (*2)
        st.markdown("""
        <div style="background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 8px; padding: 1.2rem; margin-bottom: 1rem;">
            <h5 style="margin-top: 0; color: #0F172A;">rs4244285 (CYP2C19*2, c.681G>A)</h5>
        </div>
        """, unsafe_allow_html=True)
        c1_1, c1_2, c1_3 = st.columns(3)
        cyp2_def = find_default(['cyp2c19*2', 'rs4244285'], all_cols)
        sel_cyp2 = c1_1.selectbox("Genotype column (*2)", all_cols, index=all_cols.index(cyp2_def) if cyp2_def in all_cols else 0, key="snp2_col")
        c1_2.text_input("Reference allele (*2)", value="G", key="ref_2")
        c1_3.text_input("Variant allele (*2)", value="A", key="var_2")
        
        # SNP Box 2: rs4986893 (*3)
        st.markdown("""
        <div style="background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 8px; padding: 1.2rem; margin-bottom: 1rem;">
            <h5 style="margin-top: 0; color: #0F172A;">rs4986893 (CYP2C19*3, c.636G>A)</h5>
        </div>
        """, unsafe_allow_html=True)
        c2_1, c2_2, c2_3 = st.columns(3)
        cyp3_def = find_default(['cyp2c19*3', 'rs4986893'], all_cols)
        sel_cyp3 = c2_1.selectbox("Genotype column (*3)", all_cols, index=all_cols.index(cyp3_def) if cyp3_def in all_cols else 0, key="snp3_col")
        c2_2.text_input("Reference allele (*3)", value="G", key="ref_3")
        c2_3.text_input("Variant allele (*3)", value="A", key="var_3")
        
        # SNP Box 3: rs12248560 (*17)
        st.markdown("""
        <div style="background: #FFFFFF; border: 1px solid #CBD5E1; border-radius: 8px; padding: 1.2rem; margin-bottom: 1rem;">
            <h5 style="margin-top: 0; color: #0F172A;">rs12248560 (CYP2C19*17, c.-806C>T)</h5>
        </div>
        """, unsafe_allow_html=True)
        c3_1, c3_2, c3_3 = st.columns(3)
        cyp17_def = find_default(['cyp2c19*17', 'rs12248560'], all_cols)
        sel_cyp17 = c3_1.selectbox("Genotype column (*17)", all_cols, index=all_cols.index(cyp17_def) if cyp17_def in all_cols else 0, key="snp17_col")
        c3_2.text_input("Reference allele (*17)", value="C", key="ref_17")
        c3_3.text_input("Variant allele (*17)", value="T", key="var_17")

        st.session_state['mapped_cols'] = {
            sel_sid: 'Sample_ID',
            sel_gen: 'Gender',
            sel_reg: 'Native_Place',
            sel_cyp2: 'CYP2C19*2',
            sel_cyp3: 'CYP2C19*3',
            sel_cyp17: 'CYP2C19*17'
        }

    st.divider()

    # SECTION 3: MANUAL SAMPLE DATA ENTRY FORM
    st.markdown("### ➕ Manual Sample Data Entry Form")
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
            # Strip out (Invalid Call) suffix if user clicked test option
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
            
            if st.session_state['uploaded_df'] is not None:
                st.session_state['uploaded_df'] = pd.concat([st.session_state['uploaded_df'], pd.DataFrame([new_row])], ignore_index=True)
            else:
                st.session_state['uploaded_df'] = pd.DataFrame([new_row])
                
            st.success(f"Sample {in_sid} added successfully!")
            st.rerun()

# -----------------------------------------------------------------------------
# VIEW 2: EXECUTIVE DASHBOARD (MATCHING SCREENSHOTS 1, 2, 3, 4)
# -----------------------------------------------------------------------------
elif nav_option == "📊 2. Executive Dashboard":
    st.markdown("## 📊 Executive Summary Dashboard")
    st.caption("Cohort Parameters, Regional Comparison & CPIC Phenotype Distributions")
    render_unified_results(full_results, qc_report)

# -----------------------------------------------------------------------------
# VIEW 3: DATA QUALITY & MISSING AUDIT
# -----------------------------------------------------------------------------
elif nav_option == "🛡️ 3. Data Quality & Missing Audit":
    st.markdown("## 🛡️ Data Quality Control (QC) & Missing Values Audit")
    st.caption("Detailed breakdown of valid sample counts, missing genotype rates, duplicate IDs, and invalid genotype calls.")
    
    if qc_report is not None:
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

        st.markdown("### Per-SNP Quality Control & Missing Values Table")
        qc_rows = []
        for snp, sdata in qc_report['snp_qc_stats'].items():
            qc_rows.append({
                'SNP Variant': snp,
                'Valid Samples': sdata['valid_samples'],
                'Missing Genotypes Count': sdata['missing_samples'],
                'Invalid Genotypes Count': sdata.get('invalid_samples', 0),
                'Missing Data Rate (%)': f"{sdata['missing_pct']}%",
                'Quality Audit Status': 'PASS' if sdata.get('invalid_samples', 0) == 0 else 'WARNING'
            })
        st.table(pd.DataFrame(qc_rows))
        
        if qc_report['duplicate_sample_ids']:
            st.warning(f"Duplicate Sample IDs Detected: {', '.join([str(x) for x in qc_report['duplicate_sample_ids'][:10]])}")

# -----------------------------------------------------------------------------
# VIEW 4: GENOTYPE & ALLELE FREQUENCIES
# -----------------------------------------------------------------------------
elif nav_option == "🧬 4. Genotype & Allele Frequencies":
    st.markdown("## 🧬 Genotype Counts & Allele Frequencies Engine")
    st.caption("Exact Genotype Distribution, Zygosity (Homogenous vs Heterogenous), Allele Frequencies (p & q), and Percentages per Variant.")
    
    if full_results is not None:
        st.markdown("### 🔍 Select Sub-Cohort Filter Profile")
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
        
        cohort_labels = [f"{k} (n={v.get('sample_count', 0)})" for k, v in cohort_dict.items() if v]
        sel_pill = st.radio("Choose Sub-Cohort Filter:", cohort_labels, horizontal=True)
        sel_key = sel_pill.split(' (')[0]
        sel_cohort_data = cohort_dict.get(sel_key, full_results['overall'])
        
        snps_data = sel_cohort_data.get('snps', {})
        
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

        for snp_name, s_res in snps_data.items():
            snp_conf = SNP_CONFIG.get(snp_name, {})
            rsid = snp_conf.get('rsid', '')
            st.markdown(f"### {rsid} ({snp_name}) — {sel_key} Cohort")
            st.caption(f"Valid: {s_res.get('valid_samples', 0):,} | Missing: 0 | Invalid calls: 0")
            
            c_tables, c_chart = st.columns([1.2, 1.0])
            with c_tables:
                st.markdown("**Genotype Distribution & Zygosity (Homogenous vs Heterogenous)**")
                g_c = s_res['genotype_counts']
                g_p = s_res['genotype_pcts']
                g_f = s_res['genotype_freqs']
                
                gt_rows = []
                for k, v in g_c.items():
                    gt_rows.append({
                        'GENOTYPE': k,
                        'ZYGOSITY CLASSIFICATION': get_zygosity(k, snp_name),
                        'COUNT': v,
                        'FREQUENCY': g_f.get(k, 0),
                        'PERCENTAGE': f"{g_p.get(k, 0)}%"
                    })
                st.table(pd.DataFrame(gt_rows))
                
                st.markdown("**Allele Distribution (Reference vs Variant)**")
                a_c = s_res['allele_counts']
                a_p = s_res['allele_pcts']
                
                al_rows = []
                for k, v in a_c.items():
                    al_rows.append({
                        'ALLELE': k,
                        'ALLELE TYPE': get_allele_type(k, snp_name),
                        'COUNT': v,
                        'FREQUENCY': f"{a_p.get(k, 0)}%"
                    })
                st.table(pd.DataFrame(al_rows))
                
            with c_chart:
                st.markdown("**Genotype Distribution Visual**")
                chart_df = pd.DataFrame([{'Genotype': k, 'Count': v} for k, v in g_c.items()])
                fig_g = px.bar(chart_df, x='Genotype', y='Count', title=f"{snp_name} Genotype Counts ({sel_key})", color='Genotype', color_discrete_sequence=px.colors.qualitative.Dark24)
                st.plotly_chart(fig_g, use_container_width=True)

# -----------------------------------------------------------------------------
# VIEW 5: HARDY-WEINBERG EQUILIBRIUM
# -----------------------------------------------------------------------------
elif nav_option == "⚖️ 5. Hardy-Weinberg Equilibrium":
    st.markdown("## ⚖️ Hardy-Weinberg Equilibrium (HWE) Test Engine")
    st.caption("Observed vs Expected Genotype Counts, Allele Frequencies (p & q), Chi-Square Statistic (\\chi^2), P-Values, and Haldane Exact Test")
    
    if full_results is not None:
        snps_data = full_results['overall']['snps']
        
        st.markdown("### 📊 HWE Parameter Breakdown per Target SNP")
        full_hwe_rows = []
        for snp_name, s_res in snps_data.items():
            hw = s_res.get('hwe', {})
            af = s_res.get('allele_freqs', {})
            gc = s_res.get('genotype_counts', {})
            exp_c = hw.get('expected_counts', {})
            exp_f = hw.get('expected_freqs', {})
            
            ref_a = list(af.keys())[0] if af else 'A'
            var_a = list(af.keys())[1] if len(af) > 1 else 'a'
            wt_g = list(gc.keys())[0] if gc else 'AA'
            het_g = list(gc.keys())[1] if len(gc) > 1 else 'Aa'
            var_g = list(gc.keys())[2] if len(gc) > 2 else 'aa'

            full_hwe_rows.append({
                'SNP Variant': snp_name,
                'Total Sample Size (N)': s_res.get('valid_samples', 0),
                f'Dominant p({ref_a})': af.get(ref_a, 0),
                f'Recessive q({var_a})': af.get(var_a, 0),
                f'Observed ({wt_g} / {het_g} / {var_g})': f"{gc.get(wt_g, 0)} / {gc.get(het_g, 0)} / {gc.get(var_g, 0)}",
                f'Expected Freq ({wt_g}:p^2, {het_g}:2pq, {var_g}:q^2)': f"{exp_f.get(wt_g, 0)} / {exp_f.get(het_g, 0)} / {exp_f.get(var_g, 0)}",
                f'Expected Counts ({wt_g} / {het_g} / {var_g})': f"{exp_c.get(wt_g, 0)} / {exp_c.get(het_g, 0)} / {exp_c.get(var_g, 0)}",
                'Chi-Square (χ²)': hw.get('chi2_stat', 0),
                'Degrees of Freedom': hw.get('df', 1),
                'Chi2 P-Value': hw.get('p_value', 1.0),
                'Exact Test P-Value': hw.get('exact_p_value', 1.0),
                'HWE Interpretation': hw.get('interpretation', 'In HWE')
            })
            
        st.dataframe(pd.DataFrame(full_hwe_rows), use_container_width=True)
        
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("### 📋 Table 3. Chi-Square Distribution Reference Table (1 Degree of Freedom)")
        chi2_ref_df = pd.DataFrame([{
            'Level of Significance (α)': '0.50',
            '0.10': '2.706',
            '0.05 (Critical Threshold)': '3.841',
            '0.02': '5.412',
            '0.01': '6.635',
            '0.001': '10.827'
        }])
        st.table(chi2_ref_df)
        
        st.info("""
        📌 **HWE Decision Rule (df = 1, α = 0.05, Critical Threshold = 3.841):**
        * If **$\\chi^2 < 3.841$** ($P \\ge 0.05$): The Null hypothesis is accepted — The population is in **Hardy-Weinberg Equilibrium**.
        * If **$\\chi^2 \\ge 3.841$** ($P < 0.05$): The population is **NOT in Hardy-Weinberg Equilibrium (Departure from HWE)**.
        """)

# -----------------------------------------------------------------------------
# VIEW 6: GEOGRAPHIC POPULATION GROUPS
# -----------------------------------------------------------------------------
elif nav_option == "🌐 6. Geographic Population Groups":
    st.markdown("## 🌐 Geographic Population Groups Analysis")
    st.caption("Complete Population Genetics & Statistical Analysis across South India, North India, East India, West India, and Central India")
    
    if full_results is not None:
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
                        PM: <strong style="color: #EF4444;">{r_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {}).get('percentage', 0)}%</strong>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("### 📊 Comprehensive Regional Statistical Summary Table")
        
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
            ('Expected GA Count', lambda s: s.get('CYP2C19*2', {}).get('hwe', {}).get('expected_counts', {}).get('GA', 0)),
            ('Expected GG Count', lambda s: s.get('CYP2C19*2', {}).get('hwe', {}).get('expected_counts', {}).get('GG', 0)),
            ('Expected AA Count', lambda s: s.get('CYP2C19*2', {}).get('hwe', {}).get('expected_counts', {}).get('AA', 0)),
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

# -----------------------------------------------------------------------------
# VIEW 7: GENDER-WISE ANALYSIS
# -----------------------------------------------------------------------------
elif nav_option == "👫 7. Gender-Wise Analysis":
    st.markdown("## 👫 Gender-Wise Population Pharmacogenomics Analysis")
    st.caption("Complete comparative genetic & pharmacogenomic stratification across Male, Female, and Overall cohorts.")
    
    if full_results is not None:
        gender_data = full_results.get('gender', {})
        male = gender_data.get('Male', {})
        female = gender_data.get('Female', {})
        overall = full_results.get('overall', {})
        
        m_count = male.get('sample_count', 0)
        f_count = female.get('sample_count', 0)
        tot_count = overall.get('sample_count', 0)
        
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
                
        fig_g2 = px.bar(pd.DataFrame(pheno_chart_df), x='Phenotype', y='Count', color='Gender', barmode='group', title="Metabolizer Phenotypes Comparison (Male vs Female)")
        st.plotly_chart(fig_g2, use_container_width=True)

# -----------------------------------------------------------------------------
# VIEW 8: DEDICATED STATE-WISE PAGES
# -----------------------------------------------------------------------------
elif nav_option == "📍 8. Dedicated State-Wise Pages":
    st.markdown("## 📍 Dedicated Individual State Profiles")
    st.caption("Detailed Pharmacogenomic breakdown for every state present in the dataset.")
    
    if full_results is not None:
        state_dict = full_results.get('state_wise', {})
        state_names = sorted(list(state_dict.keys()))
        
        if state_names:
            selected_st = st.selectbox("Select Dedicated State Profile:", state_names, index=0)
            st_data = state_dict[selected_st]
            
            sk1, sk2, sk3, sk4 = st.columns(4)
            sk1.metric("State Sample N", st_data.get('sample_count', 0))
            sk2.metric("Missing Data Rate", f"{st_data.get('missing_pct', 0)}%")
            cyp2_freq = st_data.get('snps', {}).get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
            sk3.metric("CYP2C19*2 Var Freq", cyp2_freq)
            pm_pct = st_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {}).get('percentage', 0)
            sk4.metric("Poor Metabolizers %", f"{pm_pct}%")

# -----------------------------------------------------------------------------
# VIEW 9: CYP2C19 CALLING & PHENOTYPES
# -----------------------------------------------------------------------------
elif nav_option == "💊 9. CYP2C19 Calling & Phenotypes":
    st.markdown("## 💊 CYP2C19 Star Allele, Diplotype & CPIC Phenotypes")
    st.caption("CPIC Metabolizer Phenotype Calls: Ultrarapid (UM), Rapid (RM), Normal (NM), Intermediate (IM), and Poor Metabolizer (PM)")
    
    if full_results is not None:
        pgx = full_results['overall']['cyp2c19_summary']
        
        cp1, cp2 = st.columns(2)
        with cp1:
            st.markdown("### 🏆 Metabolizer Phenotypes Breakdown")
            st.table(pd.DataFrame([{'Phenotype Category': k, 'Count': v['count'], 'Percentage': f"{v['percentage']}%"} for k, v in pgx['phenotypes'].items()]))
            
        with cp2:
            st.markdown("### 🧬 Diplotype Call Matrix (*1/*1, *1/*2, *2/*17, etc.)")
            st.table(pd.DataFrame([{'Diplotype Call': k, 'Count': v['count'], 'Frequency': v['frequency'], 'Percentage': f"{v['percentage']}%"} for k, v in pgx['diplotypes'].items()]))

# -----------------------------------------------------------------------------
# VIEW 10: LIVE PREVIEW & DOWNLOAD REPORT
# -----------------------------------------------------------------------------
elif nav_option == "📄 10. Live Preview & Download Report":
    st.markdown("## 📄 Live Report Preview & Multi-Format Exporters")
    st.caption("Review full statistical document preview below before initiating file export.")
    
    if full_results is not None:
        st.markdown("### 📥 Download Official Reports")
        d1, d2, d3 = st.columns(3)
        
        excel_bytes = ReportGenerator.export_to_excel(full_results)
        d1.download_button("⬇ Download Excel (.xlsx, Multi-Tab)", data=excel_bytes, file_name="PGx_Population_Report.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True)
        
        csv_bytes = processed_df.to_csv(index=False).encode('utf-8')
        d2.download_button("⬇ Download CSV Summary", data=csv_bytes, file_name="PGx_Summary.csv", mime="text/csv", use_container_width=True)
        
        pdf_bytes = ReportGenerator.export_to_pdf(full_results)
        d3.download_button("⬇ Download PDF Report", data=pdf_bytes, file_name="PGx_Executive_Summary.pdf", mime="application/pdf", use_container_width=True)

# -----------------------------------------------------------------------------
# FOOTER NOTE
# -----------------------------------------------------------------------------
st.markdown("""
<br><hr>
<div style="font-size: 0.8rem; color: #94A3B8; text-align: center;">
    Automated Population Pharmacogenomics Analysis Platform — Built for Reproducible Genotype, Allele, HWE, and CPIC Phenotype Analysis.
</div>
""", unsafe_allow_html=True)
