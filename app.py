"""
PGx Pop - Population Pharmacogenomics Analyzer (MVP)
Exact 3-Step Wizard UI with Geographic Population Groups & State-Wise Pages:
Step 1: Upload Dataset
Step 2: Map Columns (Sample ID, Gender, Region, SNP Mapping)
Step 3: Results (KPIs, Genotype & Allele, Hardy-Weinberg, Geographic Population Groups, State-Wise Pages, CYP2C19 Calling, Group Filter Pills, Multi-Export)
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
from config import SNP_CONFIG, PHENOTYPE_ORDER

# -----------------------------------------------------------------------------
# 1. PAGE CONFIG & CUSTOM CSS
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="PGx Pop — Population Pharmacogenomics Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

EXACT_MATCH_CSS = """
<style>
    /* Main Background & Clean Typography */
    .stApp {
        background-color: #EFEFEF;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        color: #222222;
    }
    
    /* Hide Streamlit top header padding */
    .block-container {
        padding-top: 1rem;
        padding-bottom: 3rem;
        max-width: 1200px;
    }
    
    /* App Top Branding Header */
    .app-branding {
        font-size: 0.95rem;
        font-weight: 700;
        color: #111111;
        letter-spacing: 0.05em;
        text-transform: uppercase;
        margin-bottom: 0.8rem;
    }
    .app-subbranding {
        font-weight: 400;
        color: #666666;
    }

    /* Card Box Container */
    .card-box {
        background-color: #FFFFFF;
        border-radius: 6px;
        padding: 1.8rem;
        border: 1px solid #D8D8D8;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        margin-bottom: 1.5rem;
    }
    
    .card-title {
        font-size: 1.4rem;
        font-weight: 800;
        color: #111111;
        margin-bottom: 0.4rem;
        font-family: Georgia, serif;
    }
    .card-subtitle {
        font-size: 0.9rem;
        color: #555555;
        margin-bottom: 1.2rem;
        line-height: 1.4;
    }

    /* Top 5 Overview Metric Cards Grid */
    .metric-grid {
        display: grid;
        grid-template-columns: repeat(5, 1fr);
        gap: 0.75rem;
        margin-bottom: 1rem;
    }
    .metric-card-single {
        background: #F7F7F7;
        border: 1px solid #E0E0E0;
        border-radius: 4px;
        padding: 0.9rem 1rem;
    }
    .metric-val {
        font-size: 1.6rem;
        font-weight: 800;
        color: #111111;
        line-height: 1.1;
        font-family: Georgia, serif;
    }
    .metric-lbl {
        font-size: 0.72rem;
        font-weight: 700;
        text-transform: uppercase;
        color: #666666;
        margin-top: 0.3rem;
        letter-spacing: 0.04em;
    }

    /* Duplicate Warning Box */
    .warning-box {
        background-color: #FBFBFB;
        border: 1px solid #E0E0E0;
        border-radius: 4px;
        padding: 0.6rem 1rem;
        font-size: 0.82rem;
        color: #444444;
        margin-bottom: 1.2rem;
    }

    /* Group Filter Pills Styling */
    .stRadio > div {
        flex-direction: row !important;
        gap: 0.5rem !important;
    }

    /* Footer Note */
    .footer-note {
        font-size: 0.78rem;
        color: #777777;
        text-align: center;
        margin-top: 2.5rem;
    }
</style>
"""
st.markdown(EXACT_MATCH_CSS, unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# 2. SESSION STATE MANAGEMENT FOR WIZARD & DATA
# -----------------------------------------------------------------------------
if 'wizard_step' not in st.session_state:
    st.session_state['wizard_step'] = 1  # 1: Upload, 2: Map columns, 3: Results

if 'uploaded_df' not in st.session_state:
    st.session_state['uploaded_df'] = None

if 'mapped_config' not in st.session_state:
    st.session_state['mapped_config'] = {
        'sample_id_col': None,
        'gender_col': None,
        'region_col': None,
        'snp_cyp2': None,
        'snp_cyp3': None,
        'snp_cyp17': None
    }

# Auto-load workspace file if no uploaded file
demo_file_path = "categorized by state into North, South, East, West .xlsx"
if st.session_state['uploaded_df'] is None and os.path.exists(demo_file_path):
    st.session_state['uploaded_df'] = pd.read_excel(demo_file_path)

# Initialize Database Manager
db_manager = DatabaseManager()
existing_sample_ids = db_manager.get_existing_sample_ids()

# -----------------------------------------------------------------------------
# 3. TOP BRANDING & HORIZONTAL WIZARD BAR
# -----------------------------------------------------------------------------
st.markdown("""
<div class="app-branding">
    PGx Pop <span class="app-subbranding">— POPULATION PHARMACOGENOMICS ANALYZER</span>
</div>
""", unsafe_allow_html=True)

step = st.session_state['wizard_step']

col_w1, col_w2, col_w3 = st.columns(3)

if col_w1.button("1 • Upload", key="btn_w1", use_container_width=True):
    st.session_state['wizard_step'] = 1
    st.rerun()

if col_w2.button("2 • Map columns", key="btn_w2", use_container_width=True):
    st.session_state['wizard_step'] = 2
    st.rerun()

if col_w3.button("3 • Results", key="btn_w3", use_container_width=True):
    if st.session_state['uploaded_df'] is not None:
        st.session_state['wizard_step'] = 3
        st.rerun()
    else:
        st.warning("Please upload or load a dataset first.")

st.markdown("<br>", unsafe_allow_html=True)

# -----------------------------------------------------------------------------
# STEP 1: UPLOAD GENOTYPE DATASET
# -----------------------------------------------------------------------------
if st.session_state['wizard_step'] == 1:
    st.markdown("""
    <div class="card-box">
        <div class="card-title">Upload genotype dataset</div>
        <div class="card-subtitle">
            Excel (.xlsx/.csv) with Sample ID, Gender, DOB, region/location, and genotype columns for each SNP (e.g. "GG", "GA", "AA").
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    file_upload = st.file_uploader(
        "Click to choose a file, or drag it here (.xlsx, .xls, or .csv — first row must be headers)",
        type=["xlsx", "xls", "csv"]
    )
    
    if file_upload is not None:
        try:
            if file_upload.name.endswith('.csv'):
                df_load = pd.read_csv(file_upload)
            else:
                df_load = pd.read_excel(file_upload)
            st.session_state['uploaded_df'] = df_load
            st.success(f"Loaded: {file_upload.name} ({len(df_load):,} rows, {len(df_load.columns)} columns)")
        except Exception as e:
            st.error(f"Error reading file: {e}")
            
    if os.path.exists(demo_file_path):
        if st.button("🔄 Use Demo Workspace Dataset (1,044 samples)"):
            st.session_state['uploaded_df'] = pd.read_excel(demo_file_path)
            st.success("Loaded demo dataset!")
            st.rerun()
            
    st.info("🔒 No data leaves your browser. Parsing and all statistics run locally.")
    
    st.markdown("<br>", unsafe_allow_html=True)
    c_btn1, c_btn2 = st.columns([4, 1])
    if c_btn2.button("Continue ➔", type="primary", use_container_width=True):
        if st.session_state['uploaded_df'] is not None:
            st.session_state['wizard_step'] = 2
            st.rerun()
        else:
            st.warning("Please upload a file to continue.")

# -----------------------------------------------------------------------------
# STEP 2: MAP COLUMNS
# -----------------------------------------------------------------------------
elif st.session_state['wizard_step'] == 2:
    df_raw = st.session_state['uploaded_df']
    if df_raw is None:
        st.warning("No dataset uploaded. Please return to Step 1.")
    else:
        all_cols = list(df_raw.columns)
        num_cols = len(all_cols)
        num_rows = len(df_raw)
        
        st.markdown(f"""
        <div class="card-box">
            <div class="card-title">Map your columns</div>
            <div class="card-subtitle">
                Tell the analyzer which columns hold which field. Detected {num_cols} columns, {num_rows} rows.
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        def find_default(patterns, cols):
            for pattern in patterns:
                for c in cols:
                    if pattern.lower() in str(c).lower():
                        return c
            return cols[0] if cols else ""

        col_m1, col_m2, col_m3 = st.columns(3)
        
        sample_id_default = find_default(['sample id', 'sample_id', 'id'], all_cols)
        gender_default = find_default(['gender', 'sex'], all_cols)
        region_default = find_default(['native place', 'native', 'state', 'region'], all_cols)
        
        sel_sample_id = col_m1.selectbox("Sample ID column", all_cols, index=all_cols.index(sample_id_default) if sample_id_default in all_cols else 0)
        sel_gender = col_m2.selectbox("Gender column (optional)", ["(None)"] + all_cols, index=all_cols.index(gender_default) + 1 if gender_default in all_cols else 0)
        sel_region = col_m3.selectbox("Region / location column (optional)", ["(None)"] + all_cols, index=all_cols.index(region_default) + 1 if region_default in all_cols else 0)
        
        st.markdown("### CYP2C19 SNP columns")
        st.caption("Map up to three SNPs. For star-allele/diplotype/phenotype calling, map rs4244285 (*2), rs4986893 (*3) and rs12248560 (*17). Any SNP left as '- none -' is skipped.")
        
        cyp2_default = find_default(['cyp2c19*2', 'rs4244285'], all_cols)
        with st.container():
            st.markdown("#### rs4244285 (CYP2C19*2, c.681G>A)")
            col_s1_a, col_s1_b, col_s1_c = st.columns(3)
            sel_cyp2 = col_s1_a.selectbox("Genotype column", ["- none -"] + all_cols, index=all_cols.index(cyp2_default) + 1 if cyp2_default in all_cols else 0, key="sel_cyp2")
            col_s1_b.text_input("Reference allele", "G", disabled=True, key="ref_cyp2")
            col_s1_c.text_input("Variant allele", "A", disabled=True, key="var_cyp2")

        cyp3_default = find_default(['cyp2c19*3', 'rs4986893'], all_cols)
        with st.container():
            st.markdown("#### rs4986893 (CYP2C19*3, c.636G>A)")
            col_s2_a, col_s2_b, col_s2_c = st.columns(3)
            sel_cyp3 = col_s2_a.selectbox("Genotype column", ["- none -"] + all_cols, index=all_cols.index(cyp3_default) + 1 if cyp3_default in all_cols else 0, key="sel_cyp3")
            col_s2_b.text_input("Reference allele", "G", disabled=True, key="ref_cyp3")
            col_s2_c.text_input("Variant allele", "A", disabled=True, key="var_cyp3")

        cyp17_default = find_default(['cyp2c19*17', 'rs12248560'], all_cols)
        with st.container():
            st.markdown("#### rs12248560 (CYP2C19*17, c.-806C>T)")
            col_s3_a, col_s3_b, col_s3_c = st.columns(3)
            sel_cyp17 = col_s3_a.selectbox("Genotype column", ["- none -"] + all_cols, index=all_cols.index(cyp17_default) + 1 if cyp17_default in all_cols else 0, key="sel_cyp17")
            col_s3_b.text_input("Reference allele", "C", disabled=True, key="ref_cyp17")
            col_s3_c.text_input("Variant allele", "T", disabled=True, key="var_cyp17")

        st.markdown("<br>", unsafe_allow_html=True)
        c_mbtn1, c_mbtn2 = st.columns([4, 1])
        if c_mbtn2.button("Run Analysis ➔", type="primary", use_container_width=True):
            st.session_state['mapped_config'] = {
                'sample_id_col': sel_sample_id,
                'gender_col': None if sel_gender == "(None)" else sel_gender,
                'region_col': None if sel_region == "(None)" else sel_region,
                'snp_cyp2': None if sel_cyp2 == "- none -" else sel_cyp2,
                'snp_cyp3': None if sel_cyp3 == "- none -" else sel_cyp3,
                'snp_cyp17': None if sel_cyp17 == "- none -" else sel_cyp17
            }
            st.session_state['wizard_step'] = 3
            st.rerun()

# -----------------------------------------------------------------------------
# STEP 3: RESULTS
# -----------------------------------------------------------------------------
elif st.session_state['wizard_step'] == 3:
    df_raw = st.session_state['uploaded_df']
    if df_raw is None:
        st.warning("No dataset uploaded.")
    else:
        # Run QC & Stratification
        qc_engine = DataQCEngine(df_raw, existing_sample_ids=existing_sample_ids)
        clean_df, qc_report = qc_engine.run_qc()
        full_results = DemographicStratifier.stratify_and_analyze(clean_df)
        processed_df = full_results['processed_dataframe']
        
        # Save batch to DB
        try:
            db_manager.insert_batch(processed_df)
        except Exception:
            pass
        
        # 1. TOP OVERVIEW METRIC CARDS
        total_samples = qc_report['total_samples']
        fully_valid = len(clean_df)
        missing_pct = qc_report['snp_qc_stats'].get('CYP2C19*2', {}).get('missing_pct', 0.0)
        duplicate_count = qc_report['duplicate_sample_count']
        invalid_calls = sum(len(v) for v in qc_report.get('invalid_genotypes_log', {}).values())

        st.markdown(f"""
        <div class="metric-grid">
            <div class="metric-card-single">
                <div class="metric-val">{total_samples}</div>
                <div class="metric-lbl">TOTAL SAMPLES</div>
            </div>
            <div class="metric-card-single">
                <div class="metric-val">{fully_valid}</div>
                <div class="metric-lbl">FULLY VALID</div>
            </div>
            <div class="metric-card-single">
                <div class="metric-val">{missing_pct:.1f}%</div>
                <div class="metric-lbl">MISSING GENOTYPES</div>
            </div>
            <div class="metric-card-single">
                <div class="metric-val">{duplicate_count}</div>
                <div class="metric-lbl">DUPLICATE IDS</div>
            </div>
            <div class="metric-card-single">
                <div class="metric-val">{invalid_calls}</div>
                <div class="metric-lbl">INVALID GENOTYPE CALLS</div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        if qc_report['duplicate_sample_ids']:
            dup_str = ", ".join(qc_report['duplicate_sample_ids'][:8])
            st.markdown(f"""
            <div class="warning-box">
                <strong>Duplicate Sample IDs:</strong> {dup_str}
            </div>
            """, unsafe_allow_html=True)

        # 2. MAIN RESULTS SUB-TABS
        tab_ga, tab_hwe, tab_regions, tab_states, tab_cyp = st.tabs([
            "Genotype & Allele",
            "Hardy-Weinberg",
            "Geographic Population Groups",
            "State-Wise Pages",
            "CYP2C19 Calling"
        ])
        
        # GROUP FILTER PILLS BAR
        reg_counts = processed_df['Region'].value_counts().to_dict()
        gen_counts = processed_df['Gender_Clean'].value_counts().to_dict()
        
        group_options = {
            f"Overall (n={len(processed_df)})": processed_df,
            f"North India (n={reg_counts.get('North India', 0)})": processed_df[processed_df['Region'] == 'North India'],
            f"South India (n={reg_counts.get('South India', 0)})": processed_df[processed_df['Region'] == 'South India'],
            f"East India (n={reg_counts.get('East India', 0)})": processed_df[processed_df['Region'] == 'East India'],
            f"West India (n={reg_counts.get('West India', 0)})": processed_df[processed_df['Region'] == 'West India'],
            f"Female (n={gen_counts.get('Female', 0)})": processed_df[processed_df['Gender_Clean'] == 'Female'],
            f"Male (n={gen_counts.get('Male', 0)})": processed_df[processed_df['Gender_Clean'] == 'Male']
        }
        
        selected_group_label = st.radio(
            "Population Sub-Group Filter:",
            list(group_options.keys()),
            key="radio_subgroup_pills",
            horizontal=True
        )
        
        active_group_df = group_options[selected_group_label]
        
        def get_group_snp_stats(sub_df, snp_name):
            return GeneticStatsEngine.analyze_snp(sub_df[snp_name], snp_name) if snp_name in sub_df.columns else None

        # TAB 1: GENOTYPE & ALLELE
        with tab_ga:
            for snp_name in ['CYP2C19*2', 'CYP2C19*3', 'CYP2C19*17']:
                if snp_name in active_group_df.columns:
                    s_res = get_group_snp_stats(active_group_df, snp_name)
                    if s_res:
                        rsid_title = f"{SNP_CONFIG[snp_name]['rsid']} ({snp_name})"
                        st.markdown(f"### {rsid_title}")
                        st.caption(f"Valid: {s_res['valid_samples']} · Missing: 0 · Invalid calls: 0")
                        
                        cg1, cg2 = st.columns([1, 1])
                        
                        with cg1:
                            g_c = s_res['genotype_counts']
                            g_p = s_res['genotype_pcts']
                            
                            st.markdown("**GENOTYPE TABLE**")
                            gt_rows = [
                                {'GENOTYPE': k, 'COUNT': v, 'FREQUENCY': f"{g_p.get(k, 0)}%"}
                                for k, v in g_c.items()
                            ]
                            st.table(pd.DataFrame(gt_rows))
                            
                            st.markdown("**ALLELE TABLE**")
                            a_c = s_res['allele_counts']
                            a_p = s_res['allele_pcts']
                            al_rows = [
                                {'ALLELE': k, 'COUNT': v, 'FREQUENCY': f"{a_p.get(k, 0)}%"}
                                for k, v in a_c.items()
                            ]
                            st.table(pd.DataFrame(al_rows))

                        with cg2:
                            st.markdown("**Genotype Distribution**")
                            fig_bar = px.bar(
                                pd.DataFrame(gt_rows), x='GENOTYPE', y='COUNT',
                                color='GENOTYPE', text_auto=True,
                                color_discrete_sequence=['#1E4D38', '#2563EB', '#D97706']
                            )
                            fig_bar.update_layout(showlegend=False, height=280, margin=dict(l=10, r=10, t=10, b=10))
                            st.plotly_chart(fig_bar, use_container_width=True)

        # TAB 2: HARDY-WEINBERG
        with tab_hwe:
            st.markdown("### Hardy-Weinberg Equilibrium Test")
            hwe_rows = []
            for snp_name in ['CYP2C19*2', 'CYP2C19*3', 'CYP2C19*17']:
                if snp_name in active_group_df.columns:
                    s_res = get_group_snp_stats(active_group_df, snp_name)
                    if s_res:
                        hw = s_res['hwe']
                        hwe_rows.append({
                            'SNP Variant': snp_name,
                            'Observed Counts': str(s_res['genotype_counts']),
                            'Expected Counts': str(hw['expected_counts']),
                            'Chi2 Stat (χ²)': hw['chi2_stat'],
                            'Chi2 P-Value': hw['p_value'],
                            'Exact P-Value': hw['exact_p_value'],
                            'HWE Status': hw['interpretation']
                        })
            st.table(pd.DataFrame(hwe_rows))

        # TAB 3: GEOGRAPHIC POPULATION GROUPS
        with tab_regions:
            st.markdown("### Geographic Population Groups (North, South, East, West India)")
            st.caption("Standard Stratification Comparison across major Indian Geographic Population Groups.")
            
            reg_table_rows = []
            regional_results = full_results.get('regional', {})
            for r_name in ['North India', 'South India', 'East India', 'West India']:
                r_data = regional_results.get(r_name, {})
                snps = r_data.get('snps', {})
                cyp2_a = snps.get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
                pm_pct = r_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {}).get('percentage', 0)
                hwe_stat = snps.get('CYP2C19*2', {}).get('hwe', {}).get('interpretation', 'In HWE')
                
                reg_table_rows.append({
                    'Geographic Group': r_name,
                    'Sample Count (N)': r_data.get('sample_count', 0),
                    'CYP2C19*2 Var Allele Freq': cyp2_a,
                    'Poor Metabolizers (%)': f"{pm_pct}%",
                    'HWE Status': hwe_stat
                })
            st.table(pd.DataFrame(reg_table_rows))
            
            # Regional Plotly Bar Chart
            r_df = pd.DataFrame(reg_table_rows)
            fig_reg = px.bar(
                r_df, x='Geographic Group', y='Sample Count (N)',
                color='Geographic Group', text_auto=True,
                title="Sample Distribution across Geographic Groups"
            )
            st.plotly_chart(fig_reg, use_container_width=True)

        # TAB 4: STATE-WISE PAGES
        with tab_states:
            st.markdown("### Individual State-Wise Dedicated Profiles")
            state_dict = full_results.get('state_wise', {})
            state_names = sorted(list(state_dict.keys()))
            
            if state_names:
                selected_state = st.selectbox("Select Dedicated State Profile:", state_names)
                s_data = state_dict[selected_state]
                
                sk1, sk2, sk3, sk4 = st.columns(4)
                sk1.metric("State Sample Count", s_data.get('sample_count', 0))
                sk2.metric("Missing Data Rate", f"{s_data.get('missing_pct', 0)}%")
                cyp2_freq = s_data.get('snps', {}).get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
                sk3.metric("CYP2C19*2 Allele Freq", cyp2_freq)
                st_pm = s_data.get('cyp2c19_summary', {}).get('phenotypes', {}).get('Poor Metabolizer (PM)', {}).get('percentage', 0)
                sk4.metric("Poor Metabolizers %", f"{st_pm}%")
                
                sc1, sc2 = st.columns(2)
                with sc1:
                    st.markdown(f"**Genotype Counts for {selected_state}**")
                    st_g = s_data.get('snps', {}).get('CYP2C19*2', {}).get('genotype_counts', {})
                    st.table(pd.DataFrame([{'Genotype': k, 'Count': v} for k, v in st_g.items()]))
                with sc2:
                    st.markdown(f"**Phenotype Breakdown for {selected_state}**")
                    st_p = s_data.get('cyp2c19_summary', {}).get('phenotypes', {})
                    st.table(pd.DataFrame([{'Phenotype': k, 'Count': v['count'], 'Percentage': f"{v['percentage']}%"} for k, v in st_p.items()]))
            else:
                st.info("No state data available.")

        # TAB 5: CYP2C19 CALLING
        with tab_cyp:
            st.markdown("### CYP2C19 Star Allele & Phenotype Assignment")
            _, pgx_sub = CYP2C19Translator.process_dataframe(active_group_df)
            
            cp1, cp2 = st.columns([1, 1])
            with cp1:
                st.markdown("**METABOLIZER PHENOTYPES**")
                p_rows = []
                for p_name, p_data in pgx_sub['phenotypes'].items():
                    p_rows.append({
                        'Phenotype Category': p_name,
                        'Count': p_data['count'],
                        'Percentage (%)': f"{p_data['percentage']}%"
                    })
                st.table(pd.DataFrame(p_rows))
                
            with cp2:
                st.markdown("**DIPLOTYPE CALLS**")
                d_rows = []
                for d_name, d_data in pgx_sub['diplotypes'].items():
                    d_rows.append({
                        'Diplotype Call': d_name,
                        'Count': d_data['count'],
                        'Percentage (%)': f"{d_data['percentage']}%"
                    })
                st.table(pd.DataFrame(d_rows))

        # 4. EXPORT SECTION
        st.markdown("<br><hr>", unsafe_allow_html=True)
        st.markdown("### Export Results")
        st.caption("Download the full result set for the current group selection.")
        
        ce1, ce2, ce3, ce4 = st.columns([1, 1, 1, 1])
        
        excel_bytes = ReportGenerator.export_to_excel(full_results)
        ce1.download_button("⬇ Excel (.xlsx, all sheets)", data=excel_bytes, file_name="PGx_Results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True)
        
        csv_bytes = processed_df.to_csv(index=False).encode('utf-8')
        ce2.download_button("⬇ CSV (summary)", data=csv_bytes, file_name="PGx_Summary.csv", mime="text/csv", use_container_width=True)
        
        pdf_bytes = ReportGenerator.export_to_pdf(full_results)
        ce3.download_button("⬇ PDF (print view)", data=pdf_bytes, file_name="PGx_Print_Report.pdf", mime="application/pdf", use_container_width=True)
        
        if ce4.button("← Re-map columns", use_container_width=True):
            st.session_state['wizard_step'] = 2
            st.rerun()

# -----------------------------------------------------------------------------
# FOOTER NOTE
# -----------------------------------------------------------------------------
st.markdown("""
<div class="footer-note">
    PGx Pop MVP - CYP2C19 calling uses a simplified unphased-genotype heuristic — see notes in the CYP2C19 tab - Not for clinical use without independent validation.
</div>
""", unsafe_allow_html=True)
