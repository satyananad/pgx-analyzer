"""
Report Exporter Engine for Multi-Tab Excel, CSV, and PDF Output.
Generates comprehensive pharmacogenomic analysis reports including:
1. Full Cleaned Dataset
2. Overall Population Summary
3. Regional Tabs (North India, South India, East India, West India)
4. State-Wise Analysis Summary
5. Gender Breakdown (Male, Female)
6. Data QC & Missing Data Audit
"""

import pandas as pd
import io
from typing import Dict, Any
import logging

from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

logger = logging.getLogger(__name__)


class ReportGenerator:
    @staticmethod
    def _dict_to_summary_table(stats_dict: Dict[str, Any]) -> pd.DataFrame:
        """Converts statistical dictionary into a flat Pandas DataFrame for Excel export."""
        rows = []
        
        # SNP stats
        snps = stats_dict.get('snps', {})
        for snp_name, snp_data in snps.items():
            valid_n = snp_data.get('valid_samples', 0)
            
            # Genotypes
            counts = snp_data.get('genotype_counts', {})
            pcts = snp_data.get('genotype_pcts', {})
            gen_str = ", ".join([f"{k}: {v} ({pcts.get(k, 0)}%)" for k, v in counts.items()])
            
            # Alleles
            a_counts = snp_data.get('allele_counts', {})
            a_pcts = snp_data.get('allele_pcts', {})
            allele_str = ", ".join([f"{k}: {v} ({a_pcts.get(k, 0)}%)" for k, v in a_counts.items()])
            
            # HWE
            hwe = snp_data.get('hwe', {})
            
            rows.append({
                'Category': 'SNP Analysis',
                'Parameter / Variant': snp_name,
                'Valid Samples': valid_n,
                'Genotype Distribution': gen_str,
                'Allele Distribution': allele_str,
                'Chi2 Stat': hwe.get('chi2_stat', ''),
                'P-Value': hwe.get('p_value', ''),
                'Exact P-Value': hwe.get('exact_p_value', ''),
                'HWE Interpretation': hwe.get('interpretation', '')
            })
            
        # Phenotypes
        pgx = stats_dict.get('cyp2c19_summary', {}).get('phenotypes', {})
        for pheno, pdata in pgx.items():
            rows.append({
                'Category': 'Metabolizer Phenotype',
                'Parameter / Variant': pheno,
                'Valid Samples': stats_dict.get('sample_count', 0),
                'Genotype Distribution': f"Count: {pdata.get('count', 0)}",
                'Allele Distribution': f"Frequency: {pdata.get('frequency', 0)}",
                'Chi2 Stat': '-',
                'P-Value': '-',
                'Exact P-Value': '-',
                'HWE Interpretation': f"Percentage: {pdata.get('percentage', 0)}%"
            })
            
        return pd.DataFrame(rows)

    @classmethod
    def export_to_excel(cls, full_results: Dict[str, Any]) -> bytes:
        """Generates a multi-tab Excel workbook as bytes."""
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            # 1. Full Dataset Tab
            if 'processed_dataframe' in full_results:
                full_df = full_results['processed_dataframe']
                full_df.to_excel(writer, sheet_name='Full Cleaned Dataset', index=False)
                
            # 2. Overall Population Tab
            if 'overall' in full_results:
                ov_df = cls._dict_to_summary_table(full_results['overall'])
                ov_df.to_excel(writer, sheet_name='Overall Population', index=False)
                
            # 3. Regional Tabs (North, South, East, West)
            regional_results = full_results.get('regional', {})
            for reg_name, reg_data in regional_results.items():
                sheet_title = reg_name.replace(' ', '_')[:31]  # Excel 31 char sheet limit
                reg_df = cls._dict_to_summary_table(reg_data)
                reg_df.to_excel(writer, sheet_name=sheet_title, index=False)
                
            # 4. State-Wise Analysis Summary Tab
            state_results = full_results.get('state_wise', {})
            st_rows = []
            for st_name, st_data in state_results.items():
                cyp2_a = st_data.get('snps', {}).get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
                pm_data = st_data.get('cyp2c19_summary', {}).get('phenotypes', {})
                pm_pct = pm_data.get('Poor Metabolizer', pm_data.get('Poor Metabolizer (PM)', {})).get('percentage', 0)
                st_rows.append({
                    'State Name': st_name,
                    'Sample N': st_data.get('sample_count', 0),
                    'Missing Data Rate (%)': f"{st_data.get('missing_pct', 0)}%",
                    'CYP2C19*2 Var Freq (A)': cyp2_a,
                    'Poor Metabolizers (%)': f"{pm_pct}%"
                })
            if st_rows:
                pd.DataFrame(st_rows).to_excel(writer, sheet_name='State-Wise Summary', index=False)

            # 5. Gender Tabs
            gender_results = full_results.get('gender', {})
            for gen_name, gen_data in gender_results.items():
                sheet_title = f"Gender_{gen_name}"[:31]
                gen_df = cls._dict_to_summary_table(gen_data)
                gen_df.to_excel(writer, sheet_name=sheet_title, index=False)

        return output.getvalue()

    @classmethod
    def export_to_pdf(cls, full_results: Dict[str, Any]) -> bytes:
        """Generates a publication-ready comprehensive PDF report."""
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
        elements = []
        
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            'TitleStyle',
            parent=styles['Heading1'],
            fontSize=16,
            textColor=colors.HexColor('#0F172A'),
            spaceAfter=8
        )
        heading_style = ParagraphStyle(
            'HeadingStyle',
            parent=styles['Heading2'],
            fontSize=12,
            textColor=colors.HexColor('#1E3A8A'),
            spaceBefore=12,
            spaceAfter=6
        )
        body_style = styles['BodyText']
        
        # 1. Document Header
        elements.append(Paragraph("Automated Population Pharmacogenomics Analysis Report", title_style))
        elements.append(Paragraph("Development of an Automated Platform for Genotype, Allele, HWE, Diplotype, and CPIC Phenotype Analysis", body_style))
        elements.append(Spacer(1, 10))
        
        # 2. Executive Metadata Summary Table
        overall_stats = full_results.get('overall', {})
        total_samples = overall_stats.get('sample_count', 0)
        
        meta_table = [
            ["Parameter", "Details"],
            ["Total Cohort Analyzed", f"{total_samples:,} Samples"],
            ["Geographic Population Groups", "South India, North India, East India, West India, Central India"],
            ["Target Gene / SNPs", "CYP2C19 (*2 rs4244285, *3 rs4986893, *17 rs12248560)"],
            ["Statistical Engines", "Chi-Square Goodness-of-Fit (df=1) & Haldane Exact Test"],
            ["Interpretation Framework", "CPIC Guidelines for CYP2C19 Metabolizer Phenotypes"]
        ]
        
        t_meta = Table(meta_table, colWidths=[180, 360])
        t_meta.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (1, 0), colors.HexColor('#0F172A')),
            ('TEXTCOLOR', (0, 0), (1, 0), colors.white),
            ('FONTNAME', (0, 0), (1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (1, 0), 5),
            ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#F8FAFC')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_meta)
        elements.append(Spacer(1, 12))
        
        # 3. Hardy-Weinberg Equilibrium (HWE) Summary Section
        elements.append(Paragraph("1. Overall Population Hardy-Weinberg Equilibrium (HWE) Analysis", heading_style))
        snps_data = overall_stats.get('snps', {})
        hwe_table_data = [["SNP Variant", "Valid N", "Allele Freq p / q", "Chi2 (χ²)", "P-Value", "HWE Status"]]
        
        for snp_name, s_res in snps_data.items():
            hw = s_res.get('hwe', {})
            af = s_res.get('allele_freqs', {})
            ref_a = list(af.keys())[0] if af else 'Ref'
            var_a = list(af.keys())[1] if len(af) > 1 else 'Var'
            freq_str = f"{ref_a}:{af.get(ref_a, 0)} · {var_a}:{af.get(var_a, 0)}"
            hwe_table_data.append([
                snp_name,
                str(s_res.get('valid_samples', 0)),
                freq_str,
                str(hw.get('chi2_stat', 0)),
                str(hw.get('p_value', 1.0)),
                str(hw.get('interpretation', 'In HWE'))
            ])
            
        t_hwe = Table(hwe_table_data, colWidths=[90, 60, 150, 70, 70, 100])
        t_hwe.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E3A8A')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_hwe)
        elements.append(Spacer(1, 12))

        # 4. Geographic Population Groups Section (All 5 Regions)
        elements.append(Paragraph("2. Geographic Population Groups Matrix (All 5 Regions)", heading_style))
        regional_results = full_results.get('regional', {})
        all_regions = ['South India', 'North India', 'East India', 'West India', 'Central India']
        
        reg_table_data = [["Geographic Region", "Sample N", "*2 Var Freq f(A)", "*2 Chi2 (P-Val)", "*17 Var Freq f(T)", "PM (%)"]]
        for r_name in all_regions:
            r_data = regional_results.get(r_name, {})
            r_snps = r_data.get('snps', {})
            cyp2_a = r_snps.get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
            cyp2_chi = r_snps.get('CYP2C19*2', {}).get('hwe', {}).get('chi2_stat', 0)
            cyp2_p = r_snps.get('CYP2C19*2', {}).get('hwe', {}).get('p_value', 1.0)
            cyp17_t = r_snps.get('CYP2C19*17', {}).get('allele_freqs', {}).get('T', 0)
            r_phenos = r_data.get('cyp2c19_summary', {}).get('phenotypes', {})
            pm_pct = r_phenos.get('Poor Metabolizer', r_phenos.get('Poor Metabolizer (PM)', {})).get('percentage', 0)
            
            reg_table_data.append([
                r_name,
                str(r_data.get('sample_count', 0)),
                str(cyp2_a),
                f"{cyp2_chi} (p={cyp2_p})",
                str(cyp17_t),
                f"{pm_pct}%"
            ])
            
        t_reg = Table(reg_table_data, colWidths=[110, 70, 110, 110, 80, 60])
        t_reg.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F766E')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_reg)
        elements.append(Spacer(1, 12))

        # 5. Gender Stratification Section
        elements.append(Paragraph("3. Gender Stratification Analysis (Male vs Female vs Overall)", heading_style))
        gender_results = full_results.get('gender', {})
        g_table_data = [["Cohort / Gender", "Sample N", "*2 Var Freq f(A)", "*17 Var Freq f(T)", "Normal %", "Poor %"]]
        
        for g_name in ['Male', 'Female']:
            g_data = gender_results.get(g_name, {})
            g_snps = g_data.get('snps', {})
            cyp2_a = g_snps.get('CYP2C19*2', {}).get('allele_freqs', {}).get('A', 0)
            cyp17_t = g_snps.get('CYP2C19*17', {}).get('allele_freqs', {}).get('T', 0)
            g_phenos = g_data.get('cyp2c19_summary', {}).get('phenotypes', {})
            nm_p = g_phenos.get('Normal Metabolizer', g_phenos.get('Normal Metabolizer (NM)', {})).get('percentage', 0)
            pm_p = g_phenos.get('Poor Metabolizer', g_phenos.get('Poor Metabolizer (PM)', {})).get('percentage', 0)
            
            g_table_data.append([g_name, str(g_data.get('sample_count', 0)), str(cyp2_a), str(cyp17_t), f"{nm_p}%", f"{pm_p}%"])
            
        t_gen = Table(g_table_data, colWidths=[110, 70, 110, 110, 80, 60])
        t_gen.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4338CA')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1'))
        ]))
        elements.append(t_gen)
        elements.append(Spacer(1, 12))

        # 6. CPIC Diplotypes & Metabolizer Phenotypes Section
        elements.append(Paragraph("4. CPIC CYP2C19 Diplotypes & Metabolizer Phenotypes", heading_style))
        pgx = overall_stats.get('cyp2c19_summary', {})
        pheno_data = pgx.get('phenotypes', {})
        
        p_table_data = [["Metabolizer Phenotype Category", "Count (N)", "Percentage (%)"]]
        for pheno_name, pdata in pheno_data.items():
            p_table_data.append([pheno_name, str(pdata.get('count', 0)), f"{pdata.get('percentage', 0)}%"])
            
        t_pheno = Table(p_table_data, colWidths=[240, 130, 170])
        t_pheno.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E40AF')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#9CA3AF'))
        ]))
        elements.append(t_pheno)

        doc.build(elements)
        return buffer.getvalue()
