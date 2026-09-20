"""
Report Exporter Engine for Multi-Tab Excel, CSV, and PDF Output.
Generates comprehensive pharmacogenomic analysis reports.
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
                
            # 4. Gender Tabs
            gender_results = full_results.get('gender', {})
            for gen_name, gen_data in gender_results.items():
                sheet_title = f"Gender_{gen_name}"[:31]
                gen_df = cls._dict_to_summary_table(gen_data)
                gen_df.to_excel(writer, sheet_name=sheet_title, index=False)

        return output.getvalue()

    @classmethod
    def export_to_pdf(cls, full_results: Dict[str, Any]) -> bytes:
        """Generates a publication-ready PDF summary report."""
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36)
        elements = []
        
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            'TitleStyle',
            parent=styles['Heading1'],
            fontSize=18,
            textColor=colors.HexColor('#1E3A8A'),
            spaceAfter=12
        )
        heading_style = ParagraphStyle(
            'HeadingStyle',
            parent=styles['Heading2'],
            fontSize=14,
            textColor=colors.HexColor('#1E40AF'),
            spaceBefore=10,
            spaceAfter=8
        )
        body_style = styles['BodyText']
        
        # Title & Subtitle
        elements.append(Paragraph("Population Pharmacogenomics Analysis Report", title_style))
        elements.append(Paragraph("Automated Analysis of Genotype, Allele, HWE, and Phenotype Distributions", body_style))
        elements.append(Spacer(1, 12))
        
        # Summary Overview Table
        overall_stats = full_results.get('overall', {})
        total_samples = overall_stats.get('sample_count', 0)
        
        table_data = [
            ["Metric", "Value"],
            ["Total Samples Analyzed", str(total_samples)],
            ["Geographic Subgroups", "North India, South India, East India, West India"],
            ["Target Pharmacogene", "CYP2C19 (*2, *3, *17)"]
        ]
        
        t = Table(table_data, colWidths=[200, 300])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (1, 0), colors.HexColor('#3B82F6')),
            ('TEXTCOLOR', (0, 0), (1, 0), colors.white),
            ('FONTNAME', (0, 0), (1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (1, 0), 6),
            ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#F3F4F6')),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#D1D5DB'))
        ]))
        elements.append(t)
        elements.append(Spacer(1, 16))
        
        # Phenotype Breakdown Section
        elements.append(Paragraph("CYP2C19 Phenotype Summary (Overall Population)", heading_style))
        pgx = overall_stats.get('cyp2c19_summary', {}).get('phenotypes', {})
        
        p_table_data = [["Metabolizer Phenotype", "Count", "Percentage (%)"]]
        for pheno_name, pdata in pgx.items():
            p_table_data.append([pheno_name, str(pdata.get('count', 0)), f"{pdata.get('percentage', 0)}%"])
            
        ptable = Table(p_table_data, colWidths=[250, 100, 150])
        ptable.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E40AF')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#9CA3AF'))
        ]))
        elements.append(ptable)
        elements.append(Spacer(1, 16))

        doc.build(elements)
        return buffer.getvalue()
