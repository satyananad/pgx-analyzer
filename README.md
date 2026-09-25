---
title: PGx Analyzer
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: streamlit
sdk_version: 1.41.1
app_file: app.py
pinned: false
---

# 🧬 Automated Population Pharmacogenomics Analysis Platform

An industry-level, production-grade automated platform designed for genotype, allele frequency, Hardy–Weinberg Equilibrium (HWE), regional/gender demographic stratification, and CYP2C19 star allele/diplotype/phenotype classification.

---

## 📐 Project Architecture & Directory Structure

```
analysis tools/
├── app.py                      # Streamlit & HTML/CSS Web Application Dashboard
├── config.py                   # Global configuration & reference rules (State mappings, SNP rules, CPIC matrix)
├── core/                       # Core analytical & statistical engine package
│   ├── __init__.py
│   ├── qc.py                   # Data Quality Control, Validation & Cleaning Engine
│   ├── stats.py                # Genotype, Allele Frequency & HWE (Chi-Square & Fisher Exact) Engine
│   ├── stratification.py       # Geographic (North/South/East/West India) & Gender Stratifier
│   └── cyp2c19.py              # CYP2C19 Star Allele, Diplotype & CPIC Phenotype Classifier
├── database/                   # Data Persistence Layer
│   ├── __init__.py
│   └── db_manager.py           # SQLite database manager for 1,000 to 10,000+ sample batch growth
├── reports/                    # Exporter & Reporting Layer
│   ├── __init__.py
│   └── report_generator.py     # Multi-Tab Excel, CSV, and Publication PDF Exporter
├── tests/                      # Automated Unit Test Suite (PyTest)
│   ├── __init__.py
│   ├── test_qc.py
│   ├── test_stats.py
│   └── test_cyp2c19.py
├── requirements.txt            # Python dependencies
└── README.md                   # Full documentation & usage guide
```

---

## ⚡ Key Features

1. **Automated Quality Control & Validation (`core/qc.py`)**:
   - Standardizes headers (`sample ID`, `Gender`, `State`, `CYP2C19*2`, etc.).
   - Normalizes gender (`Male`, `Female`, `Unknown`).
   - Detects duplicate Sample IDs across uploaded files and historical database records.
   - Cleans genotype strings and flags invalid alleles.

2. **Population Genetics Engine (`core/stats.py`)**:
   - Calculates genotype counts ($N_{AA}, N_{Aa}, N_{aa}$), frequencies, and percentages.
   - Calculates allele counts ($N_{Ref}, N_{Var}$) and frequencies ($p$ and $q$).
   - Calculates Hardy-Weinberg Equilibrium (HWE): Expected counts ($E_{ij}$), Chi-Square statistic ($\chi^2$), P-value, and exact Haldane/Wigginton test fallback for small sample sizes ($E < 5$).

3. **Demographic & Regional Stratification (`core/stratification.py`)**:
   - Categorizes samples into **North India, South India, East India, and West India**.
   - Features a 3-tier mapping engine: Direct State match $\rightarrow$ Native Place/City gazetteer lookup $\rightarrow$ Region string parsing.
   - Computes separate statistical and pharmacogenomic parameters for Overall, Regional, and Gender (Male vs Female) subgroups.

4. **CYP2C19 Star Allele & Phenotype Translation (`core/cyp2c19.py`)**:
   - Translates *2 (rs4244285), *3 (rs4986893), and *17 (rs12248560) variant combinations into maternal/paternal diplotypes (`*1/*1`, `*1/*2`, `*1/*17`, `*2/*17`, `*2/*2`, `*17/*17`, etc.).
   - Maps diplotypes to CPIC Metabolizer Phenotypes:
     - **Ultrarapid Metabolizer (UM)**
     - **Rapid Metabolizer (RM)**
     - **Normal Metabolizer (NM)**
     - **Intermediate Metabolizer (IM)**
     - **Poor Metabolizer (PM)**

5. **Scalable Database Backend (`database/db_manager.py`)**:
   - SQLite persistence layer for growing datasets (1,000 $\rightarrow$ 1,200 $\rightarrow$ 2,000 $\rightarrow$ 10,000+ samples).
   - Instant deduplication and automatic recalculation across historical + new sample batches.

6. **Interactive Dashboard & Automated Exporters (`app.py` & `reports/report_generator.py`)**:
   - Streamlit web interface with interactive tables and regional comparative bar charts.
   - One-click multi-tab Excel export (separate worksheets for Overall, North, South, East, West, Gender, and Full Dataset).
   - One-click publication-ready PDF summary report.

---

## 🚀 How to Run the Application

### 1. Execute Web Dashboard
```bash
python3 -m streamlit run app.py
```

### 2. Run Automated Unit Tests
```bash
python3 -m pytest tests/
```

---

## 🔬 Mathematical Formulas

### 1. Allele Frequencies
$$p = \frac{2 N_{AA} + N_{Aa}}{2 N}, \quad q = 1 - p$$

### 2. Expected HWE Counts
$$E_{AA} = N \cdot p^2, \quad E_{Aa} = N \cdot 2pq, \quad E_{aa} = N \cdot q^2$$

### 3. Chi-Square Goodness-of-Fit Test
$$\chi^2 = \sum \frac{(O - E)^2}{E} = \frac{(N_{AA} - E_{AA})^2}{E_{AA}} + \frac{(N_{Aa} - E_{Aa})^2}{E_{Aa}} + \frac{(N_{aa} - E_{aa})^2}{E_{aa}}$$
$$\text{P-value} = P(\chi^2_1 \ge \text{stat})$$
- **HWE Status**: *In HWE* if $P \ge 0.05$; *Departure from HWE* if $P < 0.05$.
