"""
FastAPI Server for Custom Tailwind CSS + Chart.js Medical SaaS Interface.
Provides REST APIs for file uploads, QC, HWE math, demographic stratification,
state-wise analysis, and multi-format report exports.
"""

from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.responses import HTMLResponse, Response, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
import pandas as pd
import numpy as np
import io
import os
import json
from typing import Dict, Any, Optional

from core.qc import DataQCEngine
from core.stratification import DemographicStratifier
from core.stats import GeneticStatsEngine
from core.cyp2c19 import CYP2C19Translator
from database.db_manager import DatabaseManager
from reports.report_generator import ReportGenerator
from config import SNP_CONFIG, PHENOTYPE_ORDER

app = FastAPI(title="PGx Analytics SaaS API", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

db_manager = DatabaseManager()

# Ensure static directory exists
os.makedirs("static", exist_ok=True)
app.mount("/static", StaticFiles(directory="static"), name="static")


@app.get("/", response_class=HTMLResponse)
async def serve_frontend():
    index_path = os.path.join("static", "index.html")
    if os.path.exists(index_path):
        with open(index_path, "r", encoding="utf-8") as f:
            return f.read()
    return "<h1>PGx Analytics SaaS API is running. Please upload static/index.html</h1>"


@app.get("/api/demo")
async def load_demo_data():
    """Loads demo workspace dataset and returns processed analytics JSON."""
    demo_file = "categorized by state into North, South, East, West .xlsx"
    if not os.path.exists(demo_file):
        raise HTTPException(status_code=404, detail="Demo dataset not found")
        
    df = pd.read_excel(demo_file)
    existing_ids = db_manager.get_existing_sample_ids()
    
    qc_engine = DataQCEngine(df, existing_sample_ids=existing_ids)
    clean_df, qc_report = qc_engine.run_qc()
    
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    
    processed_df = full_results['processed_dataframe'].copy().fillna("")
    
    # Store/merge new samples into database
    try:
        db_manager.insert_batch(processed_df)
    except Exception as db_err:
        print(f"Database batch insertion note: {db_err}")
    
    return {
        "status": "success",
        "qc_report": qc_report,
        "overall": full_results['overall'],
        "regional": full_results['regional'],
        "state_wise": full_results['state_wise'],
        "gender": full_results['gender'],
        "sample_data": processed_df.to_dict(orient="records")
    }


@app.post("/api/upload")
async def upload_dataset(file: UploadFile = File(...)):
    """Ingests uploaded Excel/CSV file and returns full analytical JSON."""
    contents = await file.read()
    file_bytes = io.BytesIO(contents)
    
    try:
        if file.filename.endswith('.csv'):
            df = pd.read_csv(file_bytes)
        else:
            df = pd.read_excel(file_bytes)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Invalid file format: {str(e)}")

    existing_ids = db_manager.get_existing_sample_ids()
    qc_engine = DataQCEngine(df, existing_sample_ids=existing_ids)
    clean_df, qc_report = qc_engine.run_qc()
    
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    
    processed_df = full_results['processed_dataframe'].copy().fillna("")
    
    # Store/merge new samples into database
    try:
        inserted_count, skipped_count = db_manager.insert_batch(processed_df)
        qc_report["db_inserted_samples"] = inserted_count
        qc_report["db_skipped_duplicates"] = skipped_count
    except Exception as db_err:
        print(f"Database batch insertion note: {db_err}")
    
    return {
        "status": "success",
        "filename": file.filename,
        "qc_report": qc_report,
        "overall": full_results['overall'],
        "regional": full_results['regional'],
        "state_wise": full_results['state_wise'],
        "gender": full_results['gender'],
        "sample_data": processed_df.to_dict(orient="records")
    }


@app.get("/api/export/excel")
async def export_excel():
    """Exports multi-tab Excel workbook."""
    demo_file = "categorized by state into North, South, East, West .xlsx"
    df = pd.read_excel(demo_file) if os.path.exists(demo_file) else pd.DataFrame()
    qc_engine = DataQCEngine(df)
    clean_df, _ = qc_engine.run_qc()
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    
    excel_bytes = ReportGenerator.export_to_excel(full_results)
    return Response(
        content=excel_bytes,
        media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        headers={"Content-Disposition": "attachment; filename=PGx_Population_Report.xlsx"}
    )


@app.get("/api/export/pdf")
async def export_pdf():
    """Exports PDF summary report."""
    demo_file = "categorized by state into North, South, East, West .xlsx"
    df = pd.read_excel(demo_file) if os.path.exists(demo_file) else pd.DataFrame()
    qc_engine = DataQCEngine(df)
    clean_df, _ = qc_engine.run_qc()
    full_results = DemographicStratifier.stratify_and_analyze(clean_df)
    
    pdf_bytes = ReportGenerator.export_to_pdf(full_results)
    return Response(
        content=pdf_bytes,
        media_type="application/pdf",
        headers={"Content-Disposition": "attachment; filename=PGx_Executive_Summary.pdf"}
    )


if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("PORT", 8000))
    uvicorn.run("server:app", host="0.0.0.0", port=port, reload=False)
