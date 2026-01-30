"""
GeneHealth Analysis Service API
FastAPI server for processing genetic data files
"""

import os
import asyncio
import httpx
from typing import Optional
from fastapi import FastAPI, HTTPException, Header, BackgroundTasks
from pydantic import BaseModel

from scripts.run_full_analysis import run_analysis
from scripts.comprehensive_snp_database import SNP_DATABASE
from scripts.disease_risk_analyzer import analyze_disease_risks
from scripts.generate_exhaustive_report import generate_reports

app = FastAPI(title="GeneHealth Analysis Service")

# API Key for authentication
API_KEY = os.environ.get("ANALYSIS_SERVICE_API_KEY", "dev-key")


class AnalysisRequest(BaseModel):
    jobId: str
    sourceFormat: str
    genomeContent: Optional[str] = None  # Direct content (small files)
    fileUrl: Optional[str] = None        # URL to download (large files)
    callbackUrl: str
    apiKey: Optional[str] = None         # API key for downloading


class AnalysisProgress(BaseModel):
    jobId: str
    status: str
    progress: int
    currentStep: str
    error: Optional[str] = None


class AnalysisResult(BaseModel):
    jobId: str
    status: str
    snpCount: Optional[int] = None
    findingsSummary: Optional[dict] = None
    error: Optional[str] = None


def verify_api_key(x_api_key: str = Header(...)):
    if x_api_key != API_KEY:
        raise HTTPException(status_code=401, detail="Invalid API key")
    return x_api_key


@app.get("/health")
async def health_check():
    return {"status": "ok", "snp_database_size": len(SNP_DATABASE)}


@app.post("/analyze")
async def analyze(
    request: AnalysisRequest,
    background_tasks: BackgroundTasks,
    x_api_key: str = Header(...)
):
    """Start genome analysis in background"""
    verify_api_key(x_api_key)

    # Run analysis in background
    background_tasks.add_task(
        process_analysis,
        request.jobId,
        request.sourceFormat,
        request.genomeContent,
        request.fileUrl,
        request.apiKey or x_api_key,
        request.callbackUrl
    )

    return {"status": "processing", "jobId": request.jobId}


async def process_analysis(
    job_id: str,
    source_format: str,
    genome_content: Optional[str],
    file_url: Optional[str],
    api_key: str,
    callback_url: str
):
    """Process the genome analysis and send callback"""
    try:
        # Get genome content - either directly or by downloading
        if genome_content:
            content = genome_content
        elif file_url:
            print(f"Downloading genome file from: {file_url}")
            async with httpx.AsyncClient() as client:
                response = await client.get(
                    file_url,
                    headers={"X-API-Key": api_key},
                    timeout=120.0
                )
                response.raise_for_status()
                content = response.text
                print(f"Downloaded {len(content)} bytes")
        else:
            raise ValueError("No genome content or file URL provided")

        # Run the full analysis pipeline
        result = await asyncio.to_thread(
            run_analysis,
            content,
            source_format
        )

        reports = result.get("reports", {})
        print(f"Analysis complete: {result['snp_count']} SNPs, {len(reports)} reports")
        print(f"Report types generated: {list(reports.keys())}")

        # Send success callback with reports
        async with httpx.AsyncClient() as client:
            await client.post(
                callback_url,
                json={
                    "jobId": job_id,
                    "status": "completed",
                    "snpCount": result["snp_count"],
                    "findingsSummary": result["findings_summary"],
                    "reports": result.get("reports", {})
                },
                headers={"X-Analysis-Service-Key": API_KEY},
                timeout=60.0  # Increased timeout for larger payload
            )
            print(f"Callback sent successfully for job {job_id}")

    except Exception as e:
        # Send failure callback
        async with httpx.AsyncClient() as client:
            await client.post(
                callback_url,
                json={
                    "jobId": job_id,
                    "status": "failed",
                    "error": str(e)
                },
                headers={"X-Analysis-Service-Key": API_KEY},
                timeout=30.0
            )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8080)
