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
    genomeContent: str
    callbackUrl: str


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
        request.callbackUrl
    )

    return {"status": "processing", "jobId": request.jobId}


async def process_analysis(
    job_id: str,
    source_format: str,
    genome_content: str,
    callback_url: str
):
    """Process the genome analysis and send callback"""
    try:
        # Run the full analysis pipeline
        result = await asyncio.to_thread(
            run_analysis,
            genome_content,
            source_format
        )

        # Send success callback
        async with httpx.AsyncClient() as client:
            await client.post(
                callback_url,
                json={
                    "jobId": job_id,
                    "status": "completed",
                    "snpCount": result["snp_count"],
                    "findingsSummary": result["findings_summary"]
                },
                headers={"X-Analysis-Service-Key": API_KEY},
                timeout=30.0
            )

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
