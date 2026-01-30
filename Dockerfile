FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY api_server.py .
COPY scripts/ ./scripts/

# Create data directory and copy ClinVar database
RUN mkdir -p /app/data
COPY data/clinvar_alleles.tsv.gz /app/data/
RUN gunzip /app/data/clinvar_alleles.tsv.gz

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PORT=8080
ENV CLINVAR_PATH=/app/data/clinvar_alleles.tsv

# Expose port
EXPOSE 8080

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8080/health || exit 1

# Run the application
CMD ["uvicorn", "api_server:app", "--host", "0.0.0.0", "--port", "8080"]
