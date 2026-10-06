FROM python:3.11-slim-bookworm

ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONPATH=/home/appuser

ARG USER_ID=1000
ARG GROUP_ID=1000

RUN groupadd -g ${GROUP_ID} appuser && \
    useradd --create-home --shell /bin/bash -u ${USER_ID} -g ${GROUP_ID} appuser && \
    chown -R appuser:appuser /home/appuser

WORKDIR /home/appuser/app

# Build context is the repository root.
COPY api/requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY api/ .
COPY gene_corpus /home/appuser/gene_corpus

USER appuser

EXPOSE 5000

CMD ["gunicorn", "-k", "uvicorn.workers.UvicornWorker", "-w", "2", "-b", "0.0.0.0:5000", "--timeout", "120", "--access-logfile", "-", "--error-logfile", "-", "main:app"]
