# Multi-stage: build gffread, then slim runtime with pysam + API.
FROM debian:bookworm-slim AS gffread-build

ENV DEBIAN_FRONTEND=noninteractive
# Pin known-good gffread (gclib submodule required).
ARG GFFREAD_REF=v0.12.7
ARG GCLIB_REF=v0.12.7

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        g++ \
        make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
RUN git clone --depth 1 --branch "${GFFREAD_REF}" \
        https://github.com/gpertea/gffread.git gffread \
    && git clone --depth 1 --branch "${GCLIB_REF}" \
        https://github.com/gpertea/gclib.git gclib \
    && make -C gffread \
    && test -x gffread/gffread \
    && (./gffread/gffread -h 2>&1 || true) | grep -E -- '-x' \
    && (./gffread/gffread -h 2>&1 || true) | grep -E -- '-w' \
    && (./gffread/gffread -h 2>&1 || true) | grep -E -- '-y' \
    && (./gffread/gffread -h 2>&1 || true) | grep -E -- '-W' \
    && (./gffread/gffread -h 2>&1 || true) | grep -F -- '--w-add' \
    && cp gffread/gffread /usr/local/bin/gffread


FROM python:3.11-slim-bookworm

ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONPATH=/home/appuser

ARG USER_ID=1000
ARG GROUP_ID=1000

RUN groupadd -g ${GROUP_ID} appuser && \
    useradd --create-home --shell /bin/bash -u ${USER_ID} -g ${GROUP_ID} appuser && \
    chown -R appuser:appuser /home/appuser

# Runtime libs for pysam (htslib) + ca-certs for NCBI HTTPS.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4 \
        libssl3 \
        zlib1g \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=gffread-build /usr/local/bin/gffread /usr/local/bin/gffread
RUN gffread -h >/dev/null 2>&1 || gffread -h 2>&1 | head -5

WORKDIR /home/appuser/app

# Build context is the repository root.
COPY api/requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY api/ .
COPY gene_corpus /home/appuser/gene_corpus

USER appuser

EXPOSE 5000

CMD ["gunicorn", "-k", "uvicorn.workers.UvicornWorker", "-w", "2", "-b", "0.0.0.0:5000", "--timeout", "900", "--access-logfile", "-", "--error-logfile", "-", "main:app"]
