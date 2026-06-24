#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
REF_DIR="${PROJECT_DIR}/reference"
INDEX_DIR="${REF_DIR}/salmon_index"

TRANSCRIPTOME_FASTA="${REF_DIR}/gencode.v45.transcripts.fa.gz"
THREADS=${THREADS:-8}

mkdir -p "${INDEX_DIR}"

if [[ ! -f "${TRANSCRIPTOME_FASTA}" ]]; then
  echo "[ERROR] Transcriptome FASTA not found: ${TRANSCRIPTOME_FASTA}" >&2
  exit 1
fi

salmon index \
  -t "${TRANSCRIPTOME_FASTA}" \
  -i "${INDEX_DIR}" \
  -k 31 \
  -p "${THREADS}"

echo "[INFO] Salmon index created: ${INDEX_DIR}"
