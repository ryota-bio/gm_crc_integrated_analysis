#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
META="${PROJECT_DIR}/metadata/sample_metadata.csv"
TRIM_DIR="${PROJECT_DIR}/results/fastp_trimmed"
INDEX_DIR="${PROJECT_DIR}/reference/salmon_index"
OUT_DIR="${PROJECT_DIR}/results/salmon_quant"

THREADS=${THREADS:-8}

mkdir -p "${OUT_DIR}"

if [[ ! -d "${INDEX_DIR}" ]]; then
  echo "[ERROR] Salmon index directory not found: ${INDEX_DIR}" >&2
  echo "[INFO] Please run scripts/02_salmon_index.sh first." >&2
  exit 1
fi

tail -n +2 "${META}" | while IFS=',' read -r sample_id condition time_point replicate cell_line bacteria MOI infection_time fastq_R1 fastq_R2
do
  echo "[Salmon] Quantifying ${sample_id}"

  R1="${TRIM_DIR}/${sample_id}_R1.trimmed.fastq.gz"
  R2="${TRIM_DIR}/${sample_id}_R2.trimmed.fastq.gz"

  if [[ ! -f "${R1}" ]]; then
    echo "[ERROR] Missing trimmed R1 FASTQ: ${R1}" >&2
    exit 1
  fi

  if [[ ! -f "${R2}" ]]; then
    echo "[ERROR] Missing trimmed R2 FASTQ: ${R2}" >&2
    exit 1
  fi

  salmon quant \
    -i "${INDEX_DIR}" \
    -l A \
    -1 "${R1}" \
    -2 "${R2}" \
    --validateMappings \
    --gcBias \
    --seqBias \
    -p "${THREADS}" \
    -o "${OUT_DIR}/${sample_id}"
done

echo "[INFO] Salmon quantification completed."
