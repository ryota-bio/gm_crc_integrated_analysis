#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
META="${PROJECT_DIR}/metadata/sample_metadata.csv"
RAW_DIR="${PROJECT_DIR}/data/raw_fastq"
TRIM_DIR="${PROJECT_DIR}/results/fastp_trimmed"
QC_DIR="${PROJECT_DIR}/results/fastp_qc"

THREADS=${THREADS:-4}

mkdir -p "${TRIM_DIR}" "${QC_DIR}"

echo "[INFO] Project directory: ${PROJECT_DIR}"
echo "[INFO] Metadata file: ${META}"

tail -n +2 "${META}" | while IFS=',' read -r sample_id condition time_point replicate cell_line bacteria MOI infection_time fastq_R1 fastq_R2
do
  echo "[fastp] Processing ${sample_id}"

  if [[ ! -f "${RAW_DIR}/${fastq_R1}" ]]; then
    echo "[ERROR] Missing R1 FASTQ: ${RAW_DIR}/${fastq_R1}" >&2
    exit 1
  fi

  if [[ ! -f "${RAW_DIR}/${fastq_R2}" ]]; then
    echo "[ERROR] Missing R2 FASTQ: ${RAW_DIR}/${fastq_R2}" >&2
    exit 1
  fi

  fastp \
    -i "${RAW_DIR}/${fastq_R1}" \
    -I "${RAW_DIR}/${fastq_R2}" \
    -o "${TRIM_DIR}/${sample_id}_R1.trimmed.fastq.gz" \
    -O "${TRIM_DIR}/${sample_id}_R2.trimmed.fastq.gz" \
    --html "${QC_DIR}/${sample_id}.fastp.html" \
    --json "${QC_DIR}/${sample_id}.fastp.json" \
    --thread "${THREADS}"
done

echo "[INFO] fastp QC completed."
