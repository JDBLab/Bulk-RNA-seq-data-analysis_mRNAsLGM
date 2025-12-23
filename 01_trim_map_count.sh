#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Bulk RNA-seq standard pipeline:
#   Merge lanes → Trim Galore → Bowtie2 → sort/index BAM → featureCounts → MultiQC
# ============================================================

# ---------- Read config.yaml (minimal YAML parser) ----------
CFG="config/config.yaml"

# Simple YAML getter for "key: value" lines (no quotes required; supports quoted)
yget () {
  local key="$1"
  python - <<'PY' "$CFG" "$key"
import sys, re
cfg, key = sys.argv[1], sys.argv[2]
val = None
with open(cfg, "r", encoding="utf-8") as f:
  for line in f:
    line=line.strip()
    if not line or line.startswith("#"): 
      continue
    m = re.match(rf"^{re.escape(key)}\s*:\s*(.+)$", line)
    if m:
      v = m.group(1).strip()
      # strip surrounding quotes
      if (v.startswith('"') and v.endswith('"')) or (v.startswith("'") and v.endswith("'")):
        v = v[1:-1]
      val = v
      break
print("" if val is None else val)
PY
}

FASTQ_DIR="$(yget fastq_dir)"
BOWTIE2_INDEX_PREFIX="$(yget bowtie2_index_prefix)"
GTF="$(yget gtf)"
OUTDIR="$(yget outdir)"
THREADS="$(yget threads)"
STRAND="$(yget strand)"
TRIM_Q="$(yget trim_galore_quality)"
TRIM_MINLEN="$(yget trim_galore_minlen)"

# ---------- Sanity checks ----------
if [[ -z "${FASTQ_DIR}" || -z "${BOWTIE2_INDEX_PREFIX}" || -z "${GTF}" ]]; then
  echo "ERROR: Please edit config/config.yaml and set fastq_dir, bowtie2_index_prefix, gtf."
  exit 1
fi

mkdir -p "${OUTDIR}"/{00_versions,01_merged_fastq,02_trimmed,03_bam,04_counts,05_qc}

# ---------- Record tool versions ----------
{
  echo "=== trim_galore ==="; trim_galore --version 2>/dev/null || true
  echo "=== cutadapt ==="; cutadapt --version 2>/dev/null || true
  echo "=== fastqc ==="; fastqc --version 2>/dev/null || true
  echo "=== bowtie2 ==="; bowtie2 --version 2>/dev/null | head -n 1 || true
  echo "=== samtools ==="; samtools --version 2>/dev/null | head -n 2 || true
  echo "=== featureCounts ==="; featureCounts -v 2>/dev/null || true
  echo "=== multiqc ==="; multiqc --version 2>/dev/null || true
} > "${OUTDIR}/00_versions/tool_versions.txt"

# ---------- Detect sample IDs ----------
# Sample ID is everything before the first underscore.
# Example: 1CI_S61_L001_R1_001.fastq.gz -> 1CI
mapfile -t SAMPLES < <(ls -1 "${FASTQ_DIR}"/*_R1_001.fastq.gz \
  | xargs -n1 basename \
  | sed -E 's/_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz$//' \
  | sort -u)

echo "Found samples:"
printf "  %s\n" "${SAMPLES[@]}"
echo

# ---------- Merge lanes ----------
echo "Merging lanes per sample..."
for s in "${SAMPLES[@]}"; do
  r1_glob="${FASTQ_DIR}/${s}"_S*_L*_R1_001.fastq.gz
  r2_glob="${FASTQ_DIR}/${s}"_S*_L*_R2_001.fastq.gz

  if ! ls $r1_glob >/dev/null 2>&1; then
    echo "WARNING: no R1 files found for ${s}. Skipping."
    continue
  fi

  echo "  ${s}"
  cat $r1_glob > "${OUTDIR}/01_merged_fastq/${s}_R1.fastq.gz"
  cat $r2_glob > "${OUTDIR}/01_merged_fastq/${s}_R2.fastq.gz"
done

# ---------- Trim Galore ----------
echo
echo "Running Trim Galore..."
for s in "${SAMPLES[@]}"; do
  in1="${OUTDIR}/01_merged_fastq/${s}_R1.fastq.gz"
  in2="${OUTDIR}/01_merged_fastq/${s}_R2.fastq.gz"
  [[ -f "$in1" && -f "$in2" ]] || { echo "Skipping ${s} (missing merged FASTQs)"; continue; }

  trim_galore \
    --paired \
    --quality "${TRIM_Q}" \
    --length "${TRIM_MINLEN}" \
    --stringency 3 \
    --fastqc \
    --cores 4 \
    --output_dir "${OUTDIR}/02_trimmed" \
    "$in1" "$in2"
done

# ---------- Bowtie2 mapping ----------
echo
echo "Mapping with Bowtie2..."
for s in "${SAMPLES[@]}"; do
  t1="${OUTDIR}/02_trimmed/${s}_R1_val_1.fq.gz"
  t2="${OUTDIR}/02_trimmed/${s}_R2_val_2.fq.gz"
  [[ -f "$t1" && -f "$t2" ]] || { echo "Skipping ${s} (missing trimmed FASTQs)"; continue; }

  bam="${OUTDIR}/03_bam/${s}.bam"
  log="${OUTDIR}/03_bam/${s}.bowtie2.log"

  bowtie2 \
    --very-sensitive \
    -x "${BOWTIE2_INDEX_PREFIX}" \
    -1 "$t1" -2 "$t2" \
    -p "${THREADS}" \
    2> "$log" \
  | samtools view -bS - \
  | samtools sort -@ "${THREADS}" -o "$bam"

  samtools index "$bam"
  samtools flagstat "$bam" > "${OUTDIR}/05_qc/${s}.flagstat.txt"
done

# ---------- featureCounts ----------
echo
echo "Running featureCounts..."
shopt -s nullglob
BAMS=("${OUTDIR}/03_bam/"*.bam)
if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "ERROR: No BAM files found in ${OUTDIR}/03_bam"
  exit 1
fi

featureCounts \
  -T "${THREADS}" \
  -a "${GTF}" \
  -o "${OUTDIR}/04_counts/featureCounts.txt" \
  -g gene_id \
  -t exon \
  -s "${STRAND}" \
  -p -B -C \
  "${BAMS[@]}"

# ---------- MultiQC ----------
echo
echo "Running MultiQC..."
multiqc "${OUTDIR}" -o "${OUTDIR}/05_qc" >/dev/null 2>&1 || true

echo
echo "DONE."
echo "Counts: ${OUTDIR}/04_counts/featureCounts.txt"
echo "QC:     ${OUTDIR}/05_qc/"
