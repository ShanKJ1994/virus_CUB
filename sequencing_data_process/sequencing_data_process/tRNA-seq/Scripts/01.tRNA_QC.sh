#!/bin/bash

bowtie2-build /data/fanyujie/project/tRNA/Homo_rna_index/rrna2.fasta /data/fanyujie/project/tRNA/Homo_rna_index/rrna2

RAW_DIR="/data/fanyujie/project/tRNA/raw_data/"
OUT_DIR="/data/fanyujie/project/tRNA/clean_data/"
IDX_DIR="/data/fanyujie/project/tRNA/Homo_rna_index/"
ADAPTER="GATATCGTCAAGATCGGAAGAGCACACGTCTGAA"
THREADS=16

mkdir -p "$OUT_DIR"

echo "Pipeline started at $(date)"
echo "Adapter: $ADAPTER"
echo "=========================================="

for fq in ${RAW_DIR}*.fq.gz; do
    base=$(basename "$fq" .fq.gz)
    out="${OUT_DIR}${base}/"
    mkdir -p "$out"
    log="${out}${base}.log"

    echo "[$(date)] Processing $base" | tee "$log"

    cutadapt -a "$ADAPTER" --discard-untrimmed \
        -o "${out}step1.fq" "$fq" >> "$log" 2>&1 || continue

    cutadapt -m 5 -u -2 -j "$THREADS" \
        -o "${out}step2.fq" "${out}step1.fq" >> "$log" 2>&1 || continue

    cutadapt -m 14 -u 4 -j "$THREADS" \
        -o "${out}step3.fq" "${out}step2.fq" >> "$log" 2>&1 || continue

    cutadapt -m 20 \
        -o "${out}step4.fq" "${out}step3.fq" >> "$log" 2>&1 || continue

    bowtie2 -p "$THREADS" --local --very-sensitive-local \
        -x "$IDX_DIR" \
        -U "${out}step4.fq" \
        --un-gz "${out}${base}.fq.gz" \
        -S "${out}rRNA.sam" >> "$log" 2>&1

    if [ -f "${out}${base}.fq.gz" ]; then
        echo "[SUCCESS] $base" | tee -a "$log"
    else
        echo "[FAILED] $base" | tee -a "$log"
    fi

    rm -f "${out}"step*.fq "${out}"rRNA.sam
    echo "----------------------------------------" >> "$log"
done

echo "=========================================="
echo "Done at $(date)"
