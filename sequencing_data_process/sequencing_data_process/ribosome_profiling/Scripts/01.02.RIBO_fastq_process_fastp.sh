#!/bin/bash
for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_raw/*_R1.fq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id=${sample_id%_R1.fq.gz}
  fastp -c -Q --length_required 45 --length_limit 75 \
   --overlap_len_require 30 --overlap_diff_limit 50 --overlap_diff_percent_limit 20 \
   -o /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_fastp_corrected/${sample_id}_R1.fastp_corrected.fastq.gz \
   -O /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_fastp_corrected/${sample_id}_R2.fastp_corrected.fastq.gz \
   -j /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_fastp_corrected/${sample_id}_fastp.json \
   -h /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_fastp_corrected/${sample_id}_fastp.html \
   -i $i -I ${i/_R1.fq.gz/_R2.fq.gz}
  echo ${sample_id}' completed adapter trimming and correction'
  echo ''
done