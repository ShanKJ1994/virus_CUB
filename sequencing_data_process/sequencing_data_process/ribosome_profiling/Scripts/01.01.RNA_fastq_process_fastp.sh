#!/bin/bash
for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RNA_raw/*_R1.fq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id=${sample_id%_R1.fq.gz}
  fastp -l 30 -g --overlap_len_require 20 -y -Y 10 \
   -o /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RNA_processed_fastp/${sample_id}_R1.processed.fastq.gz \
   -O /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RNA_processed_fastp/${sample_id}_R2.processed.fastq.gz \
   -i $i -I ${i/_R1.fq.gz/_R2.fq.gz}
  echo ${sample_id}' completed adapter trimming'
  echo ''
done




