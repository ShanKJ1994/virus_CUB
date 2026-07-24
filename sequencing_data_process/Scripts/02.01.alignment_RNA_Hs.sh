#! /bin/bash

#align to genome

for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RNA_processed_fastp/*_R1.processed.fastq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id="${sample_id%_R1.processed.fastq.gz}-Hs"
	out_prefix="/gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_Hs/${sample_id}_"
	STAR --runThreadN 32 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
	 --genomeDir /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/index --readFilesCommand zcat \
	 --readFilesIn $i ${i/_R1.processed.fastq.gz/_R2.processed.fastq.gz} --outFileNamePrefix $out_prefix \
	 --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

STAR --genomeLoad Remove --genomeDir /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/index --outFileNamePrefix tmp_genome_remove
rm tmp_genome_remove*

for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_Hs/*_Aligned.sortedByCoord.out.bam
do
  samtools index $i
done







