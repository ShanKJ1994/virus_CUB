#! /bin/bash

#align to genome

for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RNA_processed_fastp/*_R1.processed.fastq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id="${sample_id%_R1.processed.fastq.gz}-plasmid"
	out_prefix="/gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_plasmid/${sample_id}_"
	STAR --runThreadN 32 --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
	 --genomeDir /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/spike_index/index --readFilesCommand zcat \
	 --readFilesIn $i ${i/_R1.processed.fastq.gz/_R2.processed.fastq.gz} --outFileNamePrefix $out_prefix \
	 --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

STAR --genomeLoad Remove --genomeDir /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/spike_index/index --outFileNamePrefix tmp_genome_remove
rm tmp_genome_remove*

for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_plasmid/*_Aligned.sortedByCoord.out.bam
do
  samtools index $i
done







