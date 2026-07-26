for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_length_filtered/*.length_filtered.fastq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id="${sample_id%.length_filtered.fastq.gz}-plasmid"
  out_prefix="/gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_RIBO_plasmid/${sample_id}_"
  STAR --runThreadN 32 --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
   --genomeDir /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/spike_index/index \
   --readFilesCommand zcat \
   --alignEndsType EndToEnd \
   --readFilesIn $i --outFileNamePrefix $out_prefix \
   --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

STAR --genomeLoad Remove --genomeDir /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/spike_index/index \
  --outFileNamePrefix remove_plasmid_genome_output
rm remove_plasmid_genome_output*

for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_RIBO_plasmid/*_Aligned.sortedByCoord.out.bam
do
  samtools index $i
done