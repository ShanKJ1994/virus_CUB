for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_length_filtered/*.length_filtered.fastq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id=${sample_id%.length_filtered.fastq.gz}
  out_prefix="/gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_RIBO_Hs/${sample_id}_"
  STAR --runThreadN 32 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
   --genomeDir /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/index \
   --readFilesCommand zcat \
   --alignEndsType EndToEnd \
   --readFilesIn $i --outFileNamePrefix $out_prefix \
   --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

STAR --genomeLoad Remove --genomeDir /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/index \
  --outFileNamePrefix remove_genome_output
rm remove_genome_output*

for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_RIBO_Hs/*.bam
do
  samtools index $i
done