nohup featureCounts -T 16 -s 2 -p -B -O -M --fraction --countReadPairs --fracOverlap 0.5 \
 -a /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/gtf/Homo_sapiens.GRCh38.104.gtf \
 -t CDS -g gene_id -o featureCounts_CodonOptimizeXBB_Hs_dedup_bothMapped_rev.txt \
 /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_Hs_dedup/*.dedup.bam > featureCounts_CodonOptimizeXBB_Hs_dedup_bothMapped_rev.log 2>&1 &

nohup featureCounts -T 16 -s 2 -p -B -O -M --countReadPairs --fracOverlap 0.5 \
 -a /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/spike_index/XBB.1.5_s1to5.gtf \
 -t CDS -g gene_id -o featureCounts_CodonOptimizeXBB_plasmid_dedup_bothMapped_rev.txt \
 /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_plasmid_dedup/*.dedup.bam > featureCounts_CodonOptimizeXBB_plasmid_dedup_bothMapped_rev.log 2>&1 &
 
nohup featureCounts -T 16 -s 1 -O -M --fraction -t CDS -g gene_id --fracOverlap 0.5 \
  -a /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/gtf/Homo_sapiens.GRCh38.104.gtf \
  -o featureCounts_CodonOptimizeXBB_RIBO_Hs_dedup.txt \
  /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_RIBO_Hs_dedup/*.dedup.bam \
  > featureCounts_CodonOptimizeXBB_RIBO_Hs_dedup.log 2>&1 &

nohup featureCounts -T 16 -s 1 -O -M -t CDS -g gene_id --fracOverlap 0.5 \
  -a /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/spike_index/XBB.1.5_s1to5.gtf \
  -o featureCounts_CodonOptimizeXBB_RIBO_plasmid_dedup.txt \
  /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_RIBO_plasmid_dedup/*.dedup.bam \
  > featureCounts_CodonOptimizeXBB_RIBO_plasmid_dedup.log 2>&1 &
