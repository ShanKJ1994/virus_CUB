for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_UMI_extracted/*.umi_extracted.fastq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id=${sample_id%.umi_extracted.fastq.gz}
  bowtie2 -p 16 --no-unal --local --un-gz /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_rRNArm/${sample_id}.rRNArm.fastq.gz \
    -x /lustre/user/lulab/gaopx/genome/Hs/rRNA_NCBI_tRNA_GtRNAdb/index_bt2/rRNA_NCBI_tRNA_GtRNAdb \
    -U $i | samtools sort -O BAM -@ 16 -o /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_rRNArm/bt2_res/${sample_id}_rRNA_tRNA.bam
  echo ${sample_id}' completed rRNA removing'
  echo ''
  rm $i
done

fastqc -o /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_rRNArm/fastqc /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_rRNArm/*.rRNArm.fastq.gz \
  > /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_rRNArm/fastqc/fastqc.log