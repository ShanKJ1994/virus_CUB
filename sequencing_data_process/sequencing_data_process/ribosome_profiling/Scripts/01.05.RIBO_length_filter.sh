for i in /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_rRNArm/*.rRNArm.fastq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id=${sample_id%.rRNArm.fastq.gz}
  seqkit seq -m 26 -M 33 $i \
    -o /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_length_filtered/${sample_id}.length_filtered.fastq.gz
done

fastqc -o /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_length_filtered/fastqc /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_length_filtered/*.length_filtered.fastq.gz