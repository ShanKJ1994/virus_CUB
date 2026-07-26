process_file() {
  local i="$1"
  sample_id=${i##*/}
  sample_id=${sample_id%_Aligned.sortedByCoord.out.bam}
  prefix="/gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_Hs_dedup"
  log_file=${prefix}/${sample_id}.dedup.log
  {
    umi_tools dedup --stdin=$i \
      --umi-separator=: --paired \
      --output-stats=${prefix}/${sample_id} \
      --stdout=${prefix}/${sample_id}.dedup.bam
    echo "$(date) ${sample_id} bam file deduplicated"
    echo ""
  } >${log_file} 2>&1
}

echo "$(date) task started"
export -f process_file

ls /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/alignment_Hs/*_Aligned.sortedByCoord.out.bam | xargs -P 8 -I {} bash -c 'process_file "$@"' _ {}

echo "$(date) all finished"
