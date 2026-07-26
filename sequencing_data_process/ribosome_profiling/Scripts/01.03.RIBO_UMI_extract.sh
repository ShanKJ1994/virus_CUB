process_file() {
  local i="$1"
  sample_id=${i##*/}
  sample_id=${sample_id%_R1.fastp_corrected.fastq.gz}
  prefix="/gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_UMI_extracted"
  log_file=${prefix}/${sample_id}.umi_extract.log
  {
    echo "$(date) ${sample_id} processing started"
    echo ""
    umi_tools extract --stdin=$i \
      --bc-pattern="^(?P<umi_1>.{3})(?P<cell_1>(ATC|CGA)){s<=1}(?P<umi_2>.{3})(?P<cell_2>(AGT|TAC)){s<=1}(?P<umi_3>.{3}).*(?P<umi_4>.{3})(?P<cell_3>GTC){s<=1}(?P<umi_5>.{3})(?P<cell_4>TAG){s<=1}(?P<umi_6>.{3})$" \
      --stdout=${prefix}/${sample_id}.umi_extracted.fastq.gz \
      --extract-method=regex --quality-filter-threshold=10 --quality-encoding=phred33
    echo "$(date) ${sample_id} umi extracted"
    echo ""
  } >${log_file} 2>&1
}

echo "$(date) task started"
export -f process_file

ls /gpfs2/gaopx/CodonOptimum/XBB.1.5_S_RiboSeq/fastq/RIBO_fastp_corrected/*_R1.fastp_corrected.fastq.gz | xargs -P 16 -I {} bash -c 'process_file "$@"' _ {}

echo "$(date) all finished"
