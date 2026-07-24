print(paste0(Sys.time(),' process started'))

library(riboWaltz)
library(biomaRt)
library(ggplot2)
library(data.table)
library(rtracklayer)
library(dplyr)

print(paste0(Sys.time(),' packages loaded'))

gtf_path = '/lustre/user/lulab/gaopx/genome/GRCh38/ensemble/gtf/Homo_sapiens.GRCh38.104.gtf'
bam_folder = '/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_host/'
bam_suffix = '-host_Aligned.toTranscriptome.out.bam'
output_path = '/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/psite_offset.csv'
pic_path = '/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/pic'

longest_cds_path = '/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/longest_cds.csv'
annot_cds_path = '/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/annot_longest_cds.csv'

if (!dir.exists(pic_path)) dir.create(pic_path, recursive = TRUE)

gtf = import(gtf_path)
gtf_CDS = gtf[gtf$type == "CDS", ]
gtf_exon = gtf[gtf$type == "exon", ]
gtf_CDS$length = width(gtf_CDS)
gtf_exon$length = width(gtf_exon)
CDS_summary = mcols(gtf_CDS) %>% as_tibble() %>% group_by(transcript_id, gene_id) %>%
  summarise(cds_length = sum(length)) %>% arrange(transcript_id)
transcript_summary = mcols(gtf_exon) %>% as_tibble() %>% group_by(transcript_id, gene_id) %>%
  summarise(transcript_length = sum(length)) %>% arrange(transcript_id)
CDS_summary$transcript_length = transcript_summary$transcript_length[match(CDS_summary$transcript_id, transcript_summary$transcript_id)]
CDS_longest = CDS_summary %>% arrange(gene_id, desc(cds_length), desc(transcript_length)) %>%
  group_by(gene_id) %>% slice_head(n = 1) %>% ungroup()

write.csv(CDS_longest, longest_cds_path, row.names = F, quote = F)

annot = create_annotation(gtf_path)
annot = annot[which(annot$transcript%in%CDS_longest$transcript_id),]

write.csv(annot, annot_cds_path, row.names = F, quote = F)

print(paste0('Effective genes num with CDS: ', nrow(annot)))

bamfile_names = list.files(bam_folder,paste0(bam_suffix, '$'))
names(bamfile_names) = as.vector(sapply(bamfile_names, function(x){return(strsplit(x,'.bam')[[1]][1])}))
bamfile_names = sapply(bamfile_names, function(x){return(strsplit(x,bam_suffix)[[1]][1])})
bam_list = bamtolist(bamfolder = bam_folder, annotation = annot, name_samples = bamfile_names)
bam_list = length_filter(data = bam_list, length_filter_mode = "custom", length_range = 20:40)

psite_offset = psite(bam_list, extremity = "auto")
write.csv(psite_offset, output_path, row.names = F, quote = F)


reads_psite_list = psite_info(bam_list, psite_offset)
length_plot = rlength_distr(reads_psite_list,
                            sample = bamfile_names,
                            multisamples = 'independent',
                            transcripts = annot$transcript,
                            colour = "#333f50")

frames_plot = frame_psite(reads_psite_list, annot,
                               sample = bamfile_names,
                               transcripts = annot$transcript,
                               multisamples = "independent",
                               plot_style = "split",
                               region = "all",
                               colour = "#333f50")

metaprofile_plot = metaprofile_psite(data = reads_psite_list,
                                        annotation = annot,
                                        sample = bamfile_names,
                                        multisamples = "independent",
                                        transcripts = annot$transcript,
                                        utr5l = 30, cdsl = 100, utr3l = 30,
                                        colour = "#333f50")
frame_psite_plot = frame_psite_length(data = reads_psite_list,
                                      annotation = annot,
                                      sample = bamfile_names,
                                      transcripts = annot$transcript,
                                      multisamples = "independent",
                                      plot_style = "split",
                                      region = "all",
                                      colour = "#333f50")
for(i in bamfile_names){
  ggsave(filename=paste0(pic_path,'/',i,'.length.png'), plot=length_plot[[paste0('plot_', i)]], device='png')
  ggsave(filename=paste0(pic_path,'/',i,'.frame.png'), plot=frames_plot[[paste0('plot_', i)]], device='png')
  ggsave(filename=paste0(pic_path,'/',i,'.metaprofile.png'), plot=metaprofile_plot[[paste0('plot_', i)]], device='png')
  ggsave(filename=paste0(pic_path,'/',i,'.frame_by_length.png'), plot=frame_psite_plot[[paste0('plot_', i)]], device='png')

}

