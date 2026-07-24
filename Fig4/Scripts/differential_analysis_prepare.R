library(DESeq2)
library(xtail)
library(tidyverse)

CodonOptimizeXBB_featureCounts_RNA_all = read.csv('Fig4/Results/featureCounts_CodonOptimizeXBB_RNA_combined_dedup_bothMapped_rev.csv')
sample_info = read.csv('Fig4/Results/sample_info.csv')
CodonOptimizeXBB_featureCounts_RIBO_all = read.csv('Fig4/Results/featureCounts_CodonOptimizeXBB_RIBO_combined_dedup.csv')

u_Hs_104 = read.csv('Fig4/Results/Biomart_human_annotation_ensembl_104.csv')
u_Hs_canonical_104 = u_Hs_104[which(u_Hs_104$transcript_is_canonical==1),]
u_Hs_canonical_chr_104 = u_Hs_canonical_104[which(u_Hs_canonical_104$chromosome_name%in%c(as.character(1:22), 'X', 'Y')),]
u_Hs_canonical_chr_procod_104 = u_Hs_canonical_chr_104[which(u_Hs_canonical_chr_104$gene_biotype=='protein_coding'),]
u_Hs_canonical_104_MT = u_Hs_canonical_104[which(u_Hs_canonical_104$chromosome_name=='MT'),]

CodonOptimizeXBB_featureCounts_RNA_all$gene_name = u_Hs_canonical_104$external_gene_name[match(CodonOptimizeXBB_featureCounts_RNA_all$Geneid, u_Hs_canonical_104$ensembl_gene_id)]
CodonOptimizeXBB_featureCounts_RNA = CodonOptimizeXBB_featureCounts_RNA_all
libsize_RNA = colSums(CodonOptimizeXBB_featureCounts_RNA[,3:10])
CodonOptimizeXBB_featureCounts_RNA_RPKM = CodonOptimizeXBB_featureCounts_RNA
CodonOptimizeXBB_featureCounts_RNA_RPKM[,3:10] = t(apply(CodonOptimizeXBB_featureCounts_RNA_RPKM[,3:10], 1, function(x) 1000000*x/libsize_RNA))
CodonOptimizeXBB_featureCounts_RNA_RPKM[,3:10] = apply(CodonOptimizeXBB_featureCounts_RNA_RPKM[,3:10], 2, function(x) 1000*x/CodonOptimizeXBB_featureCounts_RNA_RPKM$Length)

CodonOptimizeXBB_featureCounts_RIBO_all$gene_name = u_Hs_canonical_104$external_gene_name[match(CodonOptimizeXBB_featureCounts_RIBO_all$Geneid, u_Hs_canonical_104$ensembl_gene_id)]
CodonOptimizeXBB_featureCounts_RIBO = CodonOptimizeXBB_featureCounts_RIBO_all
libsize_RIBO = colSums(CodonOptimizeXBB_featureCounts_RIBO[,3:10])
CodonOptimizeXBB_featureCounts_RIBO_RPKM = CodonOptimizeXBB_featureCounts_RIBO
CodonOptimizeXBB_featureCounts_RIBO_RPKM[,3:10] = t(apply(CodonOptimizeXBB_featureCounts_RIBO_RPKM[,3:10], 1, function(x) 1000000*x/libsize_RIBO))
CodonOptimizeXBB_featureCounts_RIBO_RPKM[,3:10] = apply(CodonOptimizeXBB_featureCounts_RIBO_RPKM[,3:10], 2, function(x) 1000*x/CodonOptimizeXBB_featureCounts_RIBO_RPKM$Length)

CodonOptimizeXBB_featureCounts_combined_all = cbind(CodonOptimizeXBB_featureCounts_RNA_all[,c(1,2,11,3:10)], CodonOptimizeXBB_featureCounts_RIBO_all[match(CodonOptimizeXBB_featureCounts_RNA_all$Geneid, CodonOptimizeXBB_featureCounts_RIBO_all$Geneid),3:10])
CodonOptimizeXBB_featureCounts_combined = CodonOptimizeXBB_featureCounts_combined_all[which(CodonOptimizeXBB_featureCounts_combined_all$Geneid%in%c(u_Hs_canonical_chr_procod_104$ensembl_gene_id, 'spike', 'NeoR')),]

libsize_combined = colSums(CodonOptimizeXBB_featureCounts_combined[,4:19])
CodonOptimizeXBB_featureCounts_combined_RPKM = CodonOptimizeXBB_featureCounts_combined
CodonOptimizeXBB_featureCounts_combined_RPKM[,4:19] = t(apply(CodonOptimizeXBB_featureCounts_combined_RPKM[,4:19], 1, function(x) 1000000*x/libsize_combined))
CodonOptimizeXBB_featureCounts_combined_RPKM[,4:19] = apply(CodonOptimizeXBB_featureCounts_combined_RPKM[,4:19], 2, function(x) 1000*x/CodonOptimizeXBB_featureCounts_combined_RPKM$Length)

CodonOptimizeXBB_featureCounts_combined_TE = CodonOptimizeXBB_featureCounts_combined_RPKM
treatment_types = c('Ctrl_1','Ctrl_2','XS1_1','XS1_2','XS2_1','XS2_2','XS4_1','XS4_2')
for(i in treatment_types){
  CodonOptimizeXBB_featureCounts_combined_TE[,paste0('TE_', i)] = CodonOptimizeXBB_featureCounts_combined_TE[,paste0('mono_', i)]/CodonOptimizeXBB_featureCounts_combined_TE[,paste0('mR_', i)]
}

CodonOptimizeXBB_featureCounts_combined_RPKM_filterRNACtrlRPKM1 = CodonOptimizeXBB_featureCounts_combined_RPKM[which(rowMeans(CodonOptimizeXBB_featureCounts_combined_RPKM[,4:5])>1),]
CodonOptimizeXBB_featureCounts_combined_count_filterRNACtrlRPKM1 = CodonOptimizeXBB_featureCounts_combined[which(CodonOptimizeXBB_featureCounts_combined$Geneid%in%CodonOptimizeXBB_featureCounts_combined_RPKM_filterRNACtrlRPKM1$Geneid),]
rownames(CodonOptimizeXBB_featureCounts_combined_count_filterRNACtrlRPKM1) = CodonOptimizeXBB_featureCounts_combined_count_filterRNACtrlRPKM1$Geneid

# # Generate the RPKM_TE table.
# mRNA_RPF_TE_RPKM = CodonOptimizeXBB_featureCounts_combined_TE[which(CodonOptimizeXBB_featureCounts_combined_TE$Geneid%in%c(CodonOptimizeXBB_featureCounts_combined_RPKM_filterRNACtrlRPKM1$Geneid,'spike')),]
# colnames(mRNA_RPF_TE_RPKM) = gsub('mR','mRNA',colnames(mRNA_RPF_TE_RPKM))
# colnames(mRNA_RPF_TE_RPKM) = gsub('mono','Ribo',colnames(mRNA_RPF_TE_RPKM))
# mRNA_RPF_TE_RPKM = mRNA_RPF_TE_RPKM %>% relocate(gene_name, .after = last_col()) %>% rename(external_gene_name=gene_name) %>% arrange(Geneid)
# write.csv(mRNA_RPF_TE_RPKM,'Fig4/Results/mRNA_RPF_TE_RPKM.csv', row.names = F, quote = F)





CodonOptimizeXBB_featureCounts_RNA_filter = CodonOptimizeXBB_featureCounts_combined_count_filterRNACtrlRPKM1[,3:11]
CodonOptimizeXBB_featureCounts_RIBO_filter = CodonOptimizeXBB_featureCounts_combined_count_filterRNACtrlRPKM1[,c(3,12:19)]

sample_info_RNA = sample_info[which(sample_info$Type=='RNA'),]
dds_CodonOptimizeXBB_RNA = DESeqDataSetFromMatrix(countData = floor(CodonOptimizeXBB_featureCounts_RNA_filter[,sample_info_RNA$Name]),
                                                  colData = sample_info_RNA,
                                                  design = ~ treatment)
dds_CodonOptimizeXBB_RNA = DESeq(dds_CodonOptimizeXBB_RNA)
RNAFC_XS1_vs_Ctrl = data.frame(results(dds_CodonOptimizeXBB_RNA, contrast = c('treatment','XS1','Ctrl')))
RNAFC_XS1_vs_Ctrl$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RNAFC_XS1_vs_Ctrl), u_Hs_canonical_104$ensembl_gene_id)]
RNAFC_XS2_vs_Ctrl = data.frame(results(dds_CodonOptimizeXBB_RNA, contrast = c('treatment','XS2','Ctrl')))
RNAFC_XS2_vs_Ctrl$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RNAFC_XS2_vs_Ctrl), u_Hs_canonical_104$ensembl_gene_id)]
RNAFC_XS4_vs_Ctrl = data.frame(results(dds_CodonOptimizeXBB_RNA, contrast = c('treatment','XS4','Ctrl')))
RNAFC_XS4_vs_Ctrl$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RNAFC_XS4_vs_Ctrl), u_Hs_canonical_104$ensembl_gene_id)]
RNAFC_XS4_vs_XS1 = data.frame(results(dds_CodonOptimizeXBB_RNA, contrast = c('treatment','XS4','XS1')))
RNAFC_XS4_vs_XS1$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RNAFC_XS4_vs_XS1), u_Hs_canonical_104$ensembl_gene_id)]
RNAFC_XS4_vs_XS2 = data.frame(results(dds_CodonOptimizeXBB_RNA, contrast = c('treatment','XS4','XS2')))
RNAFC_XS4_vs_XS2$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RNAFC_XS4_vs_XS2), u_Hs_canonical_104$ensembl_gene_id)]
RNAFC_XS2_vs_XS1 = data.frame(results(dds_CodonOptimizeXBB_RNA, contrast = c('treatment','XS2','XS1')))
RNAFC_XS2_vs_XS1$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RNAFC_XS2_vs_XS1), u_Hs_canonical_104$ensembl_gene_id)]

sample_info_RIBO = sample_info[which(sample_info$Type=='RIBO'),]
dds_CodonOptimizeXBB_RIBO = DESeqDataSetFromMatrix(countData = floor(CodonOptimizeXBB_featureCounts_RIBO_filter[,sample_info_RIBO$Name]),
                                                   colData = sample_info_RIBO,
                                                   design = ~ treatment)
dds_CodonOptimizeXBB_RIBO = DESeq(dds_CodonOptimizeXBB_RIBO)
RIBOFC_XS1_vs_Ctrl = data.frame(results(dds_CodonOptimizeXBB_RIBO, contrast = c('treatment','XS1','Ctrl')))
RIBOFC_XS1_vs_Ctrl$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RIBOFC_XS1_vs_Ctrl), u_Hs_canonical_104$ensembl_gene_id)]
RIBOFC_XS2_vs_Ctrl = data.frame(results(dds_CodonOptimizeXBB_RIBO, contrast = c('treatment','XS2','Ctrl')))
RIBOFC_XS2_vs_Ctrl$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RIBOFC_XS2_vs_Ctrl), u_Hs_canonical_104$ensembl_gene_id)]
RIBOFC_XS4_vs_Ctrl = data.frame(results(dds_CodonOptimizeXBB_RIBO, contrast = c('treatment','XS4','Ctrl')))
RIBOFC_XS4_vs_Ctrl$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RIBOFC_XS4_vs_Ctrl), u_Hs_canonical_104$ensembl_gene_id)]
RIBOFC_XS4_vs_XS1 = data.frame(results(dds_CodonOptimizeXBB_RIBO, contrast = c('treatment','XS4','XS1')))
RIBOFC_XS4_vs_XS1$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RIBOFC_XS4_vs_XS1), u_Hs_canonical_104$ensembl_gene_id)]
RIBOFC_XS4_vs_XS2 = data.frame(results(dds_CodonOptimizeXBB_RIBO, contrast = c('treatment','XS4','XS2')))
RIBOFC_XS4_vs_XS2$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RIBOFC_XS4_vs_XS2), u_Hs_canonical_104$ensembl_gene_id)]
RIBOFC_XS2_vs_XS1 = data.frame(results(dds_CodonOptimizeXBB_RIBO, contrast = c('treatment','XS2','XS1')))
RIBOFC_XS2_vs_XS1$gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(RIBOFC_XS2_vs_XS1), u_Hs_canonical_104$ensembl_gene_id)]

CodonOptimizeXBB_featureCounts_combined_for_xtail = CodonOptimizeXBB_featureCounts_combined_count_filterRNACtrlRPKM1
for(t in c('XS1','XS2','XS4')){
  list_RNA = sample_info[sample_info$Type == 'RNA'&sample_info$treatment%in%c('Ctrl', t), c('Name','treatment')]
  list_RIBO = sample_info[sample_info$Type == 'RIBO'&sample_info$treatment%in%c('Ctrl', t), c('Name','treatment')]
  count_table_RNA = subset(CodonOptimizeXBB_featureCounts_combined_for_xtail, select = as.vector(list_RNA$Name)) %>% floor()
  count_table_RIBO = subset(CodonOptimizeXBB_featureCounts_combined_for_xtail, select = as.vector(list_RIBO$Name)) %>% floor()
  if(all(list_RNA$treatment == list_RIBO$treatment)) condition = as.vector(list_RNA$treatment) else print('WRONG!!!')
  res_CodonOptimizeXBB_RiboSeq = xtail(count_table_RNA, count_table_RIBO, condition, baseLevel = 'Ctrl')
  assign(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsCtrl'), data.frame(resultsTable(res_CodonOptimizeXBB_RiboSeq, log2FCs = T, log2Rs = T)))
  gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(get(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsCtrl'))),u_Hs_canonical_104$ensembl_gene_id)]
  assign(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsCtrl'), cbind(get(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsCtrl')),gene_name))
}
for(t in c('XS2','XS4')){
  list_RNA = sample_info[sample_info$Type == 'RNA'&sample_info$treatment%in%c('XS1', t), c('Name','treatment')]
  list_RIBO = sample_info[sample_info$Type == 'RIBO'&sample_info$treatment%in%c('XS1', t), c('Name','treatment')]
  count_table_RNA = subset(CodonOptimizeXBB_featureCounts_combined_for_xtail, select = as.vector(list_RNA$Name)) %>% floor()
  count_table_RIBO = subset(CodonOptimizeXBB_featureCounts_combined_for_xtail, select = as.vector(list_RIBO$Name)) %>% floor()
  if(all(list_RNA$treatment == list_RIBO$treatment)) condition = as.vector(list_RNA$treatment) else print('WRONG!!!')
  res_CodonOptimizeXBB_RiboSeq = xtail(count_table_RNA, count_table_RIBO, condition, baseLevel = 'XS1')
  assign(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS1'), data.frame(resultsTable(res_CodonOptimizeXBB_RiboSeq, log2FCs = T, log2Rs = T)))
  gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(get(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS1'))),u_Hs_canonical_104$ensembl_gene_id)]
  assign(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS1'), cbind(get(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS1')),gene_name))
}
for(t in c('XS4')){
  list_RNA = sample_info[sample_info$Type == 'RNA'&sample_info$treatment%in%c('XS2', t), c('Name','treatment')]
  list_RIBO = sample_info[sample_info$Type == 'RIBO'&sample_info$treatment%in%c('XS2', t), c('Name','treatment')]
  count_table_RNA = subset(CodonOptimizeXBB_featureCounts_combined_for_xtail, select = as.vector(list_RNA$Name)) %>% floor()
  count_table_RIBO = subset(CodonOptimizeXBB_featureCounts_combined_for_xtail, select = as.vector(list_RIBO$Name)) %>% floor()
  if(all(list_RNA$treatment == list_RIBO$treatment)) condition = as.vector(list_RNA$treatment) else print('WRONG!!!')
  res_CodonOptimizeXBB_RiboSeq = xtail(count_table_RNA, count_table_RIBO, condition, baseLevel = 'XS2')
  assign(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS2'), data.frame(resultsTable(res_CodonOptimizeXBB_RiboSeq, log2FCs = T, log2Rs = T)))
  gene_name = u_Hs_canonical_104$external_gene_name[match(rownames(get(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS2'))),u_Hs_canonical_104$ensembl_gene_id)]
  assign(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS2'), cbind(get(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',t,'vsXS2')),gene_name))
}

save(list = c(ls(pattern = "^RNAFC"), ls(pattern = "^RIBOFC"), ls(pattern = "^TEFC")), file = "Fig4/Results/differental_analysis.Rdata")


