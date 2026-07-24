library(ComplexHeatmap)
library(data.table)
library(tidyverse)
library(gridExtra)

## Fonts "Roboto Mono" and "DejaVu Sans Mono" are required to draw the plots.

load('Fig4/Results/differental_analysis.Rdata') # Derived from the file differential_analysis_prepare.R
RPKM = fread("Fig4/Results/mRNA_RPF_TE_RPKM.csv")

tRNAmod_gene_list = c('TRMT1','TRIT1', 'TRMO','TRMT61A','PUS1','NSUN2', 'TRMT11','QTRT2','ELP2','ADAT3')
condition_trans = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')

get_sig = function(p){
  if(is.na(p)) return('')
  if(p<0.001) return('∗∗∗')
  if(p<0.01) return('∗∗')
  if(p<0.05) return('∗')
  return('')
}
get_sig = Vectorize(get_sig)

gene_name_trans2 = function(x){
  name_trans = c('EIF2S1'='eIF2α',
                 'EIF2AK4'='GCN2',
                 'EIF2AK1'='HRI',
                 'EIF2AK2'='PKR',
                 'EIF2AK3'='PERK',
                 'EEF1A1'='eEF1A1')
  name = ifelse(x%in%names(name_trans), name_trans[x], x)
  return(name)
}
gene_name_trans2 = Vectorize(gene_name_trans2)

get_heatmap_sig = function(mat, sig_df, legend_text, cols=NULL, cluster_rows=T, ex2 = NULL, title_position, extra_legend_par = list()){
  if(is.null(cols)) cols = circlize::colorRamp2(c(min(mat, na.rm = T), 0, max(mat, na.rm = T)), c("#2166ac", "white", "#b2182b"))
  if(all(colnames(mat)%in%names(condition_trans))) colnames(mat) = condition_trans[colnames(mat)]
  Heatmap(mat, cluster_columns = F, cluster_rows = cluster_rows, column_names_rot=0, column_names_centered = T, col = cols, border='grey', rect_gp = gpar(col = "grey", lwd = 1), show_row_names = FALSE,
          right_annotation = rowAnnotation(
            RowNames    = anno_text(rownames(mat), gp = gpar(fontsize = 12, fontfamily = "Roboto Mono"), width = max_text_width(rownames(mat))),
            ExtraColumn = anno_text(sig_df$sig, gp = gpar(fontsize = 12, col = sig_df$trend, fontfamily = "DejaVu Sans Mono")),
            ex2 = ex2,
            spacer = anno_empty(width = unit(3, "mm"), border = FALSE),
            gap = unit(2, "mm")),
          width = ncol(mat) * unit(15, "mm"),
          heatmap_legend_param = c(list(title = gt_render(legend_text)), extra_legend_par))
}

mRNA_RPKM_byCondition = RPKM %>% mutate(id = Geneid) %>%
  dplyr::select(id, mRNA_Ctrl_1:mRNA_XS4_2) %>% pivot_longer(mRNA_Ctrl_1:mRNA_XS4_2, names_to = c('type','condition','rep'), names_sep = '_') %>%
  group_by(id, condition) %>% summarise(value = mean(value), .groups = 'drop') %>%
  pivot_wider(id_cols = id, names_from = 'condition', values_from = 'value')
RPF_RPKM_byCondition = RPKM %>% mutate(id = Geneid) %>%
  dplyr::select(id, Ribo_Ctrl_1:Ribo_XS4_2) %>% pivot_longer(Ribo_Ctrl_1:Ribo_XS4_2, names_to = c('type','condition','rep'), names_sep = '_') %>%
  group_by(id, condition) %>% summarise(value = mean(value), .groups = 'drop') %>%
  pivot_wider(id_cols = id, names_from = 'condition', values_from = 'value')

mRNA_DESeq2_vsCtrl = RNAFC_XS1_vs_Ctrl%>% rownames_to_column(var = 'id') %>% dplyr::select(id, gene_name, log2FoldChange) %>% dplyr::rename(`S1/Ctrl` = log2FoldChange) %>%
  full_join(RNAFC_XS2_vs_Ctrl%>% rownames_to_column(var = 'id') %>% dplyr::select(id, gene_name, log2FoldChange) %>% dplyr::rename(`S2/Ctrl` = log2FoldChange)) %>%
  full_join(RNAFC_XS4_vs_Ctrl %>% rownames_to_column(var = 'id') %>% dplyr::select(id, gene_name, log2FoldChange) %>% dplyr::rename(`S4/Ctrl` = log2FoldChange)) %>%
  mutate(`Ctrl/Ctrl`=0, .before=`S1/Ctrl`) %>% rowwise() %>% mutate(trend = case_when(
    all(diff(c_across(`Ctrl/Ctrl`:`S4/Ctrl`)) > 0) ~ "red",  
    all(diff(c_across(`Ctrl/Ctrl`:`S4/Ctrl`)) < 0) ~ "blue",  
    TRUE ~ "#666666")) %>%
  ungroup() %>% left_join(mRNA_RPKM_byCondition, by = 'id')
RPF_DESeq2_vsCtrl = RIBOFC_XS1_vs_Ctrl%>% rownames_to_column(var = 'id') %>% dplyr::select(id, gene_name, log2FoldChange) %>% dplyr::rename(`S1/Ctrl` = log2FoldChange) %>%
  full_join(RIBOFC_XS2_vs_Ctrl%>% rownames_to_column(var = 'id') %>% dplyr::select(id, gene_name, log2FoldChange) %>% dplyr::rename(`S2/Ctrl` = log2FoldChange)) %>%
  full_join(RIBOFC_XS4_vs_Ctrl %>% rownames_to_column(var = 'id') %>% dplyr::select(id, gene_name, log2FoldChange) %>% dplyr::rename(`S4/Ctrl` = log2FoldChange)) %>%
  mutate(`Ctrl/Ctrl`=0, .before=`S1/Ctrl`) %>% rowwise() %>% mutate(trend = case_when(
    all(diff(c_across(`Ctrl/Ctrl`:`S4/Ctrl`)) > 0) ~ "red",  
    all(diff(c_across(`Ctrl/Ctrl`:`S4/Ctrl`)) < 0) ~ "blue",  
    TRUE ~ "#666666")) %>%
  ungroup() %>% left_join(RPF_RPKM_byCondition, by = 'id')

RNAFC_integrated = bind_rows(list('XS1vsCtrl'=RNAFC_XS1_vs_Ctrl%>%rownames_to_column(var='gene_id'),'XS2vsCtrl'=RNAFC_XS2_vs_Ctrl%>%rownames_to_column(var='gene_id'),'XS4vsCtrl'=RNAFC_XS4_vs_Ctrl%>%rownames_to_column(var='gene_id'),'XS2vsXS1'=RNAFC_XS2_vs_XS1%>%rownames_to_column(var='gene_id'),'XS4vsXS1'=RNAFC_XS4_vs_XS1%>%rownames_to_column(var='gene_id'),'XS4vsXS2'=RNAFC_XS4_vs_XS2%>%rownames_to_column(var='gene_id')), .id='pair')
RIBOFC_integrated = bind_rows(list('XS1vsCtrl'=RIBOFC_XS1_vs_Ctrl%>%rownames_to_column(var='gene_id'),'XS2vsCtrl'=RIBOFC_XS2_vs_Ctrl%>%rownames_to_column(var='gene_id'),'XS4vsCtrl'=RIBOFC_XS4_vs_Ctrl%>%rownames_to_column(var='gene_id'),'XS2vsXS1'=RIBOFC_XS2_vs_XS1%>%rownames_to_column(var='gene_id'),'XS4vsXS1'=RIBOFC_XS4_vs_XS1%>%rownames_to_column(var='gene_id'),'XS4vsXS2'=RIBOFC_XS4_vs_XS2%>%rownames_to_column(var='gene_id')), .id='pair')

df_Wald_mRNA_RPKM = RNAFC_integrated %>% group_by(gene_id) %>% summarise(min_p = min(pvalue, na.rm = T), min_p_pair = paste(pair[pvalue==min_p], collapse = ';'), min_padj = min(padj, na.rm = T), min_padj_pair = paste(pair[padj==min_padj], collapse = ';')) %>% rename(id=gene_id) %>% left_join(mRNA_DESeq2_vsCtrl, by = 'id') %>% mutate(sig=get_sig(min_padj), sip=get_sig(min_p))
df_Wald_RPF_RPKM = RIBOFC_integrated %>% group_by(gene_id) %>% summarise(min_p = min(pvalue, na.rm = T), min_p_pair = paste(pair[pvalue==min_p], collapse = ';'), min_padj = min(padj, na.rm = T), min_padj_pair = paste(pair[padj==min_padj], collapse = ';')) %>% rename(id=gene_id) %>% left_join(RPF_DESeq2_vsCtrl, by = 'id') %>% mutate(sig=get_sig(min_padj), sip=get_sig(min_p))

oneScale = c(-0.2,0,0.8)

df = df_Wald_mRNA_RPKM %>% filter(gene_name%in%tRNAmod_gene_list) %>% mutate(gene_name2 = gene_name_trans2(gene_name)) %>% column_to_rownames(var = 'gene_name')
ordered_rows_mRNA = tRNAmod_gene_list[tRNAmod_gene_list%in%df_Wald_mRNA_RPKM$gene_name]
df = df[ordered_rows_mRNA,]
df = df %>% tibble() %>% column_to_rownames(var = 'gene_name2')
Fig4G_mRNA=get_heatmap_sig(df%>%dplyr::select(`S1/Ctrl`:`S4/Ctrl`)%>%rename_with(~ stringr::str_split_i(., "/", 1)),
                           df%>%dplyr::select(sip, trend)%>%rename(sig=sip),"value",
                           cluster_rows = F, cols = circlize::colorRamp2(oneScale, c("#2166ac", "white", "#b2182b")),
                           extra_legend_par = list(direction = "horizontal",title_position = "topleft",at = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8)))


df = df_Wald_RPF_RPKM %>% filter(gene_name%in%tRNAmod_gene_list) %>% mutate(gene_name2 = gene_name_trans2(gene_name)) %>% column_to_rownames(var = 'gene_name')
df = df[ordered_rows_mRNA,]
df = df %>% tibble() %>% column_to_rownames(var = 'gene_name2')
Fig4G_RPF=get_heatmap_sig(df%>%dplyr::select(`S1/Ctrl`:`S4/Ctrl`)%>%rename_with(~ stringr::str_split_i(., "/", 1)), 
                          df%>%dplyr::select(sip, trend)%>%rename(sig=sip),"value", 
                          cluster_rows = F, cols = circlize::colorRamp2(oneScale, c("#2166ac", "white", "#b2182b")),
                          extra_legend_par = list(direction = "horizontal",title_position = "topleft",at = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8)))

Fig4G = arrangeGrob(grobs=lapply(list(Fig4G_mRNA,Fig4G_RPF), function(x) grid.grabExpr(draw(x,heatmap_legend_side = "bottom"))), nrow=1) 
grid.draw(Fig4G)
