library(data.table)
library(ggpubr)
library(gridExtra)
library(EnvStats)
library(patchwork)
library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(ComplexHeatmap)
library(tidyverse)
library(Biobase)
library(Mfuzz)
library(org.Hs.eg.db)
library(clusterProfiler)

load('Fig4/Results/differental_analysis.Rdata') # Derived from the file differential_analysis_prepare.R
RPKM=fread("Fig4/Results/mRNA_RPF_TE_RPKM.csv")
background_genes <-unique(filter(RPKM,(mRNA_Ctrl_1+mRNA_Ctrl_2)/2>1,grepl("ENS",Geneid))$Geneid)

# RPF
rpf <- RPKM %>% dplyr::select(starts_with("Ribo")  | contains("Geneid"))

head(rpf)

rpf_long <- rpf %>%
  pivot_longer(
    -Geneid,
    names_to = "sample",
    values_to = "rpkm"
  )

head(rpf_long)
rpf_long <- rpf_long %>%
  mutate(group = str_remove(sample, "_\\d+$"))  # 去掉结尾的 _数字

rpf_mean <- rpf_long %>%
  group_by(Geneid, group) %>%
  summarise(mean_rpkm = mean(rpkm, na.rm = TRUE), .groups = "drop")

rpf <- rpf_mean %>%
  pivot_wider(names_from = group, values_from = mean_rpkm) %>%
  column_to_rownames("Geneid")

head(rpf)

eset_rpf <- new("ExpressionSet", exprs = as.matrix(rpf))
eset_rpf <- filter.std(eset_rpf, min.std = 0) 

eset_rpf_std <- standardise(eset_rpf)
m_rpf <- 1.5

set.seed(42)
k <- 4
cl_rpf <- mfuzz(eset_rpf_std, c = k, m = m_rpf)

mfuzz.plot(
  eset_rpf_std, 
  cl = cl_rpf, 
  mfrow = c(2,2),
  time.labels = colnames(eset_rpf_std),
  new.window=F
)

rpf_cluster <- data.frame(
  GeneID = rownames(eset_rpf_std),
  Cluster = cl_rpf$cluster
)
table(rpf_cluster$Cluster)

sig_genes = c()
for(p in c('XS1_vs_Ctrl','XS2_vs_Ctrl','XS4_vs_Ctrl','XS2_vs_XS1','XS4_vs_XS1','XS4_vs_XS2')){
  df_sig = get(paste0('RIBOFC_',p)) %>% filter(pvalue<0.05) %>% rownames()
  print(paste0(p,':',length(df_sig)))
  sig_genes = unique(c(sig_genes, df_sig))
}

sig_genes_rpf = sig_genes

all_clusters <- unique(rpf_cluster$Cluster)

for (cl in all_clusters) {
  gene_list <- intersect(filter(rpf_cluster, Cluster == cl)$GeneID, sig_genes)

  go_enrichment <- clusterProfiler::enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    universe = background_genes,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    keyType = "ENSEMBL",
    ont = "ALL", 
    readable = T
  )
  
  assign(paste0('RPF_strict_GO_cluster',cl),go_enrichment)
  
  cat("luster", cl, "finished | effective gene number：", length(gene_list), "| GO term number：", nrow(go_enrichment@result), "\n")
}


# mRNA
mrna <- RPKM %>% dplyr::select(starts_with("mRNA")  | contains("Geneid"))
head(mrna)

mrna_long <- mrna %>%
  pivot_longer(
    -Geneid,
    names_to = "sample",
    values_to = "rpkm"
  )

head(mrna_long)

mrna_long <- mrna_long %>%
  mutate(group = str_remove(sample, "_\\d+$"))  

mrna_mean <- mrna_long %>%
  group_by(Geneid, group) %>%
  summarise(mean_rpkm = mean(rpkm, na.rm = TRUE), .groups = "drop")

mrna <- mrna_mean %>%
  pivot_wider(names_from = group, values_from = mean_rpkm) %>%
  column_to_rownames("Geneid")


eset_mrna <- new("ExpressionSet", exprs = as.matrix(mrna))
eset_mrna <- filter.std(eset_mrna, min.std = 0)

eset_mrna_std <- standardise(eset_mrna)
m_mrna <- 1.5
set.seed(42)
k <- 4
cl_mrna <- mfuzz(eset_mrna_std, c = k, m = m_mrna)

mfuzz.plot(
  eset_mrna_std, 
  cl = cl_mrna, 
  mfrow = c(2,2),
  time.labels = colnames(eset_mrna_std),
  new.window=F
)

mrna_cluster <- data.frame(
  GeneID = rownames(eset_mrna_std),
  Cluster = cl_mrna$cluster
)
table(mrna_cluster$Cluster)

background_genes <-unique(filter(RPKM,(mRNA_Ctrl_1+mRNA_Ctrl_2)/2>1,grepl("ENS",Geneid))$Geneid)

sig_genes = c()
for(p in c('XS1_vs_Ctrl','XS2_vs_Ctrl','XS4_vs_Ctrl','XS2_vs_XS1','XS4_vs_XS1','XS4_vs_XS2')){
  df_sig = get(paste0('RNAFC_',p)) %>% filter(pvalue<0.05) %>% rownames()
  print(paste0(p,':',length(df_sig)))
  sig_genes = unique(c(sig_genes, df_sig))
}

sig_genes_mrna = sig_genes

all_clusters <- unique(mrna_cluster$Cluster)

for (cl in all_clusters) {
  gene_list <- intersect(filter(mrna_cluster, Cluster == cl)$GeneID, sig_genes)

  go_enrichment <- clusterProfiler::enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    universe = background_genes,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    keyType = "ENSEMBL",
    ont = "ALL", 
    readable = T
  )
  
  assign(paste0('mRNA_strict_GO_cluster',cl),go_enrichment)
  
  cat("luster", cl, "finished | effective gene number：", length(gene_list), "| GO term number：", nrow(go_enrichment@result), "\n")
}

# TE
te <- RPKM %>% dplyr::select(starts_with("TE")  | contains("Geneid"))
te <- te %>%filter(!if_any(where(is.numeric), is.infinite))

head(te)

te_long <- te %>%
  pivot_longer(
    -Geneid,
    names_to = "sample",
    values_to = "rpkm"
  )

head(te_long)

te_long <- te_long %>%
  mutate(group = str_remove(sample, "_\\d+$"))

te_mean <- te_long %>%
  group_by(Geneid, group) %>%
  summarise(mean_rpkm = mean(rpkm, na.rm = TRUE), .groups = "drop")

te <- te_mean %>%
  pivot_wider(names_from = group, values_from = mean_rpkm) %>%
  column_to_rownames("Geneid")

head(te)

zero_var_rows = apply(te, 1, function(x) var(x) == 0)

eset_te <- new("ExpressionSet", exprs = as.matrix(te[which(!zero_var_rows),]))
eset_te <- filter.std(eset_te, min.std = 0)  # 去掉完全不变的基因

zero_var_rows = apply(te, 1, function(x) var(x) == 0)

eset_te_std <- standardise(eset_te)
m_te <- 1.5

set.seed(42)
k <- 4
cl_te <- mfuzz(eset_te_std, c = k, m = m_te)

mfuzz.plot(
  eset_te_std, 
  cl = cl_te, 
  mfrow = c(2,2),  # 2行3列，根据k调整
  time.labels = colnames(eset_te_std),
  new.window=F
)

te_cluster <- data.frame(
  GeneID = rownames(eset_te_std),
  Cluster = cl_te$cluster
)
table(te_cluster$Cluster)

sig_genes = c()
for(p in c('XS1vsCtrl','XS2vsCtrl','XS4vsCtrl','XS2vsXS1','XS4vsXS1','XS4vsXS2')){
  df_sig = get(paste0('TEFC_CodonOptimizeXBB_RiboSeq_',p)) %>% filter(pvalue_final<0.05) %>% rownames()
  print(paste0(p,':',length(df_sig)))
  sig_genes = unique(c(sig_genes, df_sig))
}

sig_genes_te = sig_genes

all_clusters <- unique(te_cluster$Cluster)

for (cl in all_clusters) {
  gene_list <- intersect(filter(te_cluster, Cluster == cl)$GeneID, sig_genes)
  go_enrichment <- clusterProfiler::enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    universe = te_cluster$GeneID,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    keyType = "ENSEMBL",
    ont = "ALL", 
    readable = T
  )
  
  assign(paste0('TE_strict_GO_cluster',cl),go_enrichment)
  
  cat("luster", cl, "finished | effective gene number：", length(gene_list), "| GO term number：", nrow(go_enrichment@result), "\n")
}



# Fig4D
mfuzz_plot_df = eset_rpf_std@assayData$exprs
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_rpf$cluster[match(rownames(mfuzz_plot_df), names(cl_rpf$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
cluster_re_dict = c('Cluster 4','Cluster 2','Cluster 3','Cluster 1')
mfuzz_plot_df = data.frame(mfuzz_plot_df) %>% mutate(gene_id = rownames(mfuzz_plot_df),cluster_re = cluster_re_dict[cluster])

df = cl_rpf$membership
df_rownames = rownames(df)
df = data.table(df) %>% mutate(gene_id = df_rownames)
df$membership = apply(df[,1:4],1,max)

mfuzz_plot_df$membership = df$membership[match(mfuzz_plot_df$gene_id, df$gene_id)]
mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_rpf)
mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, Ribo_Ctrl:Ribo_XS4, names_to = c('type','condition'),names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

x_extend = scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'), expand = expansion(mult = c(0.3, 0.3)))
spacing = unit(1,'cm')
theme_bold_text = theme(axis.text = element_text(face = 'bold'), axis.title = element_text(face = 'bold'))

Fig4D1=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) +
  geom_line(aes(group = gene_id),alpha = 0.05, linewidth = 0.008, color = '#FBE0DC') +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.02, color = '#EA726A', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#F0948B', outlier.shape = NA, linewidth = 0.5) +
  geom_text(aes(x=2.5, y=2,label = paste0('N = ', gene_num)), data = mfuzz_plot_df%>%group_by(cluster_re)%>%summarise(gene_num=n())) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_text(face = 'bold', size = 12)) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#EA726A",
    linewidth = 1) + ylab('Z-score of RPF') +
  theme(panel.spacing.x=spacing) + x_extend + theme_bold_text

mfuzz_plot_df = data.frame(eset_mrna_std@assayData$exprs)
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_rpf$cluster[match(rownames(mfuzz_plot_df), names(cl_rpf$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
mfuzz_plot_df$gene_id = rownames(mfuzz_plot_df)

cluster_re_dict = c('Cluster 4','Cluster 2','Cluster 3','Cluster 1')
mfuzz_plot_df = mfuzz_plot_df %>% mutate(cluster_re = cluster_re_dict[cluster])

mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_rpf)

mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, mRNA_Ctrl:mRNA_XS4, names_to = c('type','condition'), names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

Fig4D2=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.06, color = '#6495ED', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#87CEFA', outlier.shape = NA, linewidth = 0.5) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_blank()) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#6495ED",
    linewidth = 1) + 
  theme(panel.spacing.x=spacing) + ylab('Z-score of mRNA') + x_extend + theme_bold_text

mfuzz_plot_df = data.frame(eset_te_std@assayData$exprs)
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_rpf$cluster[match(rownames(mfuzz_plot_df), names(cl_rpf$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
mfuzz_plot_df$gene_id = rownames(mfuzz_plot_df)

cluster_re_dict = c('Cluster 4','Cluster 2','Cluster 3','Cluster 1')
mfuzz_plot_df = mfuzz_plot_df %>% mutate(cluster_re = cluster_re_dict[cluster])

mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_rpf)

mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, TE_Ctrl:TE_XS4, names_to = c('type','condition'), names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

Fig4D3=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.06, color = '#629C35', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#D0DD97', outlier.shape = NA, linewidth = 0.5) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_blank()) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#629C35",
    linewidth = 1) + 
  theme(panel.spacing.x=spacing) + ylab('Z-score of TE') + x_extend + theme_bold_text

Fig4D = Fig4D1/Fig4D2/Fig4D3

# Fig4EF

RPF_strict_GO_cluster4_df = RPF_strict_GO_cluster4@result
RPF_strict_GO_cluster4_df$Description = gsub('endoplasmic reticulum','ER',RPF_strict_GO_cluster4_df$Description)
RPF_strict_GO_cluster4_df = RPF_strict_GO_cluster4_df %>% separate(GeneRatio, into = c('changed_gene','changed_total_gene'), sep = '/', remove = F) %>% mutate(changed_gene=as.numeric(changed_gene), changed_total_gene=as.numeric(changed_total_gene))
selected_GO = c('GO:0008380','GO:0042254','GO:0042255','GO:0006399','GO:0034728','GO:0034976','GO:0140467','GO:0006400','GO:0070059','GO:0042149','GO:0006983')


Fig4E = ggplot(RPF_strict_GO_cluster4_df%>%filter(ID%in%selected_GO) %>% arrange(-p.adjust) %>% mutate(Description = factor(Description, levels=Description)), aes(x=changed_gene/changed_total_gene, y = Description)) +
  geom_point(aes(colour = -log10(p.adjust), size = changed_gene)) + theme_bw() +
  scale_y_discrete(labels=function(x) str_wrap(x, width=40), position = "right") + 
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.12)), breaks = c(0,0.04,0.08,0.12))+
  scale_color_gradient(low = '#F0aaaa', high = '#E72226', breaks = c(2, 4, 6, 8, 10), labels = parse(text = paste0("10^-", c(2, 4, 6, 8, 10))), limits = c(1, 10.5), oob = scales::squish) +
  scale_size_continuous(range = c(0.2,5), breaks = c(50, 100, 150)) + 
  labs(y=element_blank()) + xlab('Gene ratio') + ggtitle('RPF Cluster 1') + labs(size='Gene count', color='P.adj') +
  theme(axis.text.y = element_text(size=12, face = 'bold'), axis.text.x = element_text(size=12, face = 'bold'), plot.title = element_text(face = 'bold')) +
  guides(color = guide_colorbar(order = 1))

RPF_strict_GO_cluster2_df = RPF_strict_GO_cluster2@result
RPF_strict_GO_cluster2_df = RPF_strict_GO_cluster2_df %>% separate(GeneRatio, into = c('changed_gene','changed_total_gene'), sep = '/', remove = F) %>% mutate(changed_gene=as.numeric(changed_gene), changed_total_gene=as.numeric(changed_total_gene))
selected_GO = c('GO:0001822','GO:0040008','GO:0045785','GO:0016570','GO:0043201','GO:0006006','GO:2001237','GO:0000281')

Fig4F = ggplot(RPF_strict_GO_cluster2_df%>%filter(ID%in%selected_GO) %>% arrange(-p.adjust) %>% mutate(Description = factor(Description, levels=Description)), aes(x=changed_gene/changed_total_gene, y = Description)) +
  geom_point(aes(colour = -log10(p.adjust), size = changed_gene)) + theme_bw() +
  scale_y_discrete(labels=function(x) str_wrap(x, width=40), position = "right") + 
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.12)))+
  scale_color_gradient(low = '#F0aaaa', high = '#E72226', breaks = c(2, 4, 6, 8, 10), labels = parse(text = paste0("10^-", c(2, 4, 6, 8, 10))), limits = c(1, 10.5), oob = scales::squish) +
  scale_size_continuous(range = c(0.2,5), breaks = c(50, 100, 150)) + 
  labs(y=element_blank()) + xlab('Gene ratio') + ggtitle('RPF Cluster 2') + labs(size='Gene count', color='P.adj') +
  theme(axis.text.y = element_text(size=12, face = 'bold'), axis.text.x = element_text(size=12, face = 'bold'), plot.title = element_text(face = 'bold')) +
  guides(color = guide_colorbar(order = 1))

# FigS16A
mfuzz_plot_df = eset_mrna_std@assayData$exprs
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_mrna$cluster[match(rownames(mfuzz_plot_df), names(cl_mrna$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
cluster_re_dict = c('Cluster 1','Cluster 3','Cluster 2','Cluster 4')
mfuzz_plot_df = data.frame(mfuzz_plot_df) %>% mutate(gene_id = rownames(mfuzz_plot_df),cluster_re = cluster_re_dict[cluster])

df = cl_mrna$membership
df_rownames = rownames(df)
df = data.table(df) %>% mutate(gene_id = df_rownames)
df$membership = apply(df[,1:4],1,max)

mfuzz_plot_df$membership = df$membership[match(mfuzz_plot_df$gene_id, df$gene_id)]
mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_mrna)
mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, mRNA_Ctrl:mRNA_XS4, names_to = c('type','condition'),names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

FigS16A1=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) + geom_line(aes(group = gene_id),alpha = 0.05, linewidth = 0.008, color = '#E0FFFF') + ggrastr::geom_jitter_rast(height = 0, alpha = 0.02, color = '#6495ED', size = 0.5, pch=16) + geom_boxplot(fill = NA, color = '#87CEFA', outlier.shape = NA, linewidth = 0.5) + geom_text(aes(x=2.5, y=2,label = paste0('N = ', gene_num)), data = mfuzz_plot_df%>%group_by(cluster_re)%>%summarise(gene_num=n())) + facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') + theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_text(face = 'bold', size = 12), axis.text = element_text(face = 'bold'), axis.title = element_text(face = 'bold')) + labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
  aes(x = as.numeric(condition),y = expression, group = 1), 
  data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
  span = 1,
  se = FALSE,
  color = "#6495ED",
  linewidth = 1,          
) + theme(panel.spacing.x=unit(0.5,'cm')) + ylab('Z-score of mRNA') + scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

mfuzz_plot_df = data.frame(eset_rpf_std@assayData$exprs)
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_mrna$cluster[match(rownames(mfuzz_plot_df), names(cl_mrna$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
mfuzz_plot_df$gene_id = rownames(mfuzz_plot_df)

cluster_re_dict = c('Cluster 1','Cluster 3','Cluster 2','Cluster 4')
mfuzz_plot_df = mfuzz_plot_df %>% mutate(cluster_re = cluster_re_dict[cluster])

mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_rpf)

mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, Ribo_Ctrl:Ribo_XS4, names_to = c('type','condition'), names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

FigS16A2=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.06, color = '#EA726A', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#F0aaaa', outlier.shape = NA, linewidth = 0.5) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_blank(), axis.text = element_text(face = 'bold'), axis.title = element_text(face = 'bold')) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#EA726A",
    linewidth = 1) + 
  theme(panel.spacing.x=unit(0.5,'cm')) + ylab('Z-score of RPF') + scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))


mfuzz_plot_df = data.frame(eset_te_std@assayData$exprs)
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_mrna$cluster[match(rownames(mfuzz_plot_df), names(cl_mrna$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
mfuzz_plot_df$gene_id = rownames(mfuzz_plot_df)

cluster_re_dict = c('Cluster 1','Cluster 3','Cluster 2','Cluster 4')
mfuzz_plot_df = mfuzz_plot_df %>% mutate(cluster_re = cluster_re_dict[cluster])

mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_rpf)

mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, TE_Ctrl:TE_XS4, names_to = c('type','condition'), names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

FigS16A3=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.06, color = '#629C35', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#D0DD97', outlier.shape = NA, linewidth = 0.5) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_blank(), axis.text = element_text(face = 'bold'), axis.title = element_text(face = 'bold')) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#629C35",
    linewidth = 1) + 
  theme(panel.spacing.x=unit(0.5,'cm')) + ylab('Z-score of TE') + scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

FigS16A = FigS16A1/FigS16A2/FigS16A3

# FigS16BC
mRNA_strict_GO_cluster1_df = mRNA_strict_GO_cluster1@result
mRNA_strict_GO_cluster1_df$Description = gsub('endoplasmic reticulum','ER',mRNA_strict_GO_cluster1_df$Description)
mRNA_strict_GO_cluster1_df = mRNA_strict_GO_cluster1_df %>% separate(GeneRatio, into = c('changed_gene','changed_total_gene'), sep = '/', remove = F) %>% mutate(changed_gene=as.numeric(changed_gene), changed_total_gene=as.numeric(changed_total_gene))
selected_GO = c('GO:0042254','GO:0070059','GO:0006400','GO:0034620','GO:0006399','GO:0140467','GO:0034976')

FigS16B = ggplot(mRNA_strict_GO_cluster1_df%>%filter(ID%in%selected_GO) %>% arrange(-p.adjust) %>% mutate(Description = factor(Description, levels=Description)), aes(x=changed_gene/changed_total_gene, y = Description)) + geom_point(aes(colour = -log10(p.adjust), size = changed_gene)) + theme_bw() + scale_y_discrete(labels=function(x) str_wrap(x, width=20), position = 'right') + scale_color_gradient(low = '#87DDFF', high = '#6495ED', breaks = 2:6, labels = parse(text = paste0("10^-", 2:6)), limits = c(1.5, 6.5), oob = scales::squish) + labs(y=element_blank()) + xlab('Gene ratio') + ggtitle('Cluster 1') + labs(size='Gene count', color='P.adj') + theme(axis.text.y = element_text(size=14), axis.text.x = element_text(size=11), plot.title = element_text(face = 'bold')) + guides(color = guide_colorbar(order = 1)) + scale_size_continuous(range = c(0.2,5), breaks = c(20,40,60,80))

mRNA_strict_GO_cluster3_df = mRNA_strict_GO_cluster3@result
mRNA_strict_GO_cluster3_df = mRNA_strict_GO_cluster3_df %>% separate(GeneRatio, into = c('changed_gene','changed_total_gene'), sep = '/', remove = F) %>% mutate(changed_gene=as.numeric(changed_gene), changed_total_gene=as.numeric(changed_total_gene))
selected_GO = c('GO:0045333','GO:0006119','GO:0006754','GO:0010720','GO:0009260','GO:0000281','GO:2000147')

FigS16C = ggplot(mRNA_strict_GO_cluster3_df%>%filter(ID%in%selected_GO) %>% arrange(-p.adjust) %>% mutate(Description = factor(Description, levels=Description)), aes(x=changed_gene/changed_total_gene, y = Description)) + geom_point(aes(colour = -log10(p.adjust), size = changed_gene)) + theme_bw() + scale_y_discrete(labels=function(x) str_wrap(x, width=20), position = 'right') + scale_color_gradient(low = '#87DDFF', high = '#6495ED', , breaks = 2:6, labels = parse(text = paste0("10^-", 2:6)), limits = c(1.5, 6.5), oob = scales::squish) + labs(y=element_blank()) + xlab('Gene ratio') + ggtitle('Cluster 2') + labs(size='Gene count', color='P.adj') + theme(axis.text.y = element_text(size=14), axis.text.x = element_text(size=11), plot.title = element_text(face = 'bold')) + guides(color = guide_colorbar(order = 1)) + scale_size_continuous(range = c(0.2,5), breaks = c(20,40,60,80))

# FigS17A
mfuzz_plot_df = eset_te_std@assayData$exprs
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_te$cluster[match(rownames(mfuzz_plot_df), names(cl_te$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
cluster_re_dict = c('Cluster 4','Cluster 2','Cluster 3','Cluster 1')
mfuzz_plot_df = data.frame(mfuzz_plot_df) %>% mutate(gene_id = rownames(mfuzz_plot_df),cluster_re = cluster_re_dict[cluster])

df = cl_te$membership
df_rownames = rownames(df)
df = data.table(df) %>% mutate(gene_id = df_rownames)
df$membership = apply(df[,1:4],1,max)

mfuzz_plot_df$membership = df$membership[match(mfuzz_plot_df$gene_id, df$gene_id)]
mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_te)
mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, TE_Ctrl:TE_XS4, names_to = c('type','condition'),names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

FigS17A1=ggplot(mfuzz_plot_df_long, aes(x=gsub('X','',condition), y=expression)) +
  geom_line(aes(group = gene_id),alpha = 0.05, linewidth = 0.008, color = '#f0f9c8') +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.02, color = '#629C35', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#D0DD97', outlier.shape = NA, linewidth = 0.5) +
  geom_text(aes(x=2.5, y=2,label = paste0('N = ', gene_num)), data = mfuzz_plot_df%>%group_by(cluster_re)%>%summarise(gene_num=n())) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_text(face = 'bold', size = 12)) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#629C35",
    linewidth = 1) + 
  theme(panel.spacing.x=unit(0.5,'cm')) + ylab('Z-score of TE') + scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

mfuzz_plot_df = data.frame(eset_mrna_std@assayData$exprs)
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_te$cluster[match(rownames(mfuzz_plot_df), names(cl_te$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
mfuzz_plot_df$gene_id = rownames(mfuzz_plot_df)

cluster_re_dict = c('Cluster 4','Cluster 2','Cluster 3','Cluster 1')
mfuzz_plot_df = mfuzz_plot_df %>% mutate(cluster_re = cluster_re_dict[cluster])

mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_rpf)
mfuzz_plot_df = na.omit(mfuzz_plot_df)

mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, mRNA_Ctrl:mRNA_XS4, names_to = c('type','condition'), names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

FigS17A2=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) +
  #geom_line(aes(group = gene_id),alpha = 0.05, linewidth = 0.008, color = '#FBE0DC') +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.06, color = '#6495ED', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#87CEFA', outlier.shape = NA, linewidth = 0.5) +
  #geom_text(aes(x=2.5, y=2,label = paste0('N = ', gene_num)), data = mfuzz_plot_df%>%group_by(cluster_re)%>%summarise(gene_num=n())) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_blank()) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#6495ED",
    linewidth = 1) + 
  theme(panel.spacing.x=unit(0.5,'cm')) + ylab('Z-score of mRNA') + scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

mfuzz_plot_df = data.frame(eset_rpf_std@assayData$exprs)
mfuzz_plot_df = cbind(mfuzz_plot_df, cl_te$cluster[match(rownames(mfuzz_plot_df), names(cl_te$cluster))])
colnames(mfuzz_plot_df)[5] = 'cluster'
mfuzz_plot_df$gene_id = rownames(mfuzz_plot_df)

cluster_re_dict = c('Cluster 4','Cluster 2','Cluster 3','Cluster 1')
mfuzz_plot_df = mfuzz_plot_df %>% mutate(cluster_re = cluster_re_dict[cluster])

mfuzz_plot_df = mfuzz_plot_df %>% filter(gene_id %in% sig_genes_rpf)
mfuzz_plot_df = na.omit(mfuzz_plot_df)

mfuzz_plot_df_long = pivot_longer(mfuzz_plot_df, Ribo_Ctrl:Ribo_XS4, names_to = c('type','condition'), names_sep = '_', values_to = 'expression')
mfuzz_plot_df_long$condition = factor(mfuzz_plot_df_long$condition)

FigS17A3=ggplot(mfuzz_plot_df_long, aes(x=condition, y=expression)) +
  ggrastr::geom_jitter_rast(height = 0, alpha = 0.06, color = '#EA726A', size = 0.5, pch=16) +
  geom_boxplot(fill = NA, color = '#F0aaaa', outlier.shape = NA, linewidth = 0.5) +
  facet_wrap(.~cluster_re, nrow = 1) + theme_pubr() + guides(color='none') +
  theme(strip.background = element_rect(fill = NA,colour = NA), strip.text = element_blank()) +
  labs(x=element_blank()) + coord_cartesian(ylim = c(-2,2)) + geom_smooth(
    aes(x = as.numeric(condition),y = expression, group = 1), 
    data = mfuzz_plot_df_long%>%group_by(condition, cluster_re)%>%summarise(expression=median(expression)),
    span = 1, 
    se = FALSE,
    color = "#EA726A",
    linewidth = 1) + 
  theme(panel.spacing.x=unit(0.5,'cm')) + ylab('Z-score of RPF') + scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

FigS17A = FigS17A1/FigS17A2/FigS17A3

# FigS17B
TE_strict_GO_cluster4_df = TE_strict_GO_cluster4@result
TE_strict_GO_cluster4_df$Description = gsub('endoplasmic reticulum','ER',TE_strict_GO_cluster4_df$Description)
TE_strict_GO_cluster4_df = TE_strict_GO_cluster4_df %>% separate(GeneRatio, into = c('changed_gene','changed_total_gene'), sep = '/', remove = F) %>% mutate(changed_gene=as.numeric(changed_gene), changed_total_gene=as.numeric(changed_total_gene))
selected_GO = c('GO:0006119','GO:0045333','GO:0006754','GO:0042254','GO:0006457','GO:2001242','GO:0042797')

FigS17B = ggplot(TE_strict_GO_cluster4_df%>%filter(ID%in%selected_GO) %>% arrange(-p.adjust) %>% mutate(Description = factor(Description, levels=Description)), aes(x=changed_gene/changed_total_gene, y = Description)) +
  geom_point(aes(colour = -log10(p.adjust), size = changed_gene)) + theme_bw() +
  scale_y_discrete(labels=function(x) str_wrap(x, width=60), position = 'right') +
  scale_color_gradient(low = '#D0DD97', high = '#629C35', breaks = c(2, 4, 6, 8, 10), labels = parse(text = paste0("10^-", c(2, 4, 6, 8, 10)))) +
  labs(y=element_blank()) + xlab('Gene ratio') + ggtitle('Cluster 1') + labs(size='Gene count', color='P.adj') +
  theme(axis.text.y = element_text(size=14), axis.text.x = element_text(size=13), plot.title = element_text(face = 'bold')) +
  guides(color = guide_colorbar(order = 1))




Fig4D
Fig4E
Fig4F

FigS16A
FigS16B
FigS16C

FigS17A
FigS17B