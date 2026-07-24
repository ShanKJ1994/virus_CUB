library(data.table)
library(tidyverse)
library(ggpubr)

RPKM=fread("Fig4/Results/mRNA_RPF_TE_RPKM.csv")
Hs_ribosomal_protein_id=fread('Fig4/Results/ribosomal_protein_Hs_id_with_canonical_transcript.csv')

RPKM_host = gather(RPKM[grepl('(ENSG)|(NeoR)',Geneid)], sample, level, mRNA_Ctrl_1:TE_XS4_2)
RPKM_host$class=gsub("_\\d+","",RPKM_host$sample)
RPKM_host$type=gsub("_\\w+","",RPKM_host$sample)
RPKM_host$type=gsub("mono","Ribo",RPKM_host$type)
RPKM_host = RPKM_host %>% separate(sample, into = c(NA,'condition','rep'), remove = F)

mid=gather(RPKM[Geneid == 'spike'], sample, level, mRNA_Ctrl_1:TE_XS4_2)
mid$class=gsub("_\\d+","",mid$sample)
mid$type=gsub("_\\w+","",mid$sample)
mid = mid %>% separate(sample, into = c(NA,'condition','rep'), remove = F)
mid = mid %>% filter(condition!='Ctrl')
mid_spike = mid
mid_mRNA = filter(mid, type=='mRNA')
mid_TE = filter(mid, type=='TE')
mid_Ribo = filter(mid, type=='Ribo')

col_Human_genes = '#888888'
col_Human_RP_genes = '#A46300'

mid = mid_mRNA
mid$condition = gsub('S','',mid$condition)
mid$level = log2(mid$level)
host = RPKM_host %>% filter(type=='mRNA', condition!='Ctrl') %>% group_by(Geneid) %>% summarise(level = mean(level), .groups = 'drop')
host$level = log2(host$level)
host$x = 'Human\ngenes'
host2 = host %>% filter(Geneid %in% Hs_ribosomal_protein_id$ensembl_gene_id)
host2$x = 'Human\nRP genes'
host = bind_rows(host, host2)
host$condition = 'ref'

signif_df = tibble(condition = c('X1','X2','X4'), Letters = c('c','b','a'))
se <- function(x) sd(x)/sqrt(length(x))
mid2=mid %>% group_by(external_gene_name,class,condition) %>% dplyr::summarise(mean=mean(level),se=se(level))
mid2 = mid2 %>% left_join(signif_df[,c('condition','Letters')], by = 'condition')

Fig4A = ggplot(mid2,aes(x=condition,y=mean)) +
  geom_bar(aes(fill=condition), position="dodge", stat="identity", colour="black",width = 0.5)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.2,col='red',size=.2,position=position_dodge(.5))+
  geom_text(aes(label = Letters), vjust = -1, size = 6) +
  geom_violin(aes(x=x,y=level, fill=x), data = host, alpha = 0.6, width = 0.75, color = NA)+
  geom_boxplot(aes(x=x,y=level,color=x), data = host, width = 0.2, outlier.shape = NA)+
  theme_pubr()+coord_cartesian(ylim=c(0,14))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  labs(y='mRNA RPKM (log2)')+
  facet_grid(cols=vars(condition=='ref'), scales = 'free_x', space = 'free')+
  theme(strip.text.x.top = element_blank(), axis.title.x = element_blank()) + scale_fill_manual(values =  c(
    "Ctrl" = "#575959",
    "X1"  = "#E08251",
    "X2"  = "#9D71A6",
    "X4"  = "#4088B4",
    "Human\ngenes" = col_Human_genes,
    "Human\nRP genes" = col_Human_RP_genes)) +
  scale_color_manual(values =  c(
    "Human\ngenes" = col_Human_genes,
    "Human\nRP genes" = col_Human_RP_genes)) + guides(fill='none', color = 'none') + scale_x_discrete(labels = c('Ctrl'='Ctrl','X1'='S1','X2'='S2','X4'='S4'))

mid = mid_Ribo
mid$condition = gsub('XS','X',mid$condition)
mid$level = log2(mid$level)
host = RPKM_host %>% filter(type=='Ribo', condition!='Ctrl') %>% group_by(Geneid) %>% summarise(level = mean(level), .groups = 'drop')
host$level = log2(host$level)
host$x = 'Human\ngenes'
host2 = host %>% filter(Geneid %in% Hs_ribosomal_protein_id$ensembl_gene_id)
host2$x = 'Human\nRP genes'
host = bind_rows(host, host2)
host$condition = 'ref'

signif_df = tibble(condition = c('X1','X2','X4'), Letters = c('c','b','a'))
se <- function(x) sd(x)/sqrt(length(x))
mid2=mid %>% group_by(external_gene_name,class,condition) %>% dplyr::summarise(mean=mean(level),se=se(level))
mid2 = mid2 %>% left_join(signif_df[,c('condition','Letters')], by = 'condition')

Fig4B = ggplot(mid2,aes(x=condition,y=mean)) +
  geom_bar(aes(fill=condition), position="dodge", stat="identity", colour="black",width = 0.5)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.2,col='red',size=.2,position=position_dodge(.5))+
  geom_text(aes(label = Letters), vjust = -1, size = 6) +
  #  facet_grid(external_gene_name ~ type,scales = "free")+
  geom_violin(aes(x=x,y=level, fill=x), data = host, alpha = 0.6, width = 0.75, color = NA)+
  geom_boxplot(aes(x=x,y=level,color=x), data = host, width = 0.2, outlier.shape = NA)+
  theme_pubr()+coord_cartesian(ylim=c(0,14))+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  labs(y='RPF RPKM (log2)')+
  facet_grid(cols=vars(condition=='ref'), scales = 'free_x', space = 'free')+
  theme(strip.text.x.top = element_blank(), axis.title.x = element_blank()) + scale_fill_manual(values =  c(
    "Ctrl" = "#575959",
    "X1"  = "#E08251",
    "X2"  = "#9D71A6",
    "X4"  = "#4088B4",
    "Human\ngenes" = col_Human_genes,
    "Human\nRP genes" = col_Human_RP_genes)) +
  scale_color_manual(values =  c(
    "Human\ngenes" = col_Human_genes,
    "Human\nRP genes" = col_Human_RP_genes)) + guides(fill='none', color = 'none') + 
  scale_x_discrete(labels = c('Ctrl'='Ctrl','X1'='S1','X2'='S2','X4'='S4'))

mid = mid_TE
mid$condition = gsub('S','',mid$condition)
host = RPKM_host %>% filter(type=='TE', condition!='Ctrl') %>% group_by(Geneid) %>% summarise(level = mean(level), .groups = 'drop')

host$x = 'Human\ngenes'
host2 = host %>% filter(Geneid %in% Hs_ribosomal_protein_id$ensembl_gene_id)
host2$x = 'Human\nRP genes'
host = bind_rows(host, host2)
host$condition = 'ref'

signif_df = tibble(condition = c('X1','X2','X4'), Letters = c('b','a','a'))
se <- function(x) sd(x)/sqrt(length(x))
mid2=mid %>% group_by(external_gene_name,class,condition) %>% dplyr::summarise(mean=mean(level),se=se(level))
mid2 = mid2 %>% left_join(signif_df[,c('condition','Letters')], by = 'condition')

Fig4C=ggplot(mid2,aes(x=condition,y=mean)) +
  geom_bar(aes(fill=condition),position="dodge", stat="identity", colour="black",width = 0.5)+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=.2,col='red',size=.2,position=position_dodge(.5))+
  geom_text(aes(label = Letters), vjust = -1, size=6) +
  #  facet_grid(external_gene_name ~ type,scales = "free")+
  geom_violin(aes(x=x,y=level,fill=x), data = host, alpha = 0.6, width = 0.75, color = NA)+
  geom_boxplot(aes(x=x,y=level,color=x), data = host, width = 0.2, outlier.shape = NA)+
  theme_pubr()+scale_y_continuous(expand = expansion(mult = c(0.05, 2)))+
  labs(y='TE', x='')+
  facet_grid(cols=vars(condition=='ref'), scales = 'free_x', space = 'free')+
  ylim(0,10)+coord_cartesian(ylim = c(0,3))+
  theme(strip.text.x.top = element_blank()) + scale_fill_manual(values =  c(
    "Ctrl" = "#575959",
    "X1"  = "#E08251",
    "X2"  = "#9D71A6",
    "X4"  = "#4088B4",
    "Human\ngenes" = col_Human_genes,
    "Human\nRP genes" = col_Human_RP_genes)) +
  scale_color_manual(values =  c(
    "Human\ngenes" = col_Human_genes,
    "Human\nRP genes" = col_Human_RP_genes)) + guides(fill='none', color = 'none')+ 
  scale_x_discrete(labels = c('Ctrl'='Ctrl','X1'='S1','X2'='S2','X4'='S4'))

Fig4A
Fig4B
Fig4C
