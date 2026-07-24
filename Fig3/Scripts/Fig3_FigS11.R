library(data.table)
library(tidyverse)
library(ggpubr)

CAI_mid = fread('Fig3/Results/Human_gene_and_spike_CAI.csv')
head(CAI_mid)
filter(CAI_mid,!grepl("ENS",Tr_ID))
Fig3A = ggplot(filter(CAI_mid,Class=="Host",CodonSum>=150),aes(x=CAI,y=..count../100))+
  geom_histogram(binwidth = 0.005, fill='grey')+
  geom_vline(xintercept = filter(CAI_mid,grepl("XBB",Tr_ID))$CAI)+
  theme_classic()+
  scale_y_continuous(breaks = seq(0,15,by=5),limits = c(0,15))+
  labs(y="Human gene number (x100)")

western_24h_spike_protein_level_wide = fread('Fig3/Results/Wes_protein_quant_spike_24h.csv')
Fig3B = ggplot(western_24h_spike_protein_level_wide, aes(x = sample_id, y = spike_over_GAPDH, fill = sample_id)) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = 'black') +
  stat_summary(geom = 'errorbar', fun.data = mean_se, color = 'black', width = 0.2) +
  theme_pubr() + coord_cartesian(ylim = c(0,6.5)) + ylab('Spike protein level') + guides(fill='none') +
  scale_x_discrete(labels = c('mock'='Ctrl')) + theme(axis.title.x = element_blank()) +
  scale_fill_manual(values =  c(
    "mock" = "#eb746a",
    "S1"  = "#b19b19",
    "S2"  = "#2cad3f",
    "S3"  = "#1eb5b8",
    "S4"  = "#6792cd",
    "S5"  = "#c170aa"))


cell_viability = fread('Fig3/Results/Cell_viability.csv')
puro_S1to5_raw = fread('Fig3/Results/Wes_protein_quant_puromycin.csv')
se = function(x) sd(x)/sqrt(length(x))
cell_viability$condition = gsub('X','S',cell_viability$condition)
cell_viability_mean = cell_viability %>% group_by(time, condition) %>%
  summarise(viability = mean(100*E590_over_Ctrl), se_viability = se(100*E590_over_Ctrl), .groups = 'drop')
puro_S1to5_wide = puro_S1to5_raw %>% pivot_wider(id_cols = c('Sample', 'condition','rep'), names_from = 'Primary', values_from = 'Total Area')
puro_S1to5_wide$rate = puro_S1to5_wide$puromycin/puro_S1to5_wide$GAPDH
puro_S1to5_wide = puro_S1to5_wide %>% group_by(rep) %>% mutate(base = max(rate), relative_rate = rate/base)
puro_S1to5_wide = puro_S1to5_wide %>% group_by(rep) %>% mutate(puro_base = max(puromycin), relative_puro = puromycin/puro_base)
puro_S1to5_mean = puro_S1to5_wide %>% group_by(condition) %>% summarise(puro = mean((puromycin/1e6)), se_puro = se((puromycin/1e6)))
puro_S1to5_mean$condition = gsub('mock','Ctrl',puro_S1to5_mean$condition)
puro_S1to5_mean$condition = gsub('XS','S',puro_S1to5_mean$condition)
puro_S1to5_mean = left_join(puro_S1to5_mean, cell_viability_mean%>%filter(time==24), by = 'condition')


Fig3C = ggplot(puro_S1to5_wide, aes(x = condition, y = puromycin/1e6, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, width = 0.7, color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width = 0.1, color = 'black') +
  theme_pubr() + guides(fill='none') +
  coord_cartesian(ylim = c(5,11)) + labs(x = NULL, y = 'Protein synthesis rate (a.u.)') +
  scale_x_discrete(labels = c('mock'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) +
  scale_fill_manual(values =  c(
    "mock" = "#eb746a",
    "XS1"  = "#b19b19",
    "XS2"  = "#2cad3f",
    "XS3"  = "#1eb5b8",
    "XS4"  = "#6792cd",
    "XS5"  = "#c170aa"))


Fig3D = ggplot(cell_viability%>%filter(time==24), aes(x = condition, y = 100*E590_over_Ctrl, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, width = 0.7, color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width = 0.1, color = 'black') +
  theme_pubr() + guides(fill='none') + ggtitle('24 hours') +
  coord_cartesian(ylim = c(60,105)) + labs(x = NULL, y = 'Cell viability relative to Ctrl (%)') +
  scale_fill_manual(values =  c(
    "Ctrl" = "#eb746a",
    "S1"  = "#b19b19",
    "S2"  = "#2cad3f",
    "S3"  = "#1eb5b8",
    "S4"  = "#6792cd",
    "S5"  = "#c170aa"))

FigS11A = ggplot(cell_viability%>%filter(time==48), aes(x = condition, y = 100*E590_over_Ctrl, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, width = 0.7, color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width = 0.1, color = 'black') +
  theme_pubr() + guides(fill='none') + ggtitle('48 hours') +
  coord_cartesian(ylim = c(30,105)) + labs(x = NULL, y = 'Cell viability relative to Ctrl (%)') +
  scale_fill_manual(values =  c(
    "Ctrl" = "#eb746a",
    "S1"  = "#b19b19",
    "S2"  = "#2cad3f",
    "S3"  = "#1eb5b8",
    "S4"  = "#6792cd",
    "S5"  = "#c170aa"))

FigS11B = ggplot(cell_viability%>%filter(time==72), aes(x = condition, y = 100*E590_over_Ctrl, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, width = 0.7, color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width = 0.1, color = 'black') +
  theme_pubr() + guides(fill='none') + ggtitle('72 hours') +
  coord_cartesian(ylim = c(40,105)) + labs(x = NULL, y = 'Cell viability relative to Ctrl (%)') +
  scale_fill_manual(values =  c(
    "Ctrl" = "#eb746a",
    "S1"  = "#b19b19",
    "S2"  = "#2cad3f",
    "S3"  = "#1eb5b8",
    "S4"  = "#6792cd",
    "S5"  = "#c170aa"))

FigS11C = ggplot(cell_viability%>%filter(time==96), aes(x = condition, y = 100*E590_over_Ctrl, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, width = 0.7, color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width = 0.1, color = 'black') +
  theme_pubr() + guides(fill='none') + ggtitle('96 hours') +
  coord_cartesian(ylim = c(60,105)) + labs(x = NULL, y = 'Cell viability relative to Ctrl (%)') +
  scale_fill_manual(values =  c(
    "Ctrl" = "#eb746a",
    "S1"  = "#b19b19",
    "S2"  = "#2cad3f",
    "S3"  = "#1eb5b8",
    "S4"  = "#6792cd",
    "S5"  = "#c170aa"))

FigS11D = ggplot(cell_viability%>%filter(time==120), aes(x = condition, y = 100*E590_over_Ctrl, fill = condition)) +
  geom_bar(stat = "summary", fun = mean, width = 0.7, color = 'black') +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", width = 0.1, color = 'black') +
  theme_pubr() + guides(fill='none') + ggtitle('120 hours') +
  coord_cartesian(ylim = c(60,105)) + labs(x = NULL, y = 'Cell viability relative to Ctrl (%)') +
  scale_fill_manual(values =  c(
    "Ctrl" = "#eb746a",
    "S1"  = "#b19b19",
    "S2"  = "#2cad3f",
    "S3"  = "#1eb5b8",
    "S4"  = "#6792cd",
    "S5"  = "#c170aa"))


Fig3E = ggplot(puro_S1to5_mean%>%filter(condition!='Ctrl'), aes(x = puro, y = viability)) + geom_smooth(method = 'lm', se= F, color = 'lightblue') +
  geom_errorbar(aes(ymin = viability - se_viability, ymax = viability + se_viability), color = 'grey75', width = 0.03) +
  geom_errorbar(aes(xmin = puro - se_puro, xmax = puro + se_puro), color = 'grey75') + theme_pubr() +
  geom_point(aes(color = condition), size = 2.5) + ggrepel::geom_text_repel(aes(label = condition))  + guides(color = 'none') +
  xlab('Protein synthesis rate (a.u.)') + ylab('Cell viability relative to Ctrl (%)') + coord_cartesian(xlim=c(6.5,8), ylim = c(75,92)) +
  stat_cor(label.y = 78)

Fig3A
Fig3B
Fig3C
Fig3D
Fig3E

FigS11A
FigS11B
FigS11C
FigS11D


