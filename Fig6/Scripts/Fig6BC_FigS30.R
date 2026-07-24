library(tidyverse)
library(ggpubr)


WB_24h = fread("Fig6/Results/Wes_protein_quant_24h_eEF1A1_eIF2α.csv",header = T)
colnames(WB_24h)=c("sample_id","eEF1A1","eIF2α","GAPDH","P-eIF2α","eEF1A1_over_GAPDH","eIF2α_phospho_relative_level","eIF2α_over_GAPDH","level_int")
min_eEF1A1_over_GAPDH <- min(WB_24h$eEF1A1_over_GAPDH)
max_eEF1A1_over_GAPDH <- max(WB_24h$eEF1A1_over_GAPDH)
min_eIF2α_phospho_relative_level <- min(WB_24h$eIF2α_phospho_relative_level)
max_eIF2α_phospho_relative_level <- max(WB_24h$eIF2α_phospho_relative_level)
scale_factor_24h <- (max_eEF1A1_over_GAPDH - min_eEF1A1_over_GAPDH) / (max_eIF2α_phospho_relative_level - min_eIF2α_phospho_relative_level)
offset_24h <- min_eEF1A1_over_GAPDH - min_eIF2α_phospho_relative_level * scale_factor_24h

Fig6B = ggplot(WB_24h, aes(x = sample_id)) +
  geom_line(aes(y = eEF1A1_over_GAPDH, color = "eEF1A1_over_GAPDH", group = 2), size = 1.2, lineend = "round") +
  geom_point(aes(y = eEF1A1_over_GAPDH, color = "eEF1A1_over_GAPDH"), size = 3) +
  geom_line(aes(y = eIF2α_phospho_relative_level * scale_factor_24h + offset_24h, color = "eIF2α_phospho_relative_level", group = 1), 
            size = 1.2, lineend = "round") +
  geom_point(aes(y = eIF2α_phospho_relative_level * scale_factor_24h + offset_24h, color = "eIF2α_phospho_relative_level"), 
             size = 3) +
  scale_y_continuous(
    name = "eEF1A1_over_GAPDH",
    sec.axis = sec_axis(
      ~ (. - offset_24h) / scale_factor_24h,
      name = "eIF2α_phospho_relative_level"
    ) 
  )  +
  scale_color_manual(
    values = c("eEF1A1_over_GAPDH" = "#3366CC", "eIF2α_phospho_relative_level" = "#CC6677"),
    labels = c("eEF1A1_over_GAPDH", "eIF2α_phospho_relative_level")
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y.left = element_text(color = "#3366CC", size = 12),
    axis.title.y.right = element_text(color = "#CC6677", size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )+
  guides(color='none')

western_24h_protein_level_wide = WB_24h
condition_trans = c('mock'='Ctrl','S1'='S1','S2'='S2','S3'='S3','S4'='S4','S5'='S5')
western_24h_protein_level_wide$condition = condition_trans[western_24h_protein_level_wide$sample_id]
western_24h_protein_level_wide$text = paste0(western_24h_protein_level_wide$condition,':24h')
western_24h_protein_level_wide = western_24h_protein_level_wide %>%
  mutate(eEF1A1_over_GAPDH_over_Ctrl = eEF1A1_over_GAPDH/eEF1A1_over_GAPDH[condition=='Ctrl'], 
         eIF2α_P_overMax = eIF2α_phospho_relative_level/max(eIF2α_phospho_relative_level),
         eIF2α_P_overCtrl = eIF2α_phospho_relative_level/eIF2α_phospho_relative_level[condition=='Ctrl'])


Fig6C = ggplot(western_24h_protein_level_wide, aes(x = eEF1A1_over_GAPDH_over_Ctrl, y = eIF2α_P_overCtrl)) +
  geom_smooth(method = 'lm', se= F, color = 'lightblue') + theme_pubr() + geom_point(aes(color = condition), size = 2.5) +
  ggrepel::geom_text_repel(aes(label = condition))  + guides(color = 'none') +
  xlab('eEF1A1 protein level\nrelative to Ctrl') + ylab('eIF2α phosphorylation\nlevel relative to Ctrl') +
  stat_cor(label.y.npc = 'bottom')

WB_72h = fread("Fig6/Results/Wes_protein_quant_72h_eEF1A1_eIF2α.csv",header = T)
colnames(WB_72h)=c("sample_id","eEF1A1","eIF2α","GAPDH","P-eIF2α","eEF1A1_over_GAPDH","eIF2α_phospho_relative_level","eIF2α_over_GAPDH","level_int")
min_eEF1A1_over_GAPDH <- min(WB_72h$eEF1A1_over_GAPDH)
max_eEF1A1_over_GAPDH <- max(WB_72h$eEF1A1_over_GAPDH)
min_eIF2α_phospho_relative_level <- min(WB_72h$eIF2α_phospho_relative_level)
max_eIF2α_phospho_relative_level <- max(WB_72h$eIF2α_phospho_relative_level)
scale_factor_72h <- (max_eEF1A1_over_GAPDH - min_eEF1A1_over_GAPDH) / (max_eIF2α_phospho_relative_level - min_eIF2α_phospho_relative_level)
offset_72h <- min_eEF1A1_over_GAPDH - min_eIF2α_phospho_relative_level * scale_factor_72h

FigS30A = ggplot(WB_72h, aes(x = sample_id)) +
  geom_line(aes(y = eEF1A1_over_GAPDH, color = "eEF1A1_over_GAPDH", group = 2), size = 1.2, lineend = "round") +
  geom_point(aes(y = eEF1A1_over_GAPDH, color = "eEF1A1_over_GAPDH"), size = 3) +
  geom_line(aes(y = eIF2α_phospho_relative_level * scale_factor_72h + offset_72h, color = "eIF2α_phospho_relative_level", group = 1), 
            size = 1.2, lineend = "round") +
  geom_point(aes(y = eIF2α_phospho_relative_level * scale_factor_72h + offset_72h, color = "eIF2α_phospho_relative_level"), 
             size = 3) +
  scale_y_continuous(
    name = "eEF1A1_over_GAPDH",
    sec.axis = sec_axis(
      ~ (. - offset_72h) / scale_factor_72h,
      name = "eIF2α_phospho_relative_level"
    ) 
  )  +
  scale_color_manual(
    values = c("eEF1A1_over_GAPDH" = "#3366CC", "eIF2α_phospho_relative_level" = "#CC6677"),
    labels = c("eEF1A1_over_GAPDH", "eIF2α_phospho_relative_level")
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.y.left = element_text(color = "#3366CC", size = 12),
    axis.title.y.right = element_text(color = "#CC6677", size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )+
  guides(color='none')

western_72h_protein_level_wide = WB_72h
condition_trans = c('mock'='Ctrl','S1'='S1','S2'='S2','S3'='S3','S4'='S4','S5'='S5')
western_72h_protein_level_wide$condition = condition_trans[western_72h_protein_level_wide$sample_id]
western_72h_protein_level_wide$text = paste0(western_72h_protein_level_wide$condition,':72h')
western_72h_protein_level_wide = western_72h_protein_level_wide %>%
  mutate(eEF1A1_over_GAPDH_over_Ctrl = eEF1A1_over_GAPDH/eEF1A1_over_GAPDH[condition=='Ctrl'], 
         eIF2α_P_overMax = eIF2α_phospho_relative_level/max(eIF2α_phospho_relative_level),
         eIF2α_P_overCtrl = eIF2α_phospho_relative_level/eIF2α_phospho_relative_level[condition=='Ctrl'])


FigS30B = ggplot(western_72h_protein_level_wide, aes(x = eEF1A1_over_GAPDH_over_Ctrl, y = eIF2α_P_overCtrl)) +
  geom_smooth(method = 'lm', se= F, color = 'lightblue') + theme_pubr() + geom_point(aes(color = condition), size = 2.5) +
  ggrepel::geom_text_repel(aes(label = condition))  + guides(color = 'none') +
  xlab('eEF1A1 protein level\nrelative to Ctrl') + ylab('eIF2α phosphorylation\nlevel relative to Ctrl') +
  stat_cor(label.y.npc = 'bottom')


Fig6B
Fig6C
FigS30A
FigS30B
