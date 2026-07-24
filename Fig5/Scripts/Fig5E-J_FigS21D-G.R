library(data.table)
library(tidyverse)
library(ggpubr)

Codon_stats_extended=fread('Fig5/Results/Codon_stats_extended.csv')
condition_int_trans = c('Ctrl'=0,'XS1'=1,'XS2'=2,'XS4'=3)
Codon_stats_extended$condition_int = condition_int_trans[Codon_stats_extended$condition]

get_sig = function(p){
  if(is.na(p)) return('')
  if(p<0.001) return('∗∗∗')
  if(p<0.01) return('∗∗')
  if(p<0.05) return('∗')
  return('')
}
get_sig = Vectorize(get_sig)

summary(lm(A_site_stalling_rate ~ condition_int + Codon, data = Codon_stats_extended)) -> ab
data.frame(ab$coefficients[order(ab$coefficients[,4]),]) -> ab
ab$q <- p.adjust(ab[,4])
codon_eff = ab %>% rownames_to_column('term') %>% filter(grepl('Codon',term)) %>% mutate(Codon = gsub('Codon','',term))
codon_eff$sig = get_sig(codon_eff$q)
codon_eff = codon_eff %>% add_row(term = 'CodonAAA', Estimate = 0, Pr...t.. = 1, q = 1, Codon='AAA', sig = '')
df_deltaASR = Codon_stats_extended %>% filter(condition%in%c('Ctrl','XS4')) %>%
  pivot_wider(id_cols = c(Codon, AA3), names_from = condition, values_from = c(A_site_stalling_rate, tRNA_charging_rate))
df_deltaASR = df_deltaASR %>% left_join(Codon_stats_extended%>%filter(condition=='Ctrl')%>%dplyr::select(Codon, RSCU), by = 'Codon')
df_deltaASR = df_deltaASR %>% left_join(codon_eff%>%dplyr::select(Codon, Estimate, q, sig), by = 'Codon') 
df_deltaASR = df_deltaASR %>% arrange(-Estimate) %>% mutate(name = paste0(Codon,':',AA3), name = factor(name, levels=rev(name)))
df_deltaASR = df_deltaASR %>% arrange(-Estimate) %>% mutate(name = factor(name, levels=name))
ptA1 = ggplot(df_deltaASR, aes(y=Estimate, x = name)) + geom_col(aes(fill = Estimate > 0), width = 0.8) + theme_pubr() +
  theme(axis.text.x = element_text(family = 'Roboto Mono', angle = 90, hjust = 1, vjust = 0.5), axis.ticks.x = element_blank(), axis.line.x = element_blank())+
  scale_fill_manual(values =  c("TRUE" = "#b2182b","FALSE"  = "#2166ac")) +
  labs(y='Codon-specific effect\non A-site stalling rate', x = NULL) + guides(fill='none')
ptA_sig <- ggplot(df_deltaASR, aes(x = name)) + 
  geom_text(aes(y = 0, label = sig), family = "DejaVu Sans Mono", color = "black", size = 4, angle = 90, hjust = 1, vjust = 0.5) +
  theme_void() +
  coord_cartesian(clip = "off")
ptA2 = ggplot(df_deltaASR, aes(y='RSCU', x = name, fill = log2(RSCU))) + geom_tile(color = "white", linewidth = 0.3) + 
  scale_fill_gradient2(
    low = "#2166ac",     
    mid = "white",    
    high = "#b2182b",  
    midpoint = 0,   
    name = "RSCU (log2)" 
  ) + theme_pubr(legend = 'bottom') + labs(x=NULL, y=NULL) + theme(axis.text.x = element_blank(), axis.line.x = element_blank(),, axis.line.y = element_blank(), axis.ticks.x = element_blank())

Fig4F = (ptA1/ptA_sig/ptA2) + plot_layout(heights = c(5,2,1))


lc_line = function(slope_val, intercept_val, xmin, xmax) geom_segment(
  data = tibble(xmin=xmin, xmax=xmax),
  aes(
    x = xmin, 
    xend = xmax, 
    y = slope_val * xmin + intercept_val,
    yend = slope_val * xmax + intercept_val
  ),
  colour = "red",
  linewidth = 1,
  inherit.aes = FALSE
)

lm_coeff_CR_m1A58 = broom::tidy(lm(tRNA_charging_rate ~ condition_int + m1A58_level, data = Codon_stats_extended))
head(lm_coeff_CR_m1A58)
Fig4E = ggplot(Codon_stats_extended, aes(y=tRNA_charging_rate, x = m1A58_level)) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(y='tRNA charging rate', x='m1A58 level') +
  lc_line(lm_coeff_CR_m1A58%>%filter(term=='m1A58_level')%>%pull(estimate),
          lm_coeff_CR_m1A58%>%filter(term=='(Intercept)')%>%pull(estimate), 0, 0.5) +
  scale_color_manual(values =  c(
  "Ctrl" = "#575959",
  "XS1"  = "#E08251",
  "XS2"  = "#9D71A6",
  "XS4"  = "#4088B4"),
  labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) + coord_cartesian(ylim = c(0.4,1))


ASR_lm_coeff = broom::tidy(lm(A_site_stalling_rate~condition_int+tRNA_charging_rate*m1A58_level + RSCU, data = Codon_stats_extended))
ASR_lm_coeff %>% filter(term%in%c('tRNA_charging_rate','m1A58_level','tRNA_charging_rate:m1A58_level'))
Fig4G = ggplot(Codon_stats_extended, aes(x=tRNA_charging_rate, y = A_site_stalling_rate)) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(x='tRNA charging rate', y='A-site stalling rate') +
  lc_line(ASR_lm_coeff%>%filter(term=='tRNA_charging_rate')%>%pull(estimate), 
          ASR_lm_coeff%>%filter(term=='(Intercept)')%>%pull(estimate), 0.4, 0.9) +
  scale_color_manual(values =  c(
    "Ctrl" = "#575959",
    "XS1"  = "#E08251",
    "XS2"  = "#9D71A6",
    "XS4"  = "#4088B4"),
    labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) + coord_cartesian(xlim = c(0.4,0.9), ylim = c(0,0.1))

Fig4H = ggplot(Codon_stats_extended, aes(x=m1A58_level, y = A_site_stalling_rate)) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(x='m1A58 level', y='A-site stalling rate') +
  lc_line(ASR_lm_coeff%>%filter(term=='m1A58_level')%>%pull(estimate),
          ASR_lm_coeff%>%filter(term=='(Intercept)')%>%pull(estimate), 0.1, 0.25) +
  scale_color_manual(values =  c(
    "Ctrl" = "#575959",
    "XS1"  = "#E08251",
    "XS2"  = "#9D71A6",
    "XS4"  = "#4088B4"),
    labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) + coord_cartesian(ylim = c(0,0.1))


m1A58_3value = Codon_stats_extended %>% mutate(m1A58_group = ceiling(3*rank(m1A58_level)/n())) %>% group_by(m1A58_group) %>% summarise(m1A58_level = mean(m1A58_level)) %>% pull(m1A58_level)
charging_3value = Codon_stats_extended %>% mutate(charging_group = ceiling(3*rank(tRNA_charging_rate)/n())) %>% group_by(charging_group) %>% summarise(tRNA_charging_rate = mean(tRNA_charging_rate)) %>% pull(tRNA_charging_rate)
coeff_dict = setNames(ASR_lm_coeff$estimate, ASR_lm_coeff$term)
condition_int_mean = mean(Codon_stats_extended$condition_int)
RSCU_mean = mean(Codon_stats_extended$RSCU)
m = matrix(0, ncol = 3, nrow = 3)
colnames(m) = c('low\nm1A58 level','medium\nm1A58 level','high\nm1A58 level')
rownames(m) = c('low\ncharging rate','medium\ncharging rate','high\ncharging rate')
for(a in 1:3){
  for(b in 1:3) m[a,b] = coeff_dict['(Intercept)']+coeff_dict['tRNA_charging_rate']*charging_3value[a]+coeff_dict['m1A58_level']*m1A58_3value[b] + coeff_dict['tRNA_charging_rate:m1A58_level']*charging_3value[a]*m1A58_3value[b] + coeff_dict['condition_int']*condition_int_mean + coeff_dict['RSCU']*RSCU_mean
}
m = m[c(3,2,1),]
Fig4I = Heatmap(m, cluster_rows = F, cluster_columns = F, row_names_side = 'left', column_names_rot = 0, column_names_centered = T,row_names_centered = T, heatmap_legend_param = list(title = 'Predicted stalling rate',direction = "horizontal",title_position = "topleft"), col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50))

ASR_lm_coeff_AA = broom::tidy(lm(A_site_stalling_rate~condition_int+tRNA_charging_rate*m1A58_level + RSCU + AA, data = Codon_stats_extended))
ASR_lm_coeff_AA %>% filter(term%in%c('RSCU'))
Fig4J = ggplot(Codon_stats_extended, aes(x=RSCU, y = A_site_stalling_rate)) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(x='RSCU', y='A-site stalling rate') + 
  lc_line(ASR_lm_coeff_AA%>%filter(term=='RSCU')%>%pull(estimate), 
          ASR_lm_coeff_AA%>%filter(term=='(Intercept)')%>%pull(estimate), 0, 3) +
  scale_color_manual(values =  c(
    "Ctrl" = "#575959",
    "XS1"  = "#E08251",
    "XS2"  = "#9D71A6",
    "XS4"  = "#4088B4"),
    labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) + coord_cartesian(xlim = c(0,3), ylim = c(0,0.1))

isodecoder_table = fread('Fig5/Results/isodecoder_chargingRate_m1A58.csv')
isoacceptor_table = fread('Fig5/Results/isoacceptor_chargingRate_m1A58.csv')
lc_CR_m1A_isod = broom::tidy(lm(charging_rate ~ m1A58_level + condition_int, isodecoder_table))
head(lc_CR_m1A_isod)
lc_CR_m1A_isoa = broom::tidy(lm(charging_rate ~ m1A58_level + condition_int, isoacceptor_table))
head(lc_CR_m1A_isoa)

FigS21D = ggplot(isodecoder_table, aes(y=charging_rate, x = m1A58_level)) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(y='tRNA charging rate', x='m1A58 level') +
  lc_line(lc_CR_m1A_isod%>%filter(term=='m1A58_level')%>%pull(estimate),
          lc_CR_m1A_isod%>%filter(term=='(Intercept)')%>%pull(estimate), 0, 0.6) +
  scale_color_manual(values =  c(
  "Ctrl" = "#575959",
  "XS1"  = "#E08251",
  "XS2"  = "#9D71A6",
  "XS4"  = "#4088B4"),
  labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) + coord_cartesian(ylim = c(0,1))

FigS21E = ggplot(isoacceptor_table, aes(y=charging_rate, x = m1A58_level)) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(y='tRNA charging rate', x='m1A58 level') +
  lc_line(lc_CR_m1A_isoa%>%filter(term=='m1A58_level')%>%pull(estimate),
          lc_CR_m1A_isoa%>%filter(term=='(Intercept)')%>%pull(estimate), 0, 0.6) +
  scale_color_manual(values =  c(
  "Ctrl" = "#575959",
  "XS1"  = "#E08251",
  "XS2"  = "#9D71A6",
  "XS4"  = "#4088B4"),
  labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) + coord_cartesian(ylim = c(0.3,1))

lc_CR_ASTN = broom::tidy(lm(tRNA_charging_rate~log2(A_site_total_num) + condition_int, Codon_stats_extended))
lc_m1A_ASTN = broom::tidy(lm(m1A58_level~log2(A_site_total_num) + condition_int, Codon_stats_extended))

head(lc_m1A_ASTN)
FigS21F=ggplot(Codon_stats_extended, aes(y=m1A58_level, x = log2(A_site_total_num))) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(y='m1A58 level', x='Codon number in genome (log2)') +
  lc_line(lc_m1A_ASTN%>%filter(term=='log2(A_site_total_num)')%>%pull(estimate),
          lc_m1A_ASTN%>%filter(term=='(Intercept)')%>%pull(estimate), 14.6, 18) +
  scale_color_manual(values =  c(
  "Ctrl" = "#575959",
  "XS1"  = "#E08251",
  "XS2"  = "#9D71A6",
  "XS4"  = "#4088B4"),
  labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

head(lc_CR_ASTN)
FigS21G=ggplot(Codon_stats_extended, aes(y=tRNA_charging_rate, x = log2(A_site_total_num))) + geom_point(aes(colour = condition)) +
  theme_pubr(legend = 'right') + labs(y='tRNA charging rate', x='Codon number in genome (log2)') +
  lc_line(lc_CR_ASTN%>%filter(term=='log2(A_site_total_num)')%>%pull(estimate),
          lc_CR_ASTN%>%filter(term=='(Intercept)')%>%pull(estimate), 14.6, 18) +
  scale_color_manual(values =  c(
    "Ctrl" = "#575959",
    "XS1"  = "#E08251",
    "XS2"  = "#9D71A6",
    "XS4"  = "#4088B4"),
    labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')) + coord_cartesian(ylim = c(0.3,1))

Fig4E
Fig4F
Fig4G
Fig4H
draw(Fig4I,heatmap_legend_side = "top")
Fig4J

FigS21D
FigS21E
FigS21F
FigS21G
