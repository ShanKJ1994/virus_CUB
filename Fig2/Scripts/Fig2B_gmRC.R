library(data.table)
library(tidyverse)
library(ggh4x)
library(ggpubr)
library(patchwork)

gmRC_HV_correlation_APsite = fread('Fig2/Results/Virus_VS_Host_Occupancy3.txt')
gmRC_HV_correlation_APsite = gmRC_HV_correlation_APsite %>% group_by(Sample, Type) %>% mutate(effective_codon_num = n()) %>% ungroup()
study_trans_tmp = c('HCMV'='HCMV\n(2012)', 'zzHCMV'='HCMV\n(2015)', 'IAV'='IAV\n(2016)', 'MHV'='MHV\n(2016)', 'VACV'='VACV\n(2017)', 'zVACV'='VACV\n(2023)' ,'zzzVero_SCV2'='SCV2\n(2020)', 'SCV2_Calu3'='SCV2\n(2021a)', 'SCV2_moi01'='SCV2\n(2021b)', 'EToV'='EToV\n(2018)')
condition_levels_tmp = c('Early','1hpi','2hpi','2.5hpi','3hpi','4hpi','5hpi','8hpi','12hpi','16hpi','24hpi','36hpi','48hpi','72hpi')
gmRC_HV_correlation_AsiteMean_filter = gmRC_HV_correlation_APsite %>% filter(effective_codon_num==61, Type=='Asite_Median') %>% mutate(Sample_tmp = gsub('\\.','_',Sample)) %>% extract(
  col = Sample_tmp, 
  into = c("study", "condition"), 
  regex = "^(.*)_(.*)$",
  remove = FALSE) %>% mutate(study_short = study_trans_tmp[study]) %>% filter(!is.na(study_short)) %>%
  mutate(condition = gsub("^0([0-9])", "\\1", condition)) %>%
  mutate(condition = gsub('h$','hpi', condition)) %>%
  mutate(condition = gsub('25hpi','2.5hpi', condition)) %>%
  mutate(condition = factor(condition, levels = condition_levels_tmp)) %>% filter(Sample!='zzzVero.SCV2_24hpi')
virus_type_trans = c('HCMV'='DNA virus','IAV'='RNA virus','MHV'='RNA virus','VACV'='DNA virus', 'SCV2'='RNA virus', 'EToV'='RNA virus')
get_sig_simple = function(x){
  if(x<0.001) return('***')
  if(x<0.01) return('**')
  if(x<0.05) return('*')
  return('')
}
get_sig_simple = Vectorize(get_sig_simple)
std_cor = stat_cor(method = 'spearman', cor.coef.name = 'rho', output.type = 'text')
study_colors = c('HCMV (2012)'='#cf8b0d', 'HCMV (2015)'='#9ea020','VACV (2017)'='#df66a2','VACV (2023)'='#b274af','EToV (2018)'='#eb746a','IAV (2016)'='#43ae35','MHV (2016)'='#2cb178', 'SCV2 (2020)'='#8989c1','SCV2 (2021a)'='#2ba5df','SCV2 (2021b)'='#1eb5b8')
study_label_order = names(study_colors)
study_short_order_plot = gsub(' ','\n',study_label_order)

gmRC_HV_correlation_AsiteMean_filter_df = gmRC_HV_correlation_AsiteMean_filter %>%
  group_by(study_short, condition) %>% summarise(cor_test = list(cor.test(host, virus, method = 'spearman')), .groups = 'drop') %>%
  mutate(cor_test = purrr::map(cor_test, broom::tidy)) %>% unnest(cor_test) %>%
  separate(study_short, into = c('virus_name',NA), sep = '\n', remove = F) %>%
  mutate(virus_type = virus_type_trans[virus_name], sig = get_sig_simple(p.value), study_label = gsub('\n',' ', study_short)) %>%
  mutate(study_label = factor(study_label, levels = study_label_order), study_short = factor(study_short, levels = study_short_order_plot))



Fig2B1 = ggplot(gmRC_HV_correlation_AsiteMean_filter_df, aes(y=factor(condition, levels = rev(condition_levels_tmp)), x=estimate, fill=study_label)) +
  geom_col(width = 0.8) + geom_text(aes(label=sig), angle = 0, size = 4, hjust = -0.3, vjust = 0.6) +
  ggh4x::facet_nested(study_short~., switch = 'y', scales = "free_y", space = 'free', strip = strip_nested(text_y = element_text(angle = 0))) +
  theme_pubr() + guides(fill='none') +
  theme(strip.background = element_part_rect(side = "r", fill = NA, colour = "black", linewidth = 0.8), strip.placement = "outside", strip.text = element_text(size = 10, face = 'bold'), axis.text = element_text(face = 'bold')) +
  labs(y=NULL, x='Spearman\'s rho between\nhost and viral gmRC') + coord_cartesian(xlim = c(0,1.05)) +
  scale_fill_manual(values = study_colors)

Fig2B2 = ggplot(gmRC_HV_correlation_AsiteMean_filter%>%filter(study_short=='HCMV\n(2015)', condition=='5hpi'), aes(x = log10(host), y = log10(virus))) + 
  geom_smooth(se = F, color = 'red', linewidth=0.75, method = 'lm') + geom_point(color = '#9ea020', pch=16, size = 2) +
  std_cor + theme_pubr() + labs(x='gmRC of hosts (log10)', y='gmRC of viruses (log10)') +
  coord_cartesian(xlim = c(-0.45,0.1), ylim = c(-0.45,0.1))

Fig2B3 = ggplot(gmRC_HV_correlation_AsiteMean_filter%>%filter(study_short=='IAV\n(2016)', condition=='4hpi'), aes(x = log10(host), y = log10(virus))) +
  geom_smooth(se = F, color = 'red', linewidth=0.75, method = 'lm') + geom_point(color = '#43ae35', pch=16, size = 2) +
  std_cor + theme_pubr() + labs(x='gmRC of hosts (log10)', y='gmRC of viruses (log10)') +
  coord_cartesian(xlim = c(-0.3,0), ylim = c(-0.45,0))

Fig2B4 = ggplot(gmRC_HV_correlation_AsiteMean_filter%>%filter(study_short=='SCV2\n(2020)', condition=='5hpi'), aes(x = log10(host), y = log10(virus))) + 
  geom_smooth(se = F, color = 'red', linewidth=0.75, method = 'lm') + geom_point(color = '#8989c1', pch=16, size = 2) + 
  std_cor + theme_pubr() + labs(x='gmRC of hosts (log10)', y='gmRC of viruses (log10)') +
  coord_cartesian(xlim = c(-0.4,0.1), ylim = c(-0.4,0.1))

Fig2B = wrap_elements(Fig2B1)|wrap_elements(Fig2B2/Fig2B3/Fig2B4)

Fig2B
