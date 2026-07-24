library(ComplexHeatmap)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggpubr)

## Fonts "Roboto Mono" and "DejaVu Sans Mono" are required to draw the plots.

anticodon_codon_matrix = fread('Fig5/Results/wobble.csv')
anticodon_codon_matrix$efficiency = 1-anticodon_codon_matrix$s
condition_trans = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4')

get_sig = function(p){
  if(is.na(p)) return('')
  if(p<0.001) return('∗∗∗')
  if(p<0.01) return('∗∗')
  if(p<0.05) return('∗')
  return('')
}
get_sig = Vectorize(get_sig)

get_heatmap_sig_rate = function(mat, sig_df, legend_text, cols=NULL, cell_width=15){
  if(is.null(cols)) cols = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("#2166ac", "white", "#b2182b"))
  Heatmap(mat, cluster_columns = F, column_names_rot=0, column_names_centered = T, col = cols, border='grey', rect_gp = gpar(col = "grey", lwd = 1), heatmap_legend_param = list(title = gt_render(legend_text)), show_row_names = FALSE,
          right_annotation = rowAnnotation(
            RowNames    = anno_text(rownames(mat), gp = gpar(fontsize = 12, fontfamily = "Roboto Mono"), width = max_text_width(rownames(mat))),
            ExtraColumn = anno_text(sig_df$sig, gp = gpar(fontsize = 12, col = sig_df$trend, fontfamily = "DejaVu Sans Mono")),
            spacer = anno_empty(width = unit(3, "mm"), border = FALSE),
            gap = unit(2, "mm")),
          width = ncol(mat) * unit(cell_width, "mm"))
}

get_heatmap_sig_rate_transpose = function(mat, sig_df, legend_text, cols=NULL, cell_width=15){
  if(is.null(cols)) cols = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("#2166ac", "white", "#b2182b"))
  Heatmap(mat, cluster_rows = F, col = cols, border='grey', rect_gp = gpar(col = "grey", lwd = 1), heatmap_legend_param = list(title = gt_render(legend_text)), show_column_names = FALSE,row_names_side ='left',
          bottom_annotation = columnAnnotation(
            ColumnNames    = anno_text(format(colnames(mat), justify = "left"), gp = gpar(fontsize = 12, fontfamily = "Roboto Mono"), width = max_text_width(colnames(mat))),
            ExtraColumn = anno_text(sig_df$sig, gp = gpar(fontsize = 12, col = sig_df$trend, fontfamily = "DejaVu Sans Mono")),
            spacer = anno_empty(width = unit(3, "mm"), border = FALSE),
            gap = unit(6, "mm")),
          height = nrow(mat) * unit(cell_width, "mm"))
}

get_charging_rate_plot_2 = function(rate_df, total_df, charged_df, maxit=500, legend_text='Normalized<br/>charging<br/>rate (log<sub>2</sub>)', cell_width = 8){
  tRNA_chargeRate_bySample = rate_df
  tRNA_abundance_filter = total_df
  tRNA_chargedAbundance_filter = charged_df
  
  tRNA_chargeRate_bySample_long = tRNA_chargeRate_bySample %>% pivot_longer(Ctrl_1:S4_2, names_to = 'sample') %>% separate(sample, into = c('condition','rep'), sep = '_')
  
  tRNA_chargeRate_byCondition = tRNA_chargeRate_bySample_long %>% group_by(id, condition) %>% summarise(value = mean(value), .groups = 'drop') %>% pivot_wider(id_cols = id, names_from = 'condition', values_from = 'value')
  tRNA_chargeRate_byCondition_3samples = tRNA_chargeRate_byCondition %>% dplyr::select(-Ctrl)
  
  tRNA_chargeRate_byCondition_vsCtrl = tRNA_chargeRate_byCondition %>% mutate(`S1/Ctrl` = log2(S1/Ctrl),`S2/Ctrl` = log2(S2/Ctrl),`S4/Ctrl` = log2(S4/Ctrl)) %>% dplyr::select(id, `S1/Ctrl`, `S2/Ctrl`, `S4/Ctrl`)
  condition_int_trans = c('Ctrl'=0,'XS1'=1,'XS2'=2,'XS4'=3)
  charging_full_table = tRNA_abundance_filter %>% mutate(type = 'total_abundance') %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample') %>% dplyr::select(id,type,sample,value) %>% separate(sample, into = c('condition','rep'), remove = F)
  charging_full_table = tRNA_chargedAbundance_filter %>% mutate(type = 'charged_abundance') %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample') %>% dplyr::select(id,type,sample,value) %>% separate(sample, into = c('condition','rep'), remove = F) %>% bind_rows(charging_full_table)
  charging_full_table = charging_full_table %>% pivot_wider(names_from = 'type', values_from = 'value') %>% mutate(condition_int = condition_int_trans[condition])
  
  safe_glm_nb <- safely(
    function(df) {
      MASS::glm.nb(floor(charged_abundance) ~ condition_int + offset(log(total_abundance)), data = df, control = glm.control(maxit = maxit))
    }, 
    otherwise = NULL
  )
  charging_change_test = charging_full_table %>% nest(data = -id) %>% 
    mutate(
      model_safe = map(data, safe_glm_nb),
      model = map(model_safe, "result"),
      tidied = map(model, function(m) {
        if (is.null(m)) {
          return(tibble(term = "condition_int", statistic = NA_real_, p.value = NA_real_))
        } else {
          return(broom::tidy(m))
        }
      })
    ) %>% unnest(tidied) %>%
    filter(term == "condition_int") %>%
    dplyr::select(id, estimate, statistic, p.value) %>% 
    mutate(padj = p.adjust(p.value, method = "BH"))
  
  sig_df_4samples = tRNA_chargeRate_byCondition %>% left_join(charging_change_test, by = 'id') %>% mutate(sig = get_sig(padj)) %>% rowwise() %>% mutate(trend = case_when(
    estimate > 0 ~ "#b2182b", 
    estimate < 0 ~ "#2166ac", 
    TRUE ~ "#666666"   
  )) %>% ungroup()
  
  sig_df_vsCtrl = tRNA_chargeRate_byCondition_vsCtrl %>% left_join(sig_df_4samples, by = 'id')
  p2=get_heatmap_sig_rate(tRNA_chargeRate_byCondition_vsCtrl%>%column_to_rownames(var = 'id')%>%rename_with(~ stringr::str_split_i(., "/", 1)), sig_df_vsCtrl, legend_text, cell_width = cell_width)
  p3=get_heatmap_sig_rate_transpose(tRNA_chargeRate_byCondition_vsCtrl%>%column_to_rownames(var = 'id')%>%rename_with(~ stringr::str_split_i(., "/", 1))%>%t(), sig_df_vsCtrl, legend_text, cell_width = cell_width)
  
  
  return(list(sig_df_vsCtrl, p2, p3))
}

get_rate_plot_2_byCondition = function(rate_df, total_df, charged_df, maxit=500, legend_text='stalling<br/>rate<br/>(log2)', cell_width = 8, log2transform = T){
  tRNA_chargeRate_bySample = rate_df
  if(log2transform) tRNA_chargeRate_bySample = tRNA_chargeRate_bySample %>% mutate(across(where(is.numeric), log2))
  tRNA_abundance_filter = total_df
  tRNA_chargedAbundance_filter = charged_df
  
  tRNA_chargeRate_bySample_long = tRNA_chargeRate_bySample %>% pivot_longer(Ctrl:S4, names_to = 'condition')
  
  tRNA_chargeRate_byCondition = tRNA_chargeRate_bySample_long %>% group_by(id, condition) %>% summarise(value = mean(value), .groups = 'drop') %>% pivot_wider(id_cols = id, names_from = 'condition', values_from = 'value')
  
  
  condition_int_trans = c('Ctrl'=0,'XS1'=1,'XS2'=2,'XS4'=3)
  charging_full_table = tRNA_abundance_filter %>% mutate(type = 'total_abundance') %>% pivot_longer(Ctrl:XS4, names_to = 'condition') %>% dplyr::select(id,type,condition,value)
  charging_full_table = tRNA_chargedAbundance_filter %>% mutate(type = 'charged_abundance') %>% pivot_longer(Ctrl:XS4, names_to = 'condition') %>% dplyr::select(id,type,condition,value) %>% bind_rows(charging_full_table)
  charging_full_table = charging_full_table %>% pivot_wider(names_from = 'type', values_from = 'value') %>% mutate(condition_int = condition_int_trans[condition])
  
  safe_glm_nb <- safely(
    function(df) {
      MASS::glm.nb(floor(charged_abundance) ~ condition_int + offset(log(total_abundance)), data = df, control = glm.control(maxit = maxit))
    }, 
    otherwise = NULL
  )
  charging_change_test = charging_full_table %>% nest(data = -id) %>% 
    mutate(
      model_safe = map(data, safe_glm_nb),
      model = map(model_safe, "result"),
      tidied = map(model, function(m) {
        if (is.null(m)) {
          return(tibble(term = "condition_int", statistic = NA_real_, p.value = NA_real_))
        } else {
          return(broom::tidy(m))
        }
      })
    ) %>% unnest(tidied) %>%
    filter(term == "condition_int") %>%
    dplyr::select(id, estimate, statistic, p.value) %>% 
    mutate(padj = p.adjust(p.value, method = "BH"))
  
  sig_df_4samples = tRNA_chargeRate_byCondition %>% left_join(charging_change_test, by = 'id') %>% mutate(sig = get_sig(padj)) %>% rowwise() %>% mutate(trend = case_when(
    estimate > 0 ~ "#b2182b",  
    estimate < 0 ~ "#2166ac", 
    TRUE ~ "#666666" 
  )) %>% ungroup()
  
  sig_df_vsCtrl = tRNA_chargeRate_byCondition %>% left_join(sig_df_4samples, by = 'id')
  mat = tRNA_chargeRate_byCondition %>% column_to_rownames(var = 'id')
  cols = circlize::colorRamp2(c(min(mat), 0.5*(min(mat)+max(mat)), max(mat)), c("#2166ac", "white", "#b2182b"))
  p2=get_heatmap_sig_rate(mat, sig_df_vsCtrl, legend_text, cell_width = cell_width, cols = cols)
  p3=get_heatmap_sig_rate_transpose(mat%>%t(), sig_df_vsCtrl, legend_text, cell_width = cell_width, cols = cols)
  
  return(list(sig_df_vsCtrl, p2, p3))
}

# Fig5A
m1A58_level_misinc_byCodon_bySample = fread('Fig5/Results/m1A58_level_byCodon_bySample.csv')

m1A58_mat_rate = m1A58_level_misinc_byCodon_bySample %>% mutate(id=paste0(Codon,':',AA),sample = paste0(gsub('XS','S',condition),'_',rep)) %>% pivot_wider(id_cols = id, names_from = sample, values_from = m1A58_level)
m1A58_mat_total = m1A58_level_misinc_byCodon_bySample %>% mutate(id=paste0(Codon,':',AA),sample = paste0(condition,'_',rep)) %>% pivot_wider(id_cols = id, names_from = sample, values_from = tRNA_expression_level)
m1A58_mat_m1A58Abs = m1A58_level_misinc_byCodon_bySample %>% mutate(id=paste0(Codon,':',AA),sample = paste0(condition,'_',rep)) %>% pivot_wider(id_cols = id, names_from = sample, values_from = m1A58_abs)
plot_m1A58_level_byCodon_2 = get_charging_rate_plot_2(m1A58_mat_rate, m1A58_mat_total, m1A58_mat_m1A58Abs, legend_text='Normalized m1A58<br/>level (log<sub>2</sub>)', cell_width=8)

Fig5A = plot_m1A58_level_byCodon_2[[3]]

# Fig5BCD, FigS22AB
tRNA_abundance_raw = fread('Fig5/Results/tRNA_isoacceptor_abundance.txt')
colnames(tRNA_abundance_raw) = gsub('/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj/clean_','',colnames(tRNA_abundance_raw))
colnames(tRNA_abundance_raw) = gsub('_R1.fq.gz.unpaired_uniq.bam','',colnames(tRNA_abundance_raw))
colnames(tRNA_abundance_raw) = gsub('X','XS',colnames(tRNA_abundance_raw))
tRNA_abundance_filter = tRNA_abundance_raw %>% filter(!grepl('Escherichia_coli',Anticodon), !grepl('mito',Anticodon)) %>% mutate(isoacceptor = gsub('Homo_sapiens_tRNA-','',Anticodon))
tRNA_abundance_filter = tRNA_abundance_filter %>% filter(rowMeans(across(Ctrl_1:XS4_2))>50)
tRNA_abundance_filter$id = tRNA_abundance_filter$isoacceptor

tRNA_chargeRate_bySample = fread('Fig5/Results/tRNA_ChargeRate_bySample.tsv') %>% rename(id = Gene)
colnames(tRNA_chargeRate_bySample) = gsub('X','S',colnames(tRNA_chargeRate_bySample))
df_tRNA_chargeRate_bySample = tRNA_chargeRate_bySample %>% column_to_rownames(var = 'id')
colnames(df_tRNA_chargeRate_bySample) = gsub('S','XS',colnames(df_tRNA_chargeRate_bySample))
df_tRNA_abundance = tRNA_abundance_filter %>% column_to_rownames(var = 'isoacceptor') %>% dplyr::select(Ctrl_1:XS4_2)
df_tRNA_chargedAbundance = df_tRNA_abundance * (df_tRNA_chargeRate_bySample[rownames(df_tRNA_abundance), colnames(df_tRNA_abundance)])
tRNA_chargedAbundance_filter = df_tRNA_chargedAbundance %>% rownames_to_column(var = 'id')


tRNA_abundance_filter_forCodon = tRNA_abundance_filter %>% separate(id, into = c('AA','anticodon'), remove = F) %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample', values_to = 'count') %>% filter(AA!='iMet')
tRNA_abundance_filter_byCodon = tibble()
for(i in unique(anticodon_codon_matrix$Codon)){
  anti_i = anticodon_codon_matrix %>% filter(Codon==i) %>% dplyr::select(anticodon, efficiency)
  codon_eff_dict = anti_i$efficiency
  names(codon_eff_dict) = anti_i$anticodon
  cognate_tRNA = tRNA_abundance_filter_forCodon %>% filter(anticodon%in%anti_i$anticodon)
  cognate_tRNA = cognate_tRNA %>% mutate(Codon=i, efficiency=codon_eff_dict[anticodon])
  cognate_tRNA_toCodon = cognate_tRNA %>% group_by(AA,Codon,sample) %>% summarise(count = sum(efficiency*count), .groups = 'drop')
  tRNA_abundance_filter_byCodon = bind_rows(tRNA_abundance_filter_byCodon, cognate_tRNA_toCodon)
}
tRNA_abundance_filter_byCodon = tRNA_abundance_filter_byCodon %>% mutate(id = paste0(Codon,':',AA)) %>% pivot_wider(names_from = 'sample', values_from = 'count')

tRNA_chargedAbundance_filter_forCodon = tRNA_chargedAbundance_filter %>% separate(id, into = c('AA','anticodon'), remove = F) %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample', values_to = 'count') %>% filter(AA!='iMet')
tRNA_chargedAbundance_filter_byCodon = tibble()
for(i in unique(anticodon_codon_matrix$Codon)){
  anti_i = anticodon_codon_matrix %>% filter(Codon==i) %>% dplyr::select(anticodon, efficiency)
  codon_eff_dict = anti_i$efficiency
  names(codon_eff_dict) = anti_i$anticodon
  cognate_tRNA = tRNA_chargedAbundance_filter_forCodon %>% filter(anticodon%in%anti_i$anticodon)
  cognate_tRNA = cognate_tRNA %>% mutate(Codon=i, efficiency=codon_eff_dict[anticodon])
  cognate_tRNA_toCodon = cognate_tRNA %>% group_by(AA,Codon,sample) %>% summarise(count = sum(efficiency*count), .groups = 'drop')
  tRNA_chargedAbundance_filter_byCodon = bind_rows(tRNA_chargedAbundance_filter_byCodon, cognate_tRNA_toCodon)
}
tRNA_chargedAbundance_filter_byCodon = tRNA_chargedAbundance_filter_byCodon %>% mutate(id = paste0(Codon,':',AA)) %>% pivot_wider(names_from = 'sample', values_from = 'count')

tRNA_abundance_filter_byAA = tRNA_abundance_filter %>% separate(id, into = c('AA','anticodon'), remove = F) %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample') %>% group_by(sample,AA) %>% summarise(value = sum(value), .groups = 'drop') %>% pivot_wider(id_cols = AA, names_from = 'sample', values_from = 'value') %>% mutate(id=AA) %>% filter(AA!='SeC', AA!='iMet')
tRNA_chargedAbundance_filter_byAA = tRNA_chargedAbundance_filter %>% separate(id, into = c('AA','anticodon'), remove = F) %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample') %>% group_by(sample,AA) %>% summarise(value = sum(value), .groups = 'drop') %>% pivot_wider(id_cols = AA, names_from = 'sample', values_from = 'value') %>% mutate(id=AA) %>% filter(AA!='SeC', AA!='iMet')

mid1 = tRNA_chargedAbundance_filter_byCodon %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample', values_to = 'charged')
mid2 = tRNA_abundance_filter_byCodon %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample', values_to = 'total')
tRNA_chargeRate_byCodon = full_join(mid1, mid2, by = c('AA','Codon','id','sample')) %>% mutate(rate = charged/total, sample = gsub('XS','S',sample)) %>% pivot_wider(id_cols = id, names_from = 'sample', values_from = 'rate')

mid1 = tRNA_chargedAbundance_filter_byAA %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample', values_to = 'charged')
mid2 = tRNA_abundance_filter_byAA %>% pivot_longer(Ctrl_1:XS4_2, names_to = 'sample', values_to = 'total')
tRNA_chargeRate_byAA = full_join(mid1, mid2, by = c('AA','id','sample')) %>% mutate(rate = charged/total, sample = gsub('XS','S',sample)) %>% pivot_wider(id_cols = id, names_from = 'sample', values_from = 'rate')

isodecoder_abundance_raw = fread('Fig5/Results/tRNA_isodecoder_abundance.txt')
colnames(isodecoder_abundance_raw) = gsub('/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj/clean_','',colnames(isodecoder_abundance_raw))
colnames(isodecoder_abundance_raw) = gsub('_R1.fq.gz.unpaired_uniq.bam','',colnames(isodecoder_abundance_raw))
isodecoder_abundance_filter = isodecoder_abundance_raw %>% filter(!grepl('Escherichia_coli',isodecoder), !grepl('mito',isodecoder), !grepl('Homo_sapiens_tRX',isodecoder)) %>% mutate(gene = gsub('Homo_sapiens_tRNA-','',isodecoder), gene = gsub('/.*','',gene))
colnames(isodecoder_abundance_filter) = gsub('X','XS',colnames(isodecoder_abundance_filter))
isodecoder_abundance_filter = isodecoder_abundance_filter %>% filter(rowMeans(across(Ctrl_1:XS4_2))>10)
isodecoder_abundance_filter$id = isodecoder_abundance_filter$gene
df_isodecoder_abundance = isodecoder_abundance_filter %>% column_to_rownames(var = 'id') %>% dplyr::select(Ctrl_1:XS4_2)


isodecoder_chargeRate_bySample = fread('Fig5/Results/tRNA_isodecoder_ChargeRate.tsv')
colnames(isodecoder_chargeRate_bySample) = gsub('X','S',colnames(isodecoder_chargeRate_bySample))
isodecoder_chargeRate_bySample$id = isodecoder_chargeRate_bySample$tRNA_isodecoder
df_isodecoder_chargeRate_bySample = isodecoder_chargeRate_bySample %>% mutate(tRNA_isodecoder = gsub('-multi','',tRNA_isodecoder)) %>% column_to_rownames(var = 'tRNA_isodecoder') %>% dplyr::select(-id)
colnames(df_isodecoder_chargeRate_bySample) = gsub('S','XS',colnames(df_isodecoder_chargeRate_bySample))
charging_vaild_isodecoder_abundance = df_isodecoder_abundance[intersect(rownames(df_isodecoder_abundance), rownames(df_isodecoder_chargeRate_bySample)),] %>% dplyr::select(Ctrl_1:XS4_2)
charging_vaild_isodecoder_chargedAbundance = charging_vaild_isodecoder_abundance * df_isodecoder_chargeRate_bySample[rownames(charging_vaild_isodecoder_abundance),colnames(charging_vaild_isodecoder_abundance)]
charging_vaild_isodecoder_chargedAbundance = charging_vaild_isodecoder_chargedAbundance[which(rowMeans(charging_vaild_isodecoder_chargedAbundance)>10),]
charging_vaild_isodecoder_chargedAbundance_tibble = charging_vaild_isodecoder_chargedAbundance %>% rownames_to_column(var = 'id')

tRNA_chargeRate_bySample_long = df_tRNA_chargeRate_bySample %>% rownames_to_column(var = 'id') %>% pivot_longer(Ctrl_1:XS4_2, values_to = 'charging_rate', names_to = 'sample')
tRNA_abundance_filter_long = tRNA_abundance_filter %>% pivot_longer(Ctrl_1:XS4_2, values_to = 'abundance', names_to = 'sample') %>% dplyr::select(id, sample, abundance)
condition_int_trans = c('Ctrl'=0,'XS1'=1,'XS2'=2,'XS4'=3)
tRNA_charging_rate_total = left_join(tRNA_chargeRate_bySample_long, tRNA_abundance_filter_long, by = c('id','sample')) %>%
  separate(sample, into = c('condition','rep'), remove = F) %>% mutate(charged_abundance = charging_rate * abundance, condition_int = condition_int_trans[condition]) %>%
  group_by(sample, condition, condition_int) %>% summarise(abundance = sum(abundance), charged_abundance = sum(charged_abundance), .groups = 'drop') %>%
  mutate(charging_rate = charged_abundance/abundance)

summary(MASS::glm.nb(charged_abundance ~ condition_int + offset(log(abundance)), data = tRNA_charging_rate_total))
Fig5B = ggplot(tRNA_charging_rate_total, aes(x = condition, y = charging_rate)) + stat_summary(aes(color = condition), geom = 'point', pch=16, fun = mean, size = 2) +
  stat_summary(aes(group = 1), geom = 'line', fun = mean) + stat_summary(aes(color = condition), geom = 'errorbar', fun.data = mean_se, width = 0.1) + theme_pubr() +
  labs(x=NULL, y='Global tRNA charging rate') + guides(color = 'none') + coord_cartesian(ylim = c(0.65,0.8)) +
  scale_color_manual(values =  c(
    "Ctrl" = "#575959",
    "XS1"  = "#E08251",
    "XS2"  = "#9D71A6",
    "XS4"  = "#4088B4")) +
  scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

plot_charging_rate_byAA = get_charging_rate_plot_2(tRNA_chargeRate_byAA, tRNA_abundance_filter_byAA, tRNA_chargedAbundance_filter_byAA)
Fig5C = plot_charging_rate_byAA[[3]]

plot_charging_rate_byCodon = get_charging_rate_plot_2(tRNA_chargeRate_byCodon, tRNA_abundance_filter_byCodon, tRNA_chargedAbundance_filter_byCodon)
Fig5D = plot_charging_rate_byCodon[[3]]

plot_charging_rate_byIsoacceptor = get_charging_rate_plot_2(tRNA_chargeRate_bySample, tRNA_abundance_filter, tRNA_chargedAbundance_filter)
FigS22A = plot_charging_rate_byIsoacceptor[[3]]

plot_charging_rate_byIsodecoder = get_charging_rate_plot_2(isodecoder_chargeRate_bySample%>%filter(id%in%charging_vaild_isodecoder_chargedAbundance_tibble$id), isodecoder_abundance_filter%>%filter(id%in%charging_vaild_isodecoder_chargedAbundance_tibble$id), charging_vaild_isodecoder_chargedAbundance_tibble)
FigS22B = plot_charging_rate_byIsodecoder[[3]]

# Fig5K, FigS23C
Codon_stats_extended=fread('Fig5/Results/Codon_stats_extended.csv')
relative_SR = Codon_stats_extended %>% group_by(Codon) %>% mutate(relative_ASR = A_site_stalling_rate/A_site_stalling_rate[condition=='Ctrl'], relative_PSR = P_site_stalling_rate/P_site_stalling_rate[condition=='Ctrl']) %>% ungroup()

Fig5K = ggplot(relative_SR, aes(x=condition, y=relative_ASR)) + geom_line(aes(group = Codon, color = (Codon=='TTA')), linewidth = 0.8, alpha = 0.5) +
  geom_point(aes(color = condition), pch = 16, size = 2) + theme_pubr() +
  guides(color='none') + coord_cartesian(ylim=c(0.8, 2.5)) + labs(x=NULL, y = 'Normalized A-site stalling rate') +
  scale_color_manual(values =  c(
    "Ctrl" = "#575959",
    "XS1"  = "#E08251",
    "XS2"  = "#9D71A6",
    "XS4"  = "#4088B4",
    "TRUE" = "#e6211a",
    "FALSE" = "grey")) +
  scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

FigS23C = ggplot(relative_SR, aes(x=condition, y=relative_PSR)) + geom_line(aes(group = Codon, color = (Codon=='TTA')), linewidth = 0.8, alpha = 0.5) +
  geom_point(aes(color = condition), pch = 16, size = 2) + theme_pubr() +
  guides(color='none') + coord_cartesian(ylim=c(0.6, 2.5)) + labs(x=NULL, y = 'Normalized P-site stalling rate') +
  scale_color_manual(values =  c(
    "Ctrl" = "#575959",
    "XS1"  = "#E08251",
    "XS2"  = "#9D71A6",
    "XS4"  = "#4088B4",
    "TRUE" = "#e6211a",
    "FALSE" = "grey")) +
  scale_x_discrete(labels = c('Ctrl'='Ctrl','XS1'='S1','XS2'='S2','XS4'='S4'))

# FigS23AB


ASR_rate = Codon_stats_extended %>% dplyr::select(Codon, AA3, condition, A_site_stalling_rate) %>% mutate(id=paste0(Codon,':',AA3), condition = condition_trans[condition]) %>% pivot_wider(id_cols = id, names_from = condition, values_from = A_site_stalling_rate)
ASR_total = Codon_stats_extended %>% dplyr::select(Codon, AA3, condition, A_site_total_num) %>% mutate(id=paste0(Codon,':',AA3)) %>% pivot_wider(id_cols = id, names_from = condition, values_from = A_site_total_num)
ASR_stalled = Codon_stats_extended %>% dplyr::select(Codon, AA3, condition, A_site_stalled_num) %>% mutate(id=paste0(Codon,':',AA3)) %>% pivot_wider(id_cols = id, names_from = condition, values_from = A_site_stalled_num)
heatmap_ASR = get_rate_plot_2_byCondition(ASR_rate, ASR_total, ASR_stalled)
FigS23A = heatmap_ASR[[3]]

PSR_rate = Codon_stats_extended %>% dplyr::select(Codon, AA3, condition, P_site_stalling_rate) %>% mutate(id=paste0(Codon,':',AA3), condition = condition_trans[condition]) %>% pivot_wider(id_cols = id, names_from = condition, values_from = P_site_stalling_rate)
PSR_total = Codon_stats_extended %>% dplyr::select(Codon, AA3, condition, P_site_total_num) %>% mutate(id=paste0(Codon,':',AA3)) %>% pivot_wider(id_cols = id, names_from = condition, values_from = P_site_total_num)
PSR_stalled = Codon_stats_extended %>% dplyr::select(Codon, AA3, condition, P_site_stalled_num) %>% mutate(id=paste0(Codon,':',AA3)) %>% pivot_wider(id_cols = id, names_from = condition, values_from = P_site_stalled_num)
headmat_PSR = get_rate_plot_2_byCondition(PSR_rate, PSR_total, PSR_stalled)
FigS23B = headmat_PSR[[3]]





Fig5A
Fig5B
Fig5C
Fig5D
Fig5K

FigS22A
FigS22B

FigS23A
FigS23B
FigS23C
