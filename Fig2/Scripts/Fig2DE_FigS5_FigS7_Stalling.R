#######################################################################################################

#global

#######################################################################################################
setwd("C:/Users/Ke-jia Shan/Desktop/Virus_CUB/Fig2/Results/")
Asite=fread("Asite.txt")
#vTDI=virus/Total
ggplot(Asite, 
       aes(x = virus/Total, y = A_Peak/Asite_total/Ctrl_Astall)) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S)) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col = "red", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "s", col = "red", size = 5,label.y = 7.5 , cor.coef.name = "rho") +  # 所有点统一算相关系数
  ggpubr::stat_cor(method = "p", col = "red", size = 5,label.y = 9, cor.coef.name = "R") +  # 所有点统一算相关系数
  theme_classic() +
  scale_x_continuous(breaks = c(seq(0,1,0.2)),limits = c(0,1)) +
  labs(x = "vTBI", y = "Relative stalling rate")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/vTBI_stall_all_relative.pdf', width= 6.07 , height= 4.05 , units='in')


#######################################################################################################

#codon level

#######################################################################################################
Asite=fread("Asite_Codon.txt")

head(Asite)
ggplot(Asite, 
       aes(x = virus/(virus+Host), y = (Asite_Peak/Asite_total)/(Ctrl_Astall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S)) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  labs(x = "vTBI", y = "Fold change of A-site stalling rate")+
  #facet_wrap(~Sample_standard,scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank()#,
    #strip.text = element_blank()  # 文字也隐藏
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_Ctrl.pdf'), width= 6.07 , height= 4.05 , units='in')


ggplot(filter(Asite,Codon=="TTA"), 
       aes(x = virus/(virus+Host), y = (Asite_Peak/Asite_total)/(Ctrl_Astall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=6) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 15,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  labs(x = "vTBI", y = "Fold change of A-site stalling rate")+facet_wrap(~paste(Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_TTA.pdf'), width= 10, height= 6.1959  , units='in')

ggplot(filter(Asite,Codon=="TTA"), 
       aes(x = virus/(virus+Host), y = (Psite_Peak/Psite_total)/(Ctrl_Pstall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=6) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 15,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  labs(x = "vTBI", y = "Fold change of P-site stalling rate")+facet_wrap(~paste(Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_TTA.pdf'), width= 6.3739, height= 6.1959  , units='in')


ggplot(filter(Asite,Codon=="AGT"), 
       aes(x = virus/(virus+Host), y = (Asite_Peak/Asite_total)/(Ctrl_Astall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=6) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 18,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  scale_y_continuous(breaks = seq(0,18,3),limits = c(0,18)) +
  labs(x = "vTBI", y = "Fold change of A-site stalling rate")+facet_wrap(~paste(Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_AGT.pdf'), width= 6.3739, height= 6.1959  , units='in')


ggplot(filter(Asite,Codon=="TCA"), 
       aes(x = virus/(virus+Host), y = (Psite_Peak/Psite_total)/(Ctrl_Pstall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=6) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 18,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  scale_y_continuous(breaks = seq(0,18,3),limits = c(0,18)) +
  scale_x_continuous(breaks = seq(0,1,0.25),limits = c(0,1)) +
  labs(x = "vTBI", y = "Fold change of P-site stalling rate")+facet_wrap(~paste(Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_TCA.pdf'), width= 6.3739, height= 6.1959  , units='in')

ggplot(filter(Asite,Codon=="AAT"), 
       aes(x = virus/(virus+Host), y = (Asite_Peak/Asite_total)/(Ctrl_Astall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=6) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 11.52,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  scale_x_continuous(limits = c(0,1)) +
  labs(x = "vTBI", y = "Fold change of A-site stalling rate")+facet_wrap(~paste(Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_GAA.pdf'), width= 6.3739, height= 6.1959  , units='in')


ggplot(filter(Asite,Codon=="AAT"), 
       aes(x = virus/(virus+Host), y = (Psite_Peak/Psite_total)/(Ctrl_Pstall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=6) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 11.52,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  scale_x_continuous(limits = c(0,1)) +
  labs(x = "vTBI", y = "Fold change of P-site stalling rate")+facet_wrap(~paste(Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_GAA.pdf'), width= 6.3739, height= 6.1959  , units='in')

#======================================================================
# 1）先生成相关系数表格
#======================================================================
cor_table <- Asite %>%
  group_by(AA, Codon) %>%
  summarise(
    spearman_rho= cor.test(virus/(virus+Host), (Asite_Peak/Asite_total)/(Ctrl_Astall), method="spearman")$estimate,
    spearman_P = cor.test(virus/(virus+Host), (Asite_Peak/Asite_total)/(Ctrl_Astall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA, ":", Codon))

# 显著性 *
cor_table <- cor_table %>%
  mutate(

    p_star_spearman = case_when(
      spearman_P < 0.001 ~ "***",
      spearman_P < 0.01 ~ "**",
      spearman_P < 0.05 ~ "*", TRUE ~ ""
    )
  )


#======================================================================
# 2）排序规则：先按 AA 排序 → 每个 AA 内按 R 从大到小
#======================================================================

codon_table<-fread("C:/Users/Ke-jia Shan/Desktop/Virus_CUB/Fig1/Scripts//codon_table.txt",header=F)
colnames(codon_table)=c("Codon","AA3","AA","AAname")
cor_table=merge(cor_table,codon_table,by=c("AA","Codon"))
cor_table$AA_Codon=paste0(cor_table$AA3,":",cor_table$Codon)
# ===================== 你给的严格配色 =====================
aa_color <- data.frame(
  AA = c("L","S","T","V","N","A","G","I","F","K","Q","D","P","Y","E","R","C","H","M","W"),
  AA3 = c("Leu","Ser","Thr","Val","Asn","Ala","Gly","Ile","Phe","Lys","Gln","Asp","Pro","Tyr","Glu","Arg","Cys","His","Met","Trp"),
  color = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#CCEBC5","#FFED6F","#B2DF8A")
)

# 生成 AA3 颜色映射（按 AA3 着色，核心修改）
color_map_AA3 <- setNames(aa_color$color, aa_color$AA3)

# # ===================== 读取 & 构建排序依据 =====================
# mid <- fread("HEK293FT_RSCU.txt")
# 
# # 排序：先 AA3 → 再 RSCU 降序（你要求的排序）
# mid_sorted <- mid %>%
#   mutate(
#     Codon3 = substr(Codon, 3, 3),
#     Codon3_group = ifelse(Codon3 %in% c("C", "G"), "CG", "AT")  # 第3位 CG/AT
#   ) %>%
#   arrange(
#     AA3,                  # 第一排序：氨基酸
#     Codon3_group,         # 第二排序：CG 在前，AT 在后
#     desc(RSCU)            # 第三排序：RSCU 从大到小
#   )
# # 构建 mid1 = AA3:Codon
# mid1 <- paste0(mid_sorted$AA3, ":", mid_sorted$Codon)
# 

# ===================== 绘图：Spearman（按 AA3 着色 + 固定顺序） =====================
p_spearman <- ggplot(cor_table, aes(x = AA_Codon, y = spearman_rho, fill = AA3)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_spearman), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map_AA3) +  # 用 AA3 配色
  labs(x = "AA:Codon", y = "Spearman's rho\n(vTBI vs. fold change of A-site stalling rate)") +
  #scale_y_continuous(breaks = c(seq(0,0.6,0.1)),limits = c(0,0.601)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none")

print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_allSample_eachcodon_spearman.pdf'), width= 14 , height= 3 , units='in')



#======================================================================
# 先生成P-site 相关系数表格
#======================================================================
cor_table <- Asite %>%
  group_by(AA, Codon) %>%
  summarise(

    spearman_rho= cor.test(virus/(virus+Host), (Psite_Peak/Psite_total)/(Ctrl_Pstall), method="spearman")$estimate,
    spearman_P = cor.test(virus/(virus+Host), (Psite_Peak/Psite_total)/(Ctrl_Pstall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA, ":", Codon))

# 显著性 *
cor_table <- cor_table %>%
  mutate(

    p_star_spearman = case_when(
      spearman_P < 0.001 ~ "***",
      spearman_P < 0.01 ~ "**",
      spearman_P < 0.05 ~ "*", TRUE ~ ""
    )
  )


#======================================================================
# 2）排序规则：先按 AA 排序 → 每个 AA 内按 R 从大到小
#======================================================================

codon_table<-fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/UMAP/codon_table.txt",header=F)
colnames(codon_table)=c("Codon","AA3","AA","AAname")
cor_table=merge(cor_table,codon_table,by=c("AA","Codon"))
cor_table$AA_Codon=paste0(cor_table$AA3,":",cor_table$Codon)
# ===================== 你给的严格配色 =====================
aa_color <- data.frame(
  AA = c("L","S","T","V","N","A","G","I","F","K","Q","D","P","Y","E","R","C","H","M","W"),
  AA3 = c("Leu","Ser","Thr","Val","Asn","Ala","Gly","Ile","Phe","Lys","Gln","Asp","Pro","Tyr","Glu","Arg","Cys","His","Met","Trp"),
  color = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#CCEBC5","#FFED6F","#B2DF8A")
)

# 生成 AA3 颜色映射（按 AA3 着色，核心修改）
color_map_AA3 <- setNames(aa_color$color, aa_color$AA3)


# ===================== 给 cor_table 设置横坐标因子顺序（关键！） =====================
# 让 ggplot 横坐标严格按 mid1 排序
cor_table$AA_Codon <- factor(cor_table$AA_Codon, levels = mid1)



# ===================== 绘图：Spearman（按 AA3 着色 + 固定顺序） =====================
p_spearman <- ggplot(cor_table, aes(x = AA_Codon, y = spearman_rho, fill = AA3)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_spearman), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map_AA3) +  # 用 AA3 配色
  labs(x = "AA:Codon", y = "Spearman's rho\n(vTBI vs. fold change of P-site stalling rate)") +
  scale_y_continuous(breaks = c(seq(0,0.6,0.1)),limits = c(0,0.601)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none")


print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_allSample_eachcodon_spearman.pdf'), width= 14 , height= 3 , units='in')


# ==========================================

# AA

# ==========================================
Asite=fread("Asite_Codon.txt")
AA=Asite %>% 
  group_by(Sample_standard,Study,AA,virus_S) %>% 
  dplyr::summarise(AA_Peak=sum(Asite_Peak),AA_total=sum(Asite_total),AA_virus=sum(virus),AA_Host=sum(Host))

head(AA)
ctrl <- filter(Asite, grepl("Ctrl|mock|inf|00hpi", Sample_standard, ignore.case = TRUE)) %>% 
  group_by(Sample_standard,Study,AA,virus_S) %>% 
  dplyr::summarise(AA_Peak=sum(Asite_Peak),AA_total=sum(Asite_total),AA_virus=sum(virus),AA_Host=sum(Host))

ctrl=ctrl %>% group_by(Study,AA) %>% dplyr::summarise(Ctrl_Astall=mean(AA_Peak/AA_total))
Asite=merge(AA,ctrl,by=c("Study","AA"))

head(Asite)



ggplot(Asite, 
       aes(x = AA_virus/(AA_virus+AA_Host), y = (AA_Peak/AA_total)/(Ctrl_Astall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S)) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 25,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 20,size=5, cor.coef.name = "R") +  
  theme_classic() +
  labs(x = "vTBI", y = "Relative stalling rate")+#facet_wrap(~paste(AA,":",Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank()#,
    #strip.text = element_blank()  # 文字也隐藏
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_AA_Ctrl.pdf'), width= 6.07 , height= 4.05 , units='in')



ggplot(filter(Asite,AA=="V"), 
       aes(x = AA_virus/(AA_virus+AA_Host), y = (AA_Peak/AA_total)/(Ctrl_Astall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=5) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",label.y = 10.5,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 20,size=5, cor.coef.name = "R") +  
  theme_classic() +
  scale_x_continuous(breaks = c(seq(0,1,0.25)),limits = c(0,1)) +
  labs(x = "vTBI", y = "Relative stalling rate")+#facet_wrap(~paste(AA,":",Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
    #strip.text = element_blank()  # 文字也隐藏
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_AA_V.pdf'), width= 7.2468 , height= 6.902 , units='in')

head(Asite)
Asite$virus_stall=Asite$AA_Peak/Asite$AA_total
#colnames(Asite)[5:8]=c("Asite_Peak","Asite_total","virus","Host")


cor_table <- Asite %>%
  group_by(AA) %>%
  summarise(
   spearman_rho= cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Astall), method="spearman")$estimate,
    spearman_P = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Astall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA))

# 显著性 *
cor_table <- cor_table %>%
  mutate(

    p_star_spearman = case_when(
      spearman_P < 0.001 ~ "***",
      spearman_P < 0.01 ~ "**",
      spearman_P < 0.05 ~ "*", TRUE ~ ""
    )
  )
cor_table=merge(cor_table,unique(codon_table[,c("AA","AA3")]),by=c("AA"))

#======================================================================
# 你指定的 20 种氨基酸标准配色（已替换）
#======================================================================
aa_color <- data.frame(
  AA  = c("L","S","T","V","N","A","G","I","F","K","Q","D","P","Y","E","R","C","H","M","W"),
  AA3 = c("Leu","Ser","Thr","Val","Asn","Ala","Gly","Ile","Phe","Lys","Gln","Asp","Pro","Tyr","Glu","Arg","Cys","His","Met","Trp"),
  color = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#CCEBC5","#FFED6F","#B2DF8A")
)

# 生成 AA 颜色映射（关键：按单字母 AA 匹配）
color_map <- setNames(aa_color$color, aa_color$AA)

#======================================================================
#Spearman Rho（使用你指定的氨基酸配色）
#======================================================================
p_spearman <- ggplot(cor_table, aes(x = AA3, y = spearman_rho, fill = AA)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_spearman), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map) +  # 已替换为你的配色
  labs(x = "AA", y = "Spearman's rho\n(vTBI vs. fold change of stalling rate)") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), limits = c(0, 0.601)) +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "none")

print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_relative_eachAA_spearman.pdf'), width= 6 , height= 3 , units='in')


# ==========================================

# AA

# ==========================================
Asite=fread("Asite_Codon.txt")
AA=Asite %>% 
  group_by(Sample_standard,Study,AA,virus_S) %>% 
  dplyr::summarise(AA_Peak=sum(Psite_Peak),AA_total=sum(Psite_total),AA_virus=sum(virus),AA_Host=sum(Host))

head(AA)
ctrl <- filter(Asite, grepl("Ctrl|mock|inf|00hpi", Sample_standard, ignore.case = TRUE)) %>% 
  group_by(Sample_standard,Study,AA,virus_S) %>% 
  dplyr::summarise(AA_Peak=sum(Psite_Peak),AA_total=sum(Psite_total),AA_virus=sum(virus),AA_Host=sum(Host))

ctrl=ctrl %>% group_by(Study,AA) %>% dplyr::summarise(Ctrl_Pstall=mean(AA_Peak/AA_total))
Asite=merge(AA,ctrl,by=c("Study","AA"))

head(Asite)



ggplot(Asite, 
       aes(x = AA_virus/(AA_virus+AA_Host), y = (AA_Peak/AA_total)/(Ctrl_Pstall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S)) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",,label.y = 25,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 20,size=5, cor.coef.name = "R") +  
  theme_classic() +
  labs(x = "vTBI", y = "Fold change of stalling\nrate relative Ctrl")+#facet_wrap(~paste(AA,":",Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank()#,
    #strip.text = element_blank()  # 文字也隐藏
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_AA_Ctrl.pdf'), width= 6.07 , height= 4.05 , units='in')


ggplot(filter(Asite,AA=="V"), 
       aes(x = AA_virus/(AA_virus+AA_Host), y = (AA_Peak/AA_total)/(Ctrl_Pstall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S),size=5) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",label.y = 11.5,size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 20,size=5, cor.coef.name = "R") +  
  theme_classic() +
  scale_x_continuous(breaks = c(seq(0,1,0.25)),limits = c(0,1)) +
  labs(x = "vTBI", y = "Relative stalling rate")+#facet_wrap(~paste(AA,":",Codon),scales = "free") +
  scale_color_manual(values=c('#EB746A', 
                              #'#D08D0F', 
                              '#9EA020', '#43AE35', '#2CB178',
                              '#1EB5B8', '#2BA5DF', 
                              #'#8B8BC3', 
                              '#B274AF', '#DF66A2'))+
  theme(
    strip.background = element_blank(),
    legend.position = "none"
    #strip.text = element_blank()  # 文字也隐藏
  )
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_AA_V.pdf'), width= 7.2468 , height= 6.902 , units='in')

head(Asite)

cor_table <- Asite %>%
  group_by(AA) %>%
  summarise(
    spearman_rho= cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Pstall), method="spearman")$estimate,
    spearman_P = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Pstall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA))

# 显著性 *
cor_table <- cor_table %>%
  mutate(

    p_star_spearman = case_when(
      spearman_P < 0.001 ~ "***",
      spearman_P < 0.01 ~ "**",
      spearman_P < 0.05 ~ "*", TRUE ~ ""
    )
  )


cor_table=merge(cor_table,unique(codon_table[,c("AA","AA3")]),by=c("AA"))

#======================================================================
# 你指定的 20 种氨基酸标准配色（已替换）
#======================================================================
aa_color <- data.frame(
  AA  = c("L","S","T","V","N","A","G","I","F","K","Q","D","P","Y","E","R","C","H","M","W"),
  AA3 = c("Leu","Ser","Thr","Val","Asn","Ala","Gly","Ile","Phe","Lys","Gln","Asp","Pro","Tyr","Glu","Arg","Cys","His","Met","Trp"),
  color = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#CCEBC5","#FFED6F","#B2DF8A")
)

# 生成 AA 颜色映射（关键：按单字母 AA 匹配）
color_map <- setNames(aa_color$color, aa_color$AA)


#======================================================================
# Spearman Rho（使用你指定的氨基酸配色）
#======================================================================
p_spearman <- ggplot(cor_table, aes(x = AA3, y = spearman_rho, fill = AA)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_spearman), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map) +  # 已替换为你的配色
  labs(x = "AA", y = "Spearman's rho\n(vTBI vs. fold change of stalling rate)") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 0.6, 0.1), limits = c(0, 0.601)) +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "none")

print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_relative_eachAA_spearman.pdf'), width= 6 , height= 3 , units='in')

