codon_table<-fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/UMAP/codon_table.txt",header=F)
colnames(codon_table)=c("Codon","AA3","AA","AAname")
Ribo=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Ribo_codon.txt")
Ribo=merge(Ribo,codon_table,by=c("Codon"))
head(Ribo)
Ribo=Ribo %>% group_by(Study,AA3,Sample) %>% 
  dplyr::summarise(Virus_AA=sum(virus))  %>% 
  filter(!grepl("Ctrl|mock|inf", Sample, ignore.case = TRUE))
Ribo


Ribo=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Ribo_TBI.txt")
Asite=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/AsiteStall.csv")
Asite=merge(Asite,Ribo,by.x=c("Study","sample"),by.y=c("Study","Sample"))
#write.table(Asite,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Asite.txt",row.names = F,quote = F,sep = "\t")
Asite[Asite$Study=="SCV2_Vero",]$Study="SCV2_zzVero"
Asite[Asite$Study=="SCV2_Calu3",]$Study="SCV2_zCalu3"
Asite$virus_S="virus"
Asite[grepl("EToV",Study)]$virus_S="EToV"
Asite[grepl("HCMV",Study)]$virus_S="HCMV"
Asite[grepl("IAV",Study)]$virus_S="IAV"
Asite[grepl("MHV",Study)]$virus_S="MHV"
Asite[grepl("SCV2",Study)]$virus_S="SCV2"
Asite[grepl("VACV",Study)]$virus_S="VACV"
#去掉SCV2 uninf NC，不过加上之后pearson也显著，影响不大，就是看起来uninf有些异常
head(Asite)
ggplot(Asite,aes(x=virus/Total,y=A_Peak/Asite_total))+
  geom_point(aes(colour = Study, shape = virus_S)) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm",col="red",se=F)+
  ggpubr::stat_cor(method = "s", col = "red", size = 5,label.y = 0.15 , cor.coef.name = "rho") +  # 所有点统一算相关系数
  ggpubr::stat_cor(method = "p", col = "red", size = 5,label.y = 0.13, cor.coef.name = "R") +  # 所有点统一算相关系数
  theme_classic()+
  scale_x_continuous(breaks = c(seq(0,1,0.2)),limits = c(0,1)) +
  #scale_y_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  labs(x="vTBI",y="Stalling rate")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/vTBI_stall_all_zhongti.pdf', width= 6.07 , height= 4.05 , units='in')


unique(Asite$Study)

Asite=filter(Asite,!(grepl("inf", sample, ignore.case = TRUE) & Study=="SCV2_PMID_34433827"))
ctrl <- filter(Asite, grepl("Ctrl|mock|inf|00hpi", sample, ignore.case = TRUE))
#ctrl <- filter(Asite, grepl("Ctrl|mock|inf", sample, ignore.case = TRUE))

unique(ctrl$Study)
ctrl=ctrl %>% group_by(Study) %>% dplyr::summarise(Ctrl_Astall=mean(P_Peak/Psite_total))
head(ctrl)
head(Asite)
Asite=merge(Asite,ctrl[,c(1,2)],by=c("Study"))
Asite$virus_stall=Asite$A_Peak/Asite$Asite_total
Asite$virus_S="virus"
Asite[grepl("EToV",Study)]$virus_S="EToV"
Asite[grepl("HCMV",Study)]$virus_S="HCMV"
Asite[grepl("IAV",Study)]$virus_S="IAV"
Asite[grepl("MHV",Study)]$virus_S="MHV"
Asite[grepl("SCV2",Study)]$virus_S="SCV2"
Asite[grepl("VACV",Study)]$virus_S="VACV"
write.table(Asite,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Asite.txt",row.names = F,quote = F,sep = "\t")
Asite=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Asite.txt")
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

ggplot(filter(Asite,virus>0,sample!="SCV2_moi10_00hpi"), 
       aes(x = virus/Total, y = A_Peak/Asite_total/Ctrl_Astall)) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S)) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col = "red", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "s", col = "red", size = 5,label.y = 7.5 , cor.coef.name = "rho") +  # 所有点统一算相关系数
  ggpubr::stat_cor(method = "p", col = "red", size = 5,label.y = 9, cor.coef.name = "R") +  # 所有点统一算相关系数
  theme_classic() +
  scale_x_continuous(breaks = c(seq(0,1,0.2)),limits = c(0,1)) +
  labs(x = "vTBI", y = "Relative stalling rate")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/vTBI_stall_all_relative2.pdf', width= 6.07 , height= 4.05 , units='in')

#######################################################################################################

#codon level

#######################################################################################################
setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Codon/")
temp <- list.files(path = "./",pattern = "*.csv")
Asite=data.frame()
i=1
head(file)
for(i in 1:length(temp)){
  #temp[i]
  name=temp[i]
  name=gsub(".Peak_number.csv","",name)
  file=fread(temp[i])
  file$Study=name
  Asite=rbind(Asite,file)
  
}
head(Asite)
unique(Asite$Study)

# 核心绘图：分面 + TTA标红 + 每个Study一张图
unique(Asite$Sample)
Asite[Asite$Sample=="HMCV_5hpi",]$Sample="HMCV_05hpi"
Asite[Asite$Sample=="zVACV_uninf",]$Sample="zVACV_00uninf"
Asite[Asite$Sample=="zVACV_Early",]$Sample="zVACV_0Early"
Asite[Asite$Sample=="zVACV_uninf",]$Sample="zVACV_00uninf"
Asite[Asite$Sample=="calu3_moi10_1hpi",]$Sample="calu3_moi10_01hpi"
Asite[Asite$Sample=="calu3_moi10_2hpi",]$Sample="calu3_moi10_02hpi"
Asite[Asite$Sample=="calu3_moi10_4hpi",]$Sample="calu3_moi10_04hpi"
Asite[Asite$Sample=="calu3_moi10_0hpi",]$Sample="calu3_moi10_00hpi"

# ===================== 加载包 =====================
library(dplyr)
library(stringr)

# ===================== 第一步：建立【Study → 标准化Sample名称】的对应规则 =====================
# 从 Ribo 提取 【Study + 标准Sample】 作为标准答案
standard <- unique(data.frame(Ribo)[,c("Sample","Study")])
colnames(standard)[1] <- "Sample_standard"

# 从 Asite 提取 【Study + 原始Sample】
asite_map <- unique(data.frame(Asite)[,c("Sample","Study")])
colnames(asite_map)[1] <- "Sample_asite"

# 全自动匹配：相同Study，自动对齐时间点（hpi/uninf/mock等）
asite_map <- asite_map %>%
  rowwise() %>%
  mutate(
    # 提取关键时间特征：1hpi / 2hpi / 3hpi /4hpi /5hpi /8hpi /24hpi /48hpi /72hpi /uninf /mock /Early
    key = str_extract(Sample_asite, "[0-9.]*hpi|uninf|mock|Early|00hpi|01hpi|02hpi|04hpi|12hpi|16hpi|24hpi|36hpi"),
    
    # 如果没提取到，用原始字符
    key = ifelse(is.na(key), tolower(Sample_asite), key),
    
    # 在标准Ribo表里找到【同一个Study + 同一个时间key】的标准名称
    Sample_standard = standard$Sample_standard[
      which(standard$Study == Study & 
              grepl(key, standard$Sample_standard, ignore.case = T))[1]
    ]
  ) %>%
  ungroup()

# 查看匹配结果（检查是否正确）
asite_map
asite_map[2,4]="Ctrl"
# ===================== 第二步：把 Asite 的 Sample 全部替换成 Ribo 标准格式 =====================
Asite=merge(Asite,asite_map,by.x=c("Sample","Study"),by.y=c("Sample_asite","Study"))
head(Asite)
Asite %>% group_by(Sample,Study) %>% dplyr::summarise(P=sum(Psite_Peak),P_total=sum(Psite_total))

head(Ribo)
Asite$virus_S="virus"
Asite[grepl("EToV",Study)]$virus_S="EToV"
Asite[grepl("HCMV",Study)]$virus_S="HCMV"
Asite[grepl("IAV",Study)]$virus_S="IAV"
Asite[grepl("MHV",Study)]$virus_S="MHV"
Asite[grepl("SCV2",Study)]$virus_S="SCV2"
Asite[grepl("VACV",Study)]$virus_S="VACV"
colnames(Asite)
write.table(Asite,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Asite_Codon.txt",row.names = F,quote = F,sep = "\t")
Asite=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Asite_Codon.txt")

Ribo=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Ribo_codon.txt")
head(Ribo)
unique(Ribo$Study)

Asite2=merge(Asite[,-c(12,13)],Ribo,by.x=c("Sample_standard","Study","Codon"),by.y=c("Sample","Study","Codon"))
Asite2[Asite2$Study=="SCV2_Vero",]$Study="SCV2_zzVero"
Asite2[Asite2$Study=="SCV2_Calu3",]$Study="SCV2_zCalu3"


# ==========================================

# relative codon

# ==========================================
Asite2=filter(Asite2,!(grepl("inf", Sample_standard, ignore.case = TRUE) & Study=="SCV2_PMID_34433827"))
ctrl <- filter(Asite2, grepl("Ctrl|mock|inf|00hpi", Sample_standard, ignore.case = TRUE))
#ctrl <- filter(Asite2, grepl("Ctrl|mock|inf", Sample_standard, ignore.case = TRUE))

unique(ctrl$Study)
ctrl=ctrl %>% group_by(Study,Codon) %>% dplyr::summarise(Ctrl_Astall=mean(Asite_Peak/Asite_total),
                                                         Ctrl_Pstall=mean(Psite_Peak/Psite_total),
                                                         Ctrl_Astall_AA=mean(mean_A_Peak_Frac_AA))
head(ctrl)
head(Asite2)
Asite=merge(Asite2,ctrl,by=c("Study","Codon"))
Asite$virus_stall=Asite$A_Peak/Asite$Asite_total

head(Asite)

write.table(Asite,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Asite_Codon.txt",row.names = F,quote = F,sep = "\t")
Asite=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/Asite_Codon.txt")
Asite=filter(Asite,!(grepl("inf", Sample_standard, ignore.case = TRUE) & Study=="SCV2_PMID_34433827"))
#Asite=filter(Asite,!(grepl("00hpi", Sample_standard, ignore.case = TRUE) & Study=="SCV2_PMID_34433827"))

head(Asite)
ggplot(filter(Asite,virus>=0), 
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

head(Asite)
Asite$vTBI=Asite$virus/(Asite$virus+Asite$Host)
mid=lm((P_Peak_Frac/Ctrl_Pstall)~vTBI*Codon,data.frame(Asite))
str(Asite)
Asite$Codon=factor(Asite$Codon)
mid=lm((A_Peak_Frac/Ctrl_Astall)~vTBI*Codon,data.frame(Asite))
summary(mid)

ggplot(filter(Asite,virus>=0), 
       aes(x = virus/(virus+Host), y = (Psite_Peak/Psite_total)/(Ctrl_Pstall))) +  # 全局去掉 colour/shape
  geom_point(aes(colour = Study, shape = virus_S)) +  # 只在点上分组着色、形状
  geom_smooth(method = "glm", col="black", se = FALSE) +  # 所有点统一拟合
  ggpubr::stat_cor(method = "spearman",size=5, cor.coef.name = "rho") +
  #ggpubr::stat_cor(method = "pearson",,label.y = 30,size=5, cor.coef.name = "R") +  
  theme_classic() +
  labs(x = "vTBI", y = "Fold change of P-site stalling rate")+
  #facet_wrap(~Study,scales = "free") +
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
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_Ctrl.pdf'), width= 6.07 , height= 4.05 , units='in')

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
# 1）先生成相关系数表格（你之前的修复版）
#======================================================================
cor_table <- Asite %>%
  group_by(AA, Codon) %>%
  summarise(
    pearson_R   = cor.test(virus/(virus+Host), (Asite_Peak/Asite_total)-(Ctrl_Astall), method="pearson")$estimate,
    pearson_P   = cor.test(virus/(virus+Host), (Asite_Peak/Asite_total)-(Ctrl_Astall), method="pearson")$p.value,
    spearman_rho= cor.test(virus/(virus+Host), (Asite_Peak/Asite_total)-(Ctrl_Astall), method="spearman")$estimate,
    spearman_P = cor.test(virus/(virus+Host), (Asite_Peak/Asite_total)-(Ctrl_Astall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA, ":", Codon))

# 显著性 *
cor_table <- cor_table %>%
  mutate(
    p_star_pearson = case_when(
      pearson_P < 0.001 ~ "***",
      pearson_P < 0.01 ~ "**",
      pearson_P < 0.05 ~ "*", TRUE ~ ""
    ),
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

# ===================== 读取 & 构建排序依据 =====================
mid <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/HEK293FT_RSCU.txt")

# 排序：先 AA3 → 再 RSCU 降序（你要求的排序）
mid_sorted <- mid %>%
  mutate(
    Codon3 = substr(Codon, 3, 3),
    Codon3_group = ifelse(Codon3 %in% c("C", "G"), "CG", "AT")  # 第3位 CG/AT
  ) %>%
  arrange(
    AA3,                  # 第一排序：氨基酸
    Codon3_group,         # 第二排序：CG 在前，AT 在后
    desc(RSCU)            # 第三排序：RSCU 从大到小
  )
# 构建 mid1 = AA3:Codon
mid1 <- paste0(mid_sorted$AA3, ":", mid_sorted$Codon)

# ===================== 给 cor_table 设置横坐标因子顺序（关键！） =====================
# 让 ggplot 横坐标严格按 mid1 排序
cor_table$AA_Codon <- factor(cor_table$AA_Codon, levels = mid1)

# ===================== 绘图：Pearson（按 AA3 着色 + 固定顺序） =====================
p_pearson <- ggplot(cor_table, aes(x = AA_Codon, y = pearson_R, fill = AA3)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_pearson), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map_AA3) +  # 用 AA3 配色
  labs(x = "AA:Codon", y = "Pearson's R") +
  theme_classic() +
  scale_y_continuous(breaks = c(seq(0,0.7,0.1)),limits = c(0,0.7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none")

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

# 出图

print(p_pearson)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_allSample_eachcodon_pearson.pdf'), width= 6 , height= 3 , units='in')

print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_allSample_eachcodon_spearman.pdf'), width= 14 , height= 3 , units='in')

# ===================== 核心新增：给 cor_table 加密码子第3位CG/AT列 =====================
# 1. 提取密码子第3位单碱基
cor_table$Codon3_base <- substr(cor_table$Codon, 3, 3)
# 2. 分类：第3位是C/G标为CG，A/T标为AT
cor_table$Codon3_type <- ifelse(cor_table$Codon3_base %in% c("C", "G"), "CG", "AT")

# ===================== 箱线图：AT vs CG · spearman_rho =====================
p_box <- ggplot(cor_table, aes(x = Codon3_type, y = spearman_rho, fill = Codon3_type)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # 散点叠加上去
  scale_y_continuous(breaks = seq(0.3,0.6,0.1),limits = c(0.3,0.63))+
  # Wilcoxon 检验（两组非参数检验）
  ggsignif::geom_signif(comparisons = list(c("AT", "CG")), 
              test = "wilcox.test", 
              map_signif_level = TRUE, 
              textsize = 6, y_position = 0.6) +
  
  scale_fill_manual(values = c("AT" = "#4088B4", "CG" = "#fc8d62")) +  # 清爽配色
  labs(x = "Codon 3rd base type", y = "Spearman's rho\n(vTBI vs. fold change of A-site stalling rate)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

# 出图
print(p_box)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Asite_CG_spearman.pdf'), width= 5 , height= 5 , units='in')



#======================================================================
# 1）先生成相关系数表格（你之前的修复版）
#======================================================================
cor_table <- Asite %>%
  group_by(AA, Codon) %>%
  summarise(
    pearson_R   = cor.test(virus/(virus+Host), (Psite_Peak/Psite_total)/(Ctrl_Pstall), method="pearson")$estimate,
    pearson_P   = cor.test(virus/(virus+Host), (Psite_Peak/Psite_total)/(Ctrl_Pstall), method="pearson")$p.value,
    spearman_rho= cor.test(virus/(virus+Host), (Psite_Peak/Psite_total)/(Ctrl_Pstall), method="spearman")$estimate,
    spearman_P = cor.test(virus/(virus+Host), (Psite_Peak/Psite_total)/(Ctrl_Pstall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA, ":", Codon))

# 显著性 *
cor_table <- cor_table %>%
  mutate(
    p_star_pearson = case_when(
      pearson_P < 0.001 ~ "***",
      pearson_P < 0.01 ~ "**",
      pearson_P < 0.05 ~ "*", TRUE ~ ""
    ),
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

# ===================== 读取 & 构建排序依据 =====================
mid <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/HEK293FT_RSCU.txt")

# 排序：先 AA3 → 再 RSCU 降序（你要求的排序）
 mid_sorted <- mid %>%
  mutate(
    Codon3 = substr(Codon, 3, 3),
    Codon3_group = ifelse(Codon3 %in% c("C", "G"), "CG", "AT")  # 第3位 CG/AT
  ) %>%
  arrange(
    AA3,                  # 第一排序：氨基酸
    Codon3_group,         # 第二排序：CG 在前，AT 在后
    desc(RSCU)            # 第三排序：RSCU 从大到小
  )

# 构建 mid1 = AA3:Codon
mid1 <- paste0(mid_sorted$AA3, ":", mid_sorted$Codon)

# ===================== 给 cor_table 设置横坐标因子顺序（关键！） =====================
# 让 ggplot 横坐标严格按 mid1 排序
cor_table$AA_Codon <- factor(cor_table$AA_Codon, levels = mid1)

# ===================== 绘图：Pearson（按 AA3 着色 + 固定顺序） =====================
p_pearson <- ggplot(cor_table, aes(x = AA_Codon, y = pearson_R, fill = AA3)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_pearson), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map_AA3) +  # 用 AA3 配色
  labs(x = "AA:Codon", y = "Pearson's R") +
  theme_classic() +
  scale_y_continuous(breaks = c(seq(0,0.7,0.1)),limits = c(0,0.7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none")

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


print(p_pearson)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_allSample_eachcodon_pearson.pdf'), width= 6 , height= 3 , units='in')

print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_allSample_eachcodon_spearman.pdf'), width= 14 , height= 3 , units='in')

# ===================== 核心新增：给 cor_table 加密码子第3位CG/AT列 =====================
# 1. 提取密码子第3位单碱基
cor_table$Codon3_base <- substr(cor_table$Codon, 3, 3)
# 2. 分类：第3位是C/G标为CG，A/T标为AT
cor_table$Codon3_type <- ifelse(cor_table$Codon3_base %in% c("C", "G"), "CG", "AT")

# ===================== 箱线图：AT vs CG · spearman_rho =====================
p_box <- ggplot(cor_table, aes(x = Codon3_type, y = spearman_rho, fill = Codon3_type)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # 散点叠加上去
  scale_y_continuous(breaks = seq(0.3,0.6,0.1),limits = c(0.3,0.62))+
  # Wilcoxon 检验（两组非参数检验）
  ggsignif::geom_signif(comparisons = list(c("AT", "CG")), 
                        test = "wilcox.test", 
                        map_signif_level = TRUE, 
                        textsize = 6, y_position = 0.6) +
 scale_fill_manual(values = c("AT" = "#4088B4", "CG" = "#fc8d62")) +  # 清爽配色
  labs(x = "Codon 3rd base type", y = "Spearman's rho\n(vTBI vs. fold change of AP-site stalling rate)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

# 出图
print(p_box)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Psite_CG_spearman.pdf'), width= 3 , height= 5 , units='in')


# ==========================================

# AA

# ==========================================
AA=Asite2 %>% 
  group_by(Sample_standard,Study,AA,virus_S) %>% 
  dplyr::summarise(AA_Peak=sum(Asite_Peak),AA_total=sum(Asite_total),AA_virus=sum(virus),AA_Host=sum(Host))

head(AA)
ctrl <- filter(Asite2, grepl("Ctrl|mock|inf|00hpi", Sample_standard, ignore.case = TRUE)) %>% 
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
  ggpubr::stat_cor(method = "pearson",,label.y = 20,size=5, cor.coef.name = "R") +  
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
    pearson_R   = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Astall), method="pearson")$estimate,
    pearson_P   = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Astall), method="pearson")$p.value,
    spearman_rho= cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Astall), method="spearman")$estimate,
    spearman_P = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Astall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA))

# 显著性 *
cor_table <- cor_table %>%
  mutate(
    p_star_pearson = case_when(
      pearson_P < 0.001 ~ "***",
      pearson_P < 0.01 ~ "**",
      pearson_P < 0.05 ~ "*", TRUE ~ ""
    ),
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
# 图1：Pearson R（使用你指定的氨基酸配色）
#======================================================================
p_pearson <- ggplot(cor_table, aes(x = AA3, y = pearson_R, fill = AA)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_pearson), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map) +  # 已替换为你的配色
  labs(x = "AA", y = "Pearson's R") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "none")

#======================================================================
# 图2：Spearman Rho（使用你指定的氨基酸配色）
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

#======================================================================
# 出图
#======================================================================
print(p_pearson)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_relative_eachAA_pearson.pdf'), width= 6 , height= 3 , units='in')

print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AsiteStall_relative_eachAA_spearman.pdf'), width= 6 , height= 3 , units='in')


# ==========================================

# AA

# ==========================================
AA=Asite2 %>% 
  group_by(Sample_standard,Study,AA,virus_S) %>% 
  dplyr::summarise(AA_Peak=sum(Psite_Peak),AA_total=sum(Psite_total),AA_virus=sum(virus),AA_Host=sum(Host))

head(AA)
ctrl <- filter(Asite2, grepl("Ctrl|mock|inf|00hpi", Sample_standard, ignore.case = TRUE)) %>% 
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
    pearson_R   = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Pstall), method="pearson")$estimate,
    pearson_P   = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Pstall), method="pearson")$p.value,
    spearman_rho= cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Pstall), method="spearman")$estimate,
    spearman_P = cor.test(AA_virus/(AA_virus+AA_Host), (AA_Peak/AA_total)/(Ctrl_Pstall), method="spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(AA_Codon = paste0(AA))

# 显著性 *
cor_table <- cor_table %>%
  mutate(
    p_star_pearson = case_when(
      pearson_P < 0.001 ~ "***",
      pearson_P < 0.01 ~ "**",
      pearson_P < 0.05 ~ "*", TRUE ~ ""
    ),
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
# 图1：Pearson R（使用你指定的氨基酸配色）
#======================================================================
p_pearson <- ggplot(cor_table, aes(x = AA3, y = pearson_R, fill = AA)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = p_star_pearson), vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_map) +  # 已替换为你的配色
  labs(x = "AA", y = "Pearson's R") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "none")

#======================================================================
# 图2：Spearman Rho（使用你指定的氨基酸配色）
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

#======================================================================
# 出图
#======================================================================
print(p_pearson)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_relative_eachAA_pearson.pdf'), width= 6 , height= 3 , units='in')

print(p_spearman)
ggsave(paste0('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/PsiteStall_relative_eachAA_spearman.pdf'), width= 6 , height= 3 , units='in')

