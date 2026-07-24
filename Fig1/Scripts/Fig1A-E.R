#####################################################################

#library packages, set path and plot theme 
#~/miniconda3/envs/DEGreport/bin/R

#####################################################################
{
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  se <- function(x) sd(x)/sqrt(length(x))
  library(car)
  
  my_theme2 <- theme_classic(base_line_size = 1.5,base_size = 20)+
    theme(axis.line.x=element_line(size=1.5,color="black"),
          axis.line.y=element_line(size=1.5,color="black"),
          axis.ticks=element_line(size=0.5,color="black"),
          axis.ticks.length=unit(0.05,"inches"),
          axis.title.x = element_text(size=15,color="black"),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(angle = 90,hjust = 1,size=15,color="black"),
          axis.text.y =  element_text(size=15,color="black"),
          strip.text.x = element_text(size=15,face = "bold"),
          #strip.background = element_rect(color = "black",size=0),
          legend.position = 1,plot.title = element_text(hjust = 0.5,size=15),
          legend.text = element_text(size=10),legend.title = element_text(size=10),
          axis.line.x.top = element_line(size=0,color="black"),
          axis.line.y.right = element_line(size=0,color="black"))
  setwd("C:/Users/Ke-jia Shan/Desktop/Virus_CUB/Fig1/")
  
  
  my_theme <- theme_classic(base_line_size = 1.5,base_size = 20)+
    theme(axis.line.x=element_line(size=1.5,color="black"),
          axis.line.y=element_line(size=1.5,color="black"),
          axis.ticks=element_line(size=0.5,color="black"),
          axis.ticks.length=unit(0.05,"inches"),
          axis.title.x = element_text(size=15,color="black"),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(angle = 45,hjust = 1,size=15,color="black"),
          axis.text.y =  element_text(size=15,color="black"),
          strip.text.x = element_text(size=15,face = "bold"),
          #strip.background = element_rect(color = "black",size=0),
          legend.position = 1,plot.title = element_text(hjust = 0.5,size=15),
          legend.text = element_text(size=10),legend.title = element_text(size=10),
          axis.line.x.top = element_line(size=0,color="black"),
          axis.line.y.right = element_line(size=0,color="black"))
  
  class_cols <- c(
    "other mammals"        = "#8A8A8A",  # neutral gray
    "human"                = "#C94C47",  # soft red
    "bat"                  = "#8E55A6",  # purple
    "Liliopsida_Plant"     = "#006837FF",
    "eudicotyledons_Plant" = "#7CAE00FF",
    "other plants"         = "#b8ebb0",
    "fish"                 = "#E3B36B",  # peach/orange
    "reptile"              = "#9A6A45",  # brown
    "bird"                 = "#D88AAF",  # rose pink
    "DNA viruses"          = "#2B5C9A",  # deep blue
    "RNA viruses"          = "#57A9D9"   # cyan blue
  )
}


##########################################################

#RSCUij umap
#Figure 1A and 1B
#Note that the results of each UMAP run may vary. Despite differences in coordinates, the clustering classification remains consistent every time.

################################################
{
  RSCUij=fread("./Results/RSCU.csv")
  head(RSCUij)
  table(RSCUij$Class)
  matrix=RSCUij[,-c(1:2)]
  head(matrix)
  head(matrix[,1:61])
  
  
  #umap
  library(umap)
  #
  umap_result <- umap(matrix, n_neighbors = 20, min_dist = 0.1, metric = "euclidean", scale = FALSE)
  aa_umap <- data.frame(Class = RSCUij$Class, Species = RSCUij$Species, UMAP1 = umap_result$layout[,1], UMAP2 = umap_result$layout[,2])
  table(aa_umap$Class)
  aa_umap[grepl("DNA",aa_umap$Class),]$Class="DNA virus"
  aa_umap[grepl("RNA",aa_umap$Class),]$Class="RNA virus"
  
 # write.table(aa_umap2,"./Results/umap_plot2.txt",row.names = F,quote = F,sep = "\t")
  
  aa_umap=fread("./Results/umap_plot2.txt")
  #aa_umap=data.table(aa_umap2)
  #aa_umap[aa_umap$Class %in% c("Amborellales_Plant","Embryophyta_Plant"),]$Class="other plants"
  #aa_umap[aa_umap$Class %in% c("artiodactyl","perissodactyla","Platypus","Rabbit","MarineMammals","hedgehog","mammal","rodent","primate"),]$Class="other mammals"
  table(aa_umap$Class)
  head(aa_umap)
  #write.table(aa_umap,"umap_plot20260208.txt",row.names = F,quote = F,sep = "\t")
  
  #aa_umap=fread("umap_plot2.txt")
  aa_umap$Class=factor(aa_umap$Class,levels = c("Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                                "fish", "reptile", "bird", 
                                                "other mammals", "human","bat",
                                                "DNA virus", "RNA virus"))
  table(aa_umap$Class)
  ggplot()  +
    labs(x = "", y = "") +
    #DNA virus
    geom_point(data = filter(aa_umap,grepl("DNA virus",Class)), mapping = aes(x = UMAP1, y = UMAP2),col ="blue",
               size = 1,
               #alpha = 0.5,
               stroke = 0) +
    #RNA virus
    geom_point(data = filter(aa_umap,grepl("RNA virus",Class)), mapping = aes(x = UMAP1, y = UMAP2),col ="#57C3F3",
               size = 1,
               #alpha = 0.5,
               stroke = 0) +
    #host
    geom_point(data = filter(aa_umap,!grepl("virus",Class),!grepl("bat",Class),!grepl("human",Class)), mapping = aes(x = UMAP1, y = UMAP2, col = Class),
               size = 1,
               #alpha = 0.5,
               stroke = 0) +
    scale_color_manual(values =c("#00FF00","#006400","#b8ebb0",
                                 "#B09C85FF","black","#FFC0CB",
                                 "#ADB6B6FF","red","#631879FF",
                                 "blue","#57C3F3")) +
    scale_color_manual(values =class_cols) +
    #bat
    geom_point(data = filter(aa_umap,grepl("bat",Class)), mapping = aes(x = UMAP1, y = UMAP2),col ="#631879FF",
               size = 1,
               #alpha = 0.5,
               stroke = 0) +
    theme_classic(base_line_size = 1, base_size = 20)+
    #human
    geom_point(data = filter(aa_umap,grepl("human",Class)), mapping = aes(x = UMAP1, y = UMAP2),col ="red",
               size = 1,
               #alpha = 0.5,
               stroke = 0)+
    theme(
      legend.position = "none",                     # 去掉图例
      axis.line.x = element_line(size = 1, color = "black"),
      axis.line.y = element_line(size = 1, color = "black")#,
      #axis.text = element_blank()   # 一次性移除 X 和 Y 轴所有刻度标签
    )

}    
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/1A_UMAP.pdf', width= 4 , height= 4 , units='in')


aa_umap=fread("./Results/umap_plot2.txt")
#The file EnsemblHost.virus.CUB.txt is uploaded as a gz compressed archive on GitHub due to file size upload limitations.
Virus_CUB<-fread("./Results/EnsemblHost.virus.CUB.txt",header = F,sep = "\t")
colnames(Virus_CUB)=c("Host_class","Host","Organism_name","species","Codon","Num")
Virus_CUB$class="Virus"
Virus_CUB=filter(Virus_CUB,Organism_name %in% aa_umap$Species)
head(Virus_CUB)
aa_umap[grepl("virus",Class) & !Species %in% Virus_CUB$Organism_name]
table(aa_umap$Class)
unique(Virus_CUB$Organism_name)
head(aa_umap)
#data
Host_Virus=merge(filter(aa_umap,grepl("virus",Class)),Virus_CUB[,c(1,3)],by.x="Species",by.y="Organism_name")
colnames(Host_Virus)[2]="Virus_Class"

Host_Virus[Host_Virus$Host_class %in% c("Amborellales_Plant","Embryophyta_Plant","Chlorophyta_Plant"),]$Host_class="other plants"
Host_Virus[Host_Virus$Host_class %in% c("artiodactyl","perissodactyla","Platypus","Rabbit","MarineMammals","hedgehog","mammal","rodent","primate"),]$Host_class="other mammals"
head(Host_Virus)
unique(Host_Virus$Host_class)

unique(Host_Virus$Species)
unique(Host_Virus$Host_class)
Host_Virus$Class=factor(Host_Virus$Host_class,levels = c("Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                                         "fish", "reptile", "bird",
                                                         "other mammals", "human","bat"))
ggplot()  +
  labs(x = "", y = "") +
  geom_point(data = filter(Host_Virus,UMAP1<=5), 
             mapping = aes(x = UMAP1, y = UMAP2,col =Class),
             size=2,
             alpha = 0.02,
             stroke = 0) +
  # scale_color_manual(values =c("#00FF00","#55A84F","#b8ebb0",
  #                              "#B09C85FF","black","#FFC0CB",
  #                              "#ADB6B6FF","red","#631879FF",
  #                              "blue","#57C3F3")) +
  scale_color_manual(values =class_cols) +
  #scale_shape_manual(values =c(3,4)) +
  ylim(-4,6)+
  xlim(-7.5,5)+
  theme_classic(base_line_size = 1.5, base_size = 20)+
  theme(
    legend.position = "none",                     # 去掉图例
    axis.line.x = element_line(size = 1.5, color = "black"),
    axis.line.y = element_line(size = 1.5, color = "black")#,
    #axis.text = element_blank()   # 一次性移除 X 和 Y 轴所有刻度标签
  )
#ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/1A_UMAP_Host.pdf', width= 5 , height= 5 , units='in')
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/1A_UMAP_Host.jpg', width= 7 , height= 7 , units='in')


#################################################################

#RSCUij heatmap
#S1

################################################################

CUB=fread("./Results/RSCU.csv")
head(CUB)
aa_umap=fread("./Results/umap_plot2.txt")
head(aa_umap)
aa2=merge(aa_umap[,1:2],CUB[,-c(1)],by="Species")

aa2$Class=factor(aa2$Class,levels = c("Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                      "fish", "reptile", "bird", 
                                      "other mammals", "human","bat",
                                      "RNA virus", "DNA virus"))



#gap between different area
gaps = table(aa2$Class)
gaps_split = cumsum(gaps)


row_ann = data.frame(class= as.character(aa2$Class))
row.names(row_ann) =aa2$Species
head(row_ann)
row.names(aa2) = aa2$Species
LD_matx=aa2
head(LD_matx[1:4,1:4])
ncol(LD_matx)


library(pheatmap)
#dd_order = hclust(dist(t(abs(dd))), method = 'ward.D2')
row_ann=row_ann %>%
  mutate(class= factor(class, levels = c("Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                         "fish", "reptile", "bird", 
                                         "other mammals", "human","bat",
                                         "RNA virus", "DNA virus"))) %>% arrange(class)

head(row_ann)
#row_ann=row_ann[,-c("category")]
ann_colors=list(class=c(Liliopsida_Plant='#006400',eudicotyledons_Plant='#00FF00',`other plants`='#53A85F',
                        fish='#E18727FF',reptile='black',bird='#FFC0CB',`other mammals`='#ADB6B6FF',
                        human='red', bat='#631879FF',`RNA virus`='#57C3F3',`DNA virus`='#0000FF')) 

dd = as.matrix(LD_matx[,3:(dim(LD_matx)[2])])
rownames(dd)=LD_matx$Species
pheatmap(
  (dd[rownames(row_ann),]),#dd_order$order
  cluster_rows = F, 
  cluster_cols = F,
  annotation_colors=ann_colors,
  #treeheight_row = 0,
  #colnames_row = ,
  #annotation_col = col_ann[,-1],
  annotation_row = row_ann,
  show_colnames = T, 
  show_rownames = F, 
  angle_col = 90, 
  fontsize_col = 10,
  gaps_row = gaps_split[1:(length(gaps)-1)], 
  border_color = 'white', 
  #filename = 'C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/heatmap.pdf', 
  width = 10, height = 15, legend = T#, color = colorRampPalette(c("black","blue", "purple","white", "orange","red"))(100)
)

head(dd)


#################################################################

#Figure 1E  CG fraction

################################################################

test=fread("./Results/CUB.txt")
test$First_base <- substr(test$Codon, 1, 1) #substr("ATG", 1, 1)
test$Second_base <- substr(test$Codon, 2, 2)#substr("ATG", 2, 2)
test$Third_base <- substr(test$Codon, 3, 3) #substr("ATG", 3, 3)
head(test)

length(unique(test$Organism_name))

Third_base = test %>% group_by(Organism_name,Third_base)  %>% dplyr::summarise(BaseNum=sum(Num))
length(unique(Third_base$Organism_name))
Third_base <- Third_base %>%
  mutate(ID = str_replace(Organism_name, "^.*:(GCF_[^:]+).*$", "\\1"))
head(Third_base)

aa_umap=fread("./Results/umap_plot2.txt")
aa_umap <- aa_umap %>%
  mutate(ID = str_replace(Species, "^.*:(GCF_[^:]+).*$", "\\1"))

nrow(aa_umap)
head(aa_umap)
#Third_base[!Third_base$ID %in% aa_umap$ID,]
aa_umap[!aa_umap$ID %in% Third_base$ID,]

aa2=merge(aa_umap[,c(1,5,2)],Third_base,by.x="ID",by.y="ID")
head(aa2)

Third_base2 = aa2 %>% group_by(Species,Class)  %>% dplyr::summarise(Total=sum(BaseNum))
nrow(Third_base2)
aa2=merge(aa2,Third_base2,by=c("Species","Class"))
aa2$CG="CG"
aa2[aa2$Third_base %in% c("T","A")]$CG="AT"
head(aa2)
Third_base2 = aa2 %>% group_by(Species,Class,CG)  %>% dplyr::summarise(CG3_Frac=sum(BaseNum/Total))
head(Third_base2)
length(unique(Third_base2$Species))

Third_base2$Class=gsub(" ","",Third_base2$Class)
Third_base2[Third_base2$Class=="human",]$Class="othermammals"
Third_base2[Third_base2$Class=="bat",]$Class="othermammals"

Third_base2$Class=factor(Third_base2$Class,levels = c("othermammals", "bird",
                                                      "reptile", "fish", 
                                                      "Liliopsida_Plant", "eudicotyledons_Plant","otherplants",
                                                      "RNAvirus", "DNAvirus"))


data=filter(Third_base2,CG=="CG")
# 获取所有不同的 Class 水平
classes <- unique(data$Class)

# 初始化一个空的数据框来存储结果
results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# 进行两两比较
for (i in 1:(length(classes) - 1)) {
  for (j in (i + 1):length(classes)) {
    group1 <- data$CG3_Frac[data$Class == classes[i]]
    group2 <- data$CG3_Frac[data$Class == classes[j]]
    test_result <- wilcox.test(group1, group2)
    results <- rbind(results, data.frame(
      Group1 = as.character(classes[i]),
      Group2 = as.character(classes[j]),
      p_value = test_result$p.value
    ))
  }
}
results
# 创建一个对称的 p 值矩阵
p_matrix <- matrix(1, nrow = length(classes), ncol = length(classes))
rownames(p_matrix) <- levels(classes)
colnames(p_matrix) <- levels(classes)
for (k in 1:nrow(results)) {
  p_matrix[results$Group1[k], results$Group2[k]] <- results$p_value[k]
  p_matrix[results$Group2[k], results$Group1[k]] <- results$p_value[k]
}

# 计算每个组的中位数
medians <- tapply(data$CG3_Frac, data$Class, median)

# 按中位数大小对组别进行排序（从大到小）
sorted_classes <- names(sort(medians, decreasing = TRUE))  # 这里添加了decreasing = TRUE

# 重新排序p矩阵的行和列，使其与中位数排序一致
p_matrix_sorted <- p_matrix[sorted_classes, sorted_classes]

# 使用排序后的矩阵生成字母标记，确保字母顺序反映中位数大小（从大到小）
letters <- multcompView::multcompLetters(p_matrix_sorted, compare="<", threshold = 0.05)$Letters

# 保持原始类别名称与字母的对应关系
letters_df <- data.frame(
  Class = names(letters), 
  Letter = as.character(letters),
  stringsAsFactors = FALSE
)

# 合并字母标记到原始数据
data <- merge(data, letters_df, by = "Class")
nrow(data)
# 绘制箱线图并标记字母
ggplot(data, aes(x = Class, y = CG3_Frac,col=as.factor(Class))) +
  #geom_violin() +
  geom_boxplot() +
  #geom_jitter(mapping=aes(col=as.factor(Class)),alpha=0.2)+
  geom_text(data = letters_df, aes(x = Class, y = max(data$CG3_Frac), label = Letter),color="black", vjust = -1) +
  theme_classic()+
  scale_y_continuous(breaks = c(seq(0.1,1,0.1)),limits = c(0.1,1)) +
  scale_color_manual(values =c("#ADB6B6FF","#FFC0CB",
                               "black","#B09C85FF",
                               #"#631879FF",
                               "#006400","#00FF00","#b8ebb0",
                               "#57C3F3","blue"))+
  theme(legend.position = "none")
#my_theme

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/CG3.pdf', width= 5.135 , height= 3.325 , units='in')


#################################################################

#host VS virus CAI
#figure 1C

################################################################

{
  setwd("C:/Users/Ke-jia Shan/Desktop/Virus_CUB/Fig1/Results/")
  aa_umap=fread("umap_plot2.txt")[,1:2]
  CAI_mid=fread("Host_CAI_median.txt")
  Virus_CAI=fread("Virus_CAI.txt")
  #Virus_CAI=fread("Virus_CAI.cov.txt")
  
  head(Virus_CAI)
  length(unique(filter(Virus_CAI,Organism_name %in% aa_umap$Species)$Organism_name))
  
  Host_Virus_CAI=merge(Virus_CAI,CAI_mid,by=c("Host"))
  
  #Host class 
  aa_umap[aa_umap$Class %in% c("Amborellales_Plant","Embryophyta_Plant"),]$Class="other plants"
  aa_umap[aa_umap$Class %in% c("artiodactyl","perissodactyla","Platypus","Rabbit","MarineMammals","hedgehog","mammal","rodent","primate"),]$Class="other mammals"
  table(aa_umap$Class)
  aa_umap[grepl("DNA",Class)]$Class="DNA virus"
  aa_umap[grepl("RNA",Class)]$Class="RNA virus"
  head(aa_umap)
  
  head(Host_Virus_CAI)
  Host_Virus_CAI=merge(Host_Virus_CAI,aa_umap,by.x="Host",by.y="Species")
  colnames(Host_Virus_CAI)[8]="Host_Class"
  Host_Virus_CAI=merge(Host_Virus_CAI,aa_umap,by.x="Organism_name",by.y="Species")
  colnames(Host_Virus_CAI)[9]="Virus_Class"
  
  head(Host_Virus_CAI)
  
  #unique(sub(pattern = "_.*", replacement = "", x = Host_Virus_CAI$Host))
  Host_Virus_CAI$Host_Class=factor(Host_Virus_CAI$Host_Class,levels = c( "human","bat",
                                                                         "Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                                                         "fish", "reptile", "bird", 
                                                                         "other mammals",
                                                                         "RNA virus", "DNA virus"))
  
  
}
# # 通用格式：去掉str中第一个"_"及之后的所有字符
library(dplyr)
library(stringr)
library(purrr)
Host_Virus_CAI$Host_genus <- sub(pattern = "_.*", replacement = "", x = Host_Virus_CAI$Host)
head(Host_Virus_CAI)
Host_Virus_CAI <- Host_Virus_CAI %>%
  mutate(
    Species = Host %>%
      str_split(pattern = "\\s+|_|\\.") %>%  # 按空格、下划线、点分割
      map_chr(~ str_c(.x[1], .x[2], sep = "_")) %>%  # 前两个元素用_拼接
      str_to_lower()  # 全部转为小写
  )

mid <- Host_Virus_CAI %>%
  group_by(Host_genus,Organism_name) %>%
  slice_max(CAI, with_ties = FALSE) %>%  # 只留1行最大CAI
  select(Organism_name, Host_Class, Virus_Class, Host_genus, CAI, Host, median,CodonSum) %>%
  ungroup()
write.table(mid,"max_CAI_Virus_host2.txt",row.names = F,quote = F,sep = "\t")


mid=fread("max_CAI_Virus_host2.txt")
mid$Host_Class <- factor(mid$Host_Class,
                         levels = c("other mammals","human","bat",
                                    "Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                    "fish", "reptile", "bird"))
ggplot(data = mid, 
       mapping = aes(x = median, y = CAI, color = Host_Class)) +
  
  # 先画底层点
  geom_point(data = filter(mid, Host_Class %in% c("other mammals","Liliopsida_Plant", "eudicotyledons_Plant","other plants", "reptile")),
             size = 2,
             alpha = 0.7,
             stroke = 0) +
  geom_point(data = filter(mid, Host_Class %in% c("fish", "bird")),
             size = 2,
             alpha = 0.7,
             stroke = 0) +  #  bat 放后面置顶
  
  # 全局拟合线：只画 1 条！关键：不继承 color 分组
  geom_smooth(
    aes(x = median, y = CAI),  
    method = "glm", 
    se = FALSE, 
    color = "blue",  # 统一黑色
    linewidth = 1,
    inherit.aes = FALSE  # 核心：不跟着颜色分组
  ) +
  
  # 置顶 bat + human
  geom_point(data = filter(mid, Host_Class == "bat"),
             size = 2,
             alpha = 0.7,
             stroke = 0) +
  geom_point(data = filter(mid, Host_Class == "human"),
             size = 2,
             alpha = 0.7,
             stroke = 0) +
  
  # 对角线
  geom_abline(intercept = 0, slope = 1, col = "grey75", linetype = "dashed") +
  
  scale_y_continuous(breaks = seq(0.2,0.8,0.1), limit=c(0.27,0.8)) +
  scale_x_continuous(breaks = seq(0.4,0.8,0.1), limit=c(0.49,0.8)) +
  
  # 颜色完全不变
  # scale_color_manual(values = c(
  #   "other mammals"       = "#ADB6B6FF",
  #   "human"               = "red",
  #   "bat"                 = "#631879FF",
  #   "Liliopsida_Plant"    = "#006400",
  #   "eudicotyledons_Plant"= "#00FF00",
  #   "other plants"        = "#b8ebb0",
  #   "fish"                = "#B09C85FF",
  #   "reptile"             = "black",
  #   "bird"                = "#FFC0CB"
  # )) +
  # 颜色完全不变
  scale_color_manual(values = class_cols)+
  #labs(x = "Host CAI median", y = "Viral CAI") +
  labs(x = "", y = "") +
  theme_classic(base_line_size = 1.5, base_size = 20) +
  theme(
    axis.text = element_blank(),
    #axis.ticks = element_blank(),
    legend.position = "none"
  )

global_model <- glm(CAI ~ median, data = mid)

# 输出结果
summary(global_model)
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Host_Viral_maxCAI.pdf', width= 7 , height= 5 , units='in')

cor.test(mid$CAI,
         mid$median,method = "s")

#median
head(Host_Virus_CAI)
mid <- Host_Virus_CAI %>%
  group_by(Host_genus, Organism_name) %>%
  summarise(
    MedianVirus = median(CAI, na.rm = TRUE),
    median = first(median),
    Host_Class = first(Host_Class),
    Virus_Class = first(Virus_Class),
    Host = first(Host),
    CodonSum = first(CodonSum),
    .groups = "drop"
  )

head(mid)
mid$Host_Class <- factor(mid$Host_Class,
                         levels = c("other mammals","human","bat",
                                    "Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                    "fish", "reptile", "bird"))
ggplot(data = mid, 
       mapping = aes(x = median, y = MedianVirus, color = Host_Class)) +
  
  # 先画底层点
  geom_point(data = filter(mid, Host_Class %in% c("other mammals","Liliopsida_Plant", "eudicotyledons_Plant","other plants", "reptile")),size=0.7) +
  geom_point(data = filter(mid, Host_Class %in% c("fish", "bird")),size=0.7) +  #  bat 放后面置顶
  
  # 全局拟合线：只画 1 条！关键：不继承 color 分组
  geom_smooth(
    aes(x = median, y = MedianVirus),  
    method = "glm", 
    se = FALSE, 
    color = "blue",  # 统一黑色
    linewidth = 1,
    inherit.aes = FALSE  # 核心：不跟着颜色分组
  ) +
  
  # 置顶 bat + human
  geom_point(data = filter(mid, Host_Class == "bat"),size=0.7) +
  geom_point(data = filter(mid, Host_Class == "human"),size=0.7) +
  
  # 对角线
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
  
  scale_y_continuous(breaks = seq(0.2,0.8,0.1), limit=c(0.27,0.8)) +
  scale_x_continuous(breaks = seq(0.4,0.8,0.1), limit=c(0.49,0.8)) +
  ggpubr::stat_cor(method = "spearman",col="red",size=5,cor.coef.name = "rho")+
  # 颜色完全不变
  #scale_color_manual(values = class_cols)+
  labs(x = "Host CAI median", y = "Viral CAI (genus median)") +
  #labs(x = "", y = "") +
  theme_classic(base_line_size = 1.5, base_size = 20) +
  theme(
    #axis.text = element_blank(),
    #axis.ticks = element_blank(),
    legend.position = "none"
  )



#################################################################

#figure 1D

################################################################
table(violin_plot[,c(2,4)])
violin_plot=fread("Fig1D_violin_plot.txt")
#install.packages("gghalves")
#library(gghalves)  
ggplot(violin_plot,
       aes(x = paste0(Host_Class,"\n",Host), y = CAI, fill = class)) +
  
  # RNA 左半边
  geom_half_violin(
    data = . %>% filter(class == "host"),
    side = "l", alpha = .7, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "host"),
    side = "l", width = .2, outlier.shape = NA, color = "black"
  ) +
  
  # DNA 右半边
  geom_half_violin(
    data = . %>% filter(class == "virus"),
    side = "r", alpha = .7, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "virus"),
    side = "r", width = .2, outlier.shape = NA, color = "black"
  ) +
  
  scale_fill_manual(values = c("#E41A1C","#377EB8")) +
  #labs(x = "Host class", y = "CAI", fill = "Virus class") +
  theme_classic()
#theme(axis.text.x = element_text(angle = 45,hjust = 1,color="black"))

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Host_Viral_medianCAI.pdf', width= 14 , height= 8 , units='in')


# ggplot(violin_plot,
#        aes(x =  class, y = CAI, fill = class)) +
# 
#   geom_boxplot( ) +
#   scale_fill_manual(values = c("#E41A1C","#377EB8")) +
#   #labs(x = "Host class", y = "CAI", fill = "Virus class") +
#   theme_classic()+
#   facet_grid(~Host_Class)+
#   ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = list(c("host","virus")),label = 'p.signif')





