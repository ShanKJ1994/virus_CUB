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
  setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/")
  
  
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

########################################################

#Class

################################################################


ClassAnimal<- fread("Animal_class.txt",header = F,sep = "\t")
colnames(ClassAnimal)=c("name","long","class","species")
head(ClassAnimal)

ClassPlant<- fread("Plant_class.txt",header = T,sep = "\t")
colnames(ClassPlant)=c("species","class")
head(ClassPlant)

Class=rbind(ClassAnimal[,c(4,3)],ClassPlant)
head(Class)  


########################################################

#different species codon usage

################################################################


Class<- fread("Animal_class.txt",header = F,sep = "\t")
colnames(Class)=c("name","long","class","species")
head(Class)

CUB<- fread("Animal.ribosomalprotein.CUB.txt",header = F,sep = "\t")
head(CUB)
CUB<-CUB[(!CUB$V1 %in% c("TAA","TAG","TGA")),]

# 定义函数处理每个元素
process_string <- function(str) {
  split_string <- unlist(strsplit(str, "[_\\.]"))
  first_two_elements <- head(split_string, 2)
  result_string <- paste(first_two_elements, collapse = "_")
  return(result_string)
}

# 对向量中的每个元素应用函数
result_vector <- lapply(CUB$V3, process_string)


CUB$species <- unlist(result_vector)
nrow(CUB)
CUB=merge(CUB,unique(Class[,3:4]),by=c("species"))
CUB<-CUB[(!CUB$V2 %in% c("TAA","TAG","TGA")),]
nrow(CUB)
head(CUB)

Class<- fread("Plant_class.txt",header = T,sep = "\t")
colnames(Class)=c("species","class")
head(Class)

CUB2<- fread("Plant.ribosomalprotein.CUB.txt",header = F,sep = "\t")
CUB2<-CUB2[(!CUB2$V1 %in% c("TAA","TAG","TGA")),]
CUB2$species <- unlist(lapply(CUB2$V3, process_string))
nrow(CUB2)
head(CUB2)

CUB2=merge(CUB2,unique(Class),by=c("species"))
CUB2<-CUB2[(!CUB2$V2 %in% c("TAA","TAG","TGA")),]
nrow(CUB2)
head(CUB2)
Host_CUB=rbind(CUB,CUB2)
colnames(Host_CUB)[2:4]=c("Codon","Num","Organism_name")
head(Host_CUB)

#Virus_CUB<-fread("EnsemblHost.virus.CUB.cov.1a.txt",header = F,sep = "\t")
Virus_CUB<-fread("EnsemblHost.virus.CUB.txt",header = F,sep = "\t")
colnames(Virus_CUB)=c("Host_class","Host","Organism_name","species","Codon","Num")
filter(Virus_CUB,grepl("GCF_009858895.2",Organism_name))[(Codon %in% c("TAA","TAG","TGA"))]
Virus_CUB$class="Virus"
Virus_CUB<-Virus_CUB[(!Virus_CUB$Codon %in% c("TAA","TAG","TGA")),]
head(Virus_CUB)

count=unique(Virus_CUB[,c("Organism_name","Codon","Num","class")]) %>%  
  group_by(Organism_name) %>% dplyr::summarise(count=sum(Num)) %>%  
  filter(count>=100)



#################################################################

#Wij

################################################################

CUB_df<-rbind(unique(filter(Virus_CUB,Organism_name %in% count$Organism_name )[,c("Organism_name","Codon","Num","class")]),Host_CUB[,c("Organism_name","Codon","Num","class")])

CUB<-reshape2::dcast(CUB_df,class+Organism_name~Codon,value.var = "Num")
CUB=na.omit(CUB)
colnames(CUB)[2]<-"Species"
CUB$Sum<-rowSums(CUB[,-c(1:2)])
CUB=filter(CUB,Sum>=100)

test<-filter(CUB_df,!Codon %in% c("TGA","TAG","TAA"),Organism_name %in% CUB$Species)
codon_table<-fread("codon_table.txt",header=F)
colnames(codon_table)=c("Codon","AA3","AA","AAname")
test=merge(test,codon_table[,c(1,3)],by=c("Codon"))
nrow(test)
mid<-test %>% group_by(AA,Organism_name) %>% dplyr::summarise(Mean=mean(Num))
RSCUij<-merge(mid,test,by=c("AA","Organism_name"))
RSCUij$RSCUij=RSCUij$Num/RSCUij$Mean
nrow(RSCUij)
head(RSCUij)
mid<-RSCUij %>% group_by(AA,Organism_name) %>% dplyr::summarise(Max=max(RSCUij))
RSCUij<-merge(mid,RSCUij,by=c("AA","Organism_name"))
RSCUij$Wij=log10(RSCUij$RSCUij/RSCUij$Max)
nrow(RSCUij)
write.table(RSCUij[,-c(3,4)],"WijTable.txt",row.names = F,quote = F,sep = "\t")
#write.table(RSCUij[,-c(3,4)],"WijTable2.cov.txt",row.names = F,quote = F,sep = "\t")

head(RSCUij)
175192/61
length(unique(RSCUij$Organism_name))


#################################################################

#viral CAI

################################################################
#Wij<-read.table("WijTable2.cov.txt",header = T,sep="\t")#filter stop codons, log10
Wij<-read.table("WijTable.txt",header = T,sep="\t")#filter stop codons, log10
Wij$Host= sub("_.*", "", Wij$Organism_name)
Wij=filter(Wij,class!="Virus")
head(Wij)
#Virus_CUB<-fread("EnsemblHost.virus.CUB.cov.1a.txt",header = F,sep = "\t")
Virus_CUB<-fread("EnsemblHost.virus.CUB.txt",header = F,sep = "\t")
colnames(Virus_CUB)=c("Host_class","Host","Organism_name","species","Codon","Num")
Virus_CUB$class="Virus"
head(Virus_CUB)
CUB_virus<-filter(Virus_CUB,!Codon %in% c("TGA","TAG","TAA"),Organism_name %in% CUB$Species)
head(CUB_virus)

CAI<-data.frame()
head(CUB)
#filter(Host_Virus_CAI,grepl("GCF_000859885",Organism_name))

for (i in unique(CUB_virus$Host)) {
  CUB2=filter(Wij,Host==i)
  CUB_virus3=filter(CUB_virus,Host==i)
  for (j in unique(CUB2$Organism_name)) {
    CUB3<-filter(CUB2,Organism_name==j)
    for (virus in unique(CUB_virus3$Organism_name)) {
      CUB_virus2=filter(CUB_virus3,Organism_name==virus)
      codon_table<-merge(CUB_virus2,CUB3[,c("Codon","Wij")],by=c("Codon"))
      codon_table$CAI1<-(codon_table$Wij)*(codon_table$Num)
      CAI_mid<-codon_table %>% group_by(Organism_name) %>%  dplyr::summarise(CAI=10^(sum(CAI1)/sum(Num)),Host=j,CodonSum=sum(Num))
      CAI<-rbind(CAI_mid,CAI)
      print(paste(i,j,virus,sep="  "))
    }
  }
}
#write.table(CAI,"Virus_CAI.cov.txt",row.names = F,quote = F,sep = "\t")
write.table(CAI,"Virus_CAI.txt",row.names = F,quote = F,sep = "\t")
#################################################################

#host CAI

################################################################

Wij<-read.table("WijTable.txt",header = T,sep="\t")[,c(2,3,7)]#filter stop codons, log10
head(Wij)
CAI<-data.frame()
for (i in unique(Wij$Organism_name)) {
  CUB2=filter(Wij,Organism_name==i)
  print(i)
  print(paste0("LongestCDS_CUB/",i,".CUB.txt"))
  species=fread(paste0("LongestCDS_CUB/",i,".CUB.txt"),header=F,sep="\t")
  colnames(species)=c("Codon","Num","ID")
  codon_table<-merge(CUB2,species,by=c("Codon"))
  codon_table$CAI1<-(codon_table$Wij)*(codon_table$Num)
  CAI_mid<-codon_table %>% group_by(ID) %>%  dplyr::summarise(CAI=10^(sum(CAI1)/sum(Num)),Host=i,CodonSum=sum(Num))
  CAI<-rbind(CAI_mid,CAI)
  
}
write.table(CAI,"Host_CAI.txt",row.names = F,quote = F,sep = "\t")

CAI_mid=CAI %>% group_by(Host)  %>%  dplyr::summarise(median=median(CAI),low=quantile(CAI,0.025),high=quantile(CAI,0.975))
write.table(CAI_mid,"Host_CAI_median.txt",row.names = F,quote = F,sep = "\t")


##########################################################

#RSCUij umap
#Figure 1A and 1B

################################################
{
  RSCUij=fread("WijTable.txt")
  #RSCUij=fread("WijTable2.cov.txt")
  
  length(unique(RSCUij$Organism_name))
  CUB<-reshape2::dcast(RSCUij,class+Organism_name~Codon,value.var = "RSCUij")
  nrow(CUB)
  CUB=na.omit(CUB)
  colnames(CUB)[2]<-"Species"
  head(CUB)
  unique(CUB$class)
  CUB <- CUB %>%
    mutate(ID = str_replace(Species, "^.*:(GCF_[^:]+).*$", "\\1"))
  
 
  matrix=CUB[,-c(1:2)]
  head(matrix)
  head(matrix[,1:61])
  
  
  #umap
  library(umap)
  umap_result <- umap(matrix, n_neighbors = 20, min_dist = 0.1, metric = "euclidean", scale = FALSE)
  aa_umap <- data.frame(Class = CUB$class, Species = CUB$Species, UMAP1 = umap_result$layout[,1], UMAP2 = umap_result$layout[,2])
  
  
  # 将sDNA列按冒号分割
  split_data <- strsplit(as.character(aa_umap[aa_umap$Class=="Virus",]$Species), ":")
  # 提取第2列数据
  aa_umap[aa_umap$Class=="Virus",]$Class <-  unlist(sapply(split_data, function(x) x[3]))
  
  unknown=fread("UnknownTypeVirus.txt",header=F)
  head(unknown)
  merged_data <- merge(aa_umap, unknown, by.x = "Species", by.y = "V1")
  merged_data <- merged_data[,c(5,1,3,4)]
  colnames(merged_data)[1]="Class"
  
  aa_umap2=rbind(aa_umap[aa_umap$Class!="unknown",],merged_data)
  aa_umap2$Class <- gsub("\\(\\+\\)|\\(-\\)|\\(\\+/-\\)", "", aa_umap2$Class)
  unique(aa_umap2$Class)
  table(aa_umap2$Class)
  #write.table(aa_umap2,"umap_plot2.txt",row.names = F,quote = F,sep = "\t")
  
  aa_umap=fread("umap_plot2.txt")
  #aa_umap=data.table(aa_umap2)
  aa_umap[aa_umap$Class %in% c("Amborellales_Plant","Embryophyta_Plant"),]$Class="other plants"
  aa_umap[aa_umap$Class %in% c("artiodactyl","perissodactyla","Platypus","Rabbit","MarineMammals","hedgehog","mammal","rodent","primate"),]$Class="other mammals"
  table(aa_umap$Class)
  aa_umap[grepl("DNA",Class)]$Class="DNA virus"
  aa_umap[grepl("RNA",Class)]$Class="RNA virus"
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
               size = 0.5,
               #alpha = 0.5,
               stroke = 0) +
    #RNA virus
    geom_point(data = filter(aa_umap,grepl("RNA virus",Class)), mapping = aes(x = UMAP1, y = UMAP2),col ="#57C3F3",
               size = 0.5,
               #alpha = 0.5,
               stroke = 0) +
    #host
    geom_point(data = filter(aa_umap,!grepl("virus",Class),!grepl("bat",Class),!grepl("human",Class)), mapping = aes(x = UMAP1, y = UMAP2, col = Class),
               size = 0.5,
               #alpha = 0.5,
               stroke = 0) +
    # scale_color_manual(values =c("#00FF00","#006400","#b8ebb0",
    #                              "#B09C85FF","black","#FFC0CB",
    #                              "#ADB6B6FF","red","#631879FF",
    #                              "blue","#57C3F3")) +
    scale_color_manual(values =class_cols) +
    #bat
    geom_point(data = filter(aa_umap,grepl("bat",Class)), mapping = aes(x = UMAP1, y = UMAP2),col ="#631879FF",
               size = 0.5,
               #alpha = 0.5,
               stroke = 0) +
    theme_classic(base_line_size = 1, base_size = 20)+
    #human
    geom_point(data = filter(aa_umap,grepl("human",Class)), mapping = aes(x = UMAP1, y = UMAP2),col ="red",
               size = 0.5,
               #alpha = 0.5,
               stroke = 0)+
    theme(
      legend.position = "none",                     # 去掉图例
      axis.line.x = element_line(size = 1, color = "black"),
      axis.line.y = element_line(size = 1, color = "black")#,
      #axis.text = element_blank()   # 一次性移除 X 和 Y 轴所有刻度标签
    )
  #MPXV
  #geom_point(data = filter(aa_umap,grepl("GCF_000857045",Species)), mapping = aes(x = UMAP1, y = UMAP2), shape = 3, col = "black") +
  #MPXV
  #geom_point(data = filter(aa_umap,grepl("GCF_014621545",Species)), mapping = aes(x = UMAP1, y = UMAP2), shape = 3, col = "black") +
  #SCV2
  #geom_point(data = filter(aa_umap,grepl("oronavirus",Species)), mapping = aes(x = UMAP1, y = UMAP2), shape = 8, col = "red") +
  #geom_point(data = filter(aa_umap,grepl("GCF_009858895",Species)), mapping = aes(x = UMAP1, y = UMAP2), col = "red") +
}    
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/1A_UMAP.pdf', width= 4 , height= 4 , units='in')


aa_umap=fread("umap_plot2.txt")
Virus_CUB<-fread("EnsemblHost.virus.CUB.txt",header = F,sep = "\t")
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
Host_Virus=filter(Host_Virus,!Host_class %in% c("fly","elegans","yeast"))


#mid=Host_Virus %>% group_by(Species) %>% dplyr::summarise(count=n()) %>% filter(count>1)
#2441
#2784-332

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

RSCUij=fread("WijTable.txt")
head(RSCUij)
RSCUij$Codon=paste0(RSCUij$AA,":",RSCUij$Codon)
CUB<-reshape2::dcast(RSCUij,class+Organism_name~Codon,value.var = "RSCUij")
nrow(CUB)
CUB=na.omit(CUB)
colnames(CUB)[2]<-"Species"
head(CUB)

aa_umap=fread("umap_plot2.txt")
head(aa_umap)
aa2=merge(aa_umap[,1:2],CUB[,-1],by="Species")
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

CUB_df<-rbind(unique(filter(Virus_CUB,Organism_name %in% count$Organism_name )[,c("Organism_name","Codon","Num","class")]),
              Host_CUB[,c("Organism_name","Codon","Num","class")])
unique(CUB_df$class)
CUB<-reshape2::dcast(CUB_df,class+Organism_name~Codon,value.var = "Num")
#CUB=na.omit(CUB)
colnames(CUB)[2]<-"Species"
CUB$Sum<-rowSums(CUB[,-c(1:2)])
#CUB=filter(CUB,Sum>=100)
nrow(CUB[rowSums(CUB[,-c(1:2)])>=100,])
nrow(CUB)

test<-filter(CUB_df,!Codon %in% c("TGA","TAG","TAA"),Organism_name %in% CUB$Species)
codon_table<-fread("codon_table.txt",header=F)
colnames(codon_table)=c("Codon","AA3","AA","AAname")
test=merge(test,codon_table[,c(1,3)],by=c("Codon"))

test$First_base <- substr(test$Codon, 1, 1) #substr("ATG", 1, 1)
test$Second_base <- substr(test$Codon, 2, 2)#substr("ATG", 2, 2)
test$Third_base <- substr(test$Codon, 3, 3) #substr("ATG", 3, 3)
head(test)

Third_base = test %>% group_by(Organism_name,Third_base)  %>% dplyr::summarise(BaseNum=sum(Num))
length(unique(Third_base$Organism_name))
Third_base <- Third_base %>%
  mutate(ID = str_replace(Organism_name, "^.*:(GCF_[^:]+).*$", "\\1"))
head(Third_base)

aa_umap=fread("umap_plot2.txt")
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
letters <- multcompView::multcompLetters(p_matrix_sorted, compare="<=", threshold = 0.05)$Letters

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

First_base = test %>% group_by(Organism_name,First_base)  %>% dplyr::summarise(Base1Num=sum(Num))
Second_base = test %>% group_by(Organism_name,Second_base)  %>% dplyr::summarise(Base2Num=sum(Num))
Third_base = test %>% group_by(Organism_name,Third_base)  %>% dplyr::summarise(Base3Num=sum(Num))
base=merge(First_base,Second_base,by.x=c("Organism_name","First_base"),by.y=c("Organism_name","Second_base"))
head(base)
base=merge(base,Third_base,by.x=c("Organism_name","First_base"),by.y=c("Organism_name","Third_base"))
base$total=base$Base1Num+base$Base2Num+base$Base3Num
base = base %>% group_by(Organism_name,First_base)  %>% dplyr::summarise(BaseNum=sum(total))
colnames(base)[2]="base"

base <- base %>%
  mutate(ID = str_replace(Organism_name, "^.*:(GCF_[^:]+).*$", "\\1"))
head(base)

aa_umap=fread("umap_plot2.txt")
aa_umap <- aa_umap %>%
  mutate(ID = str_replace(Species, "^.*:(GCF_[^:]+).*$", "\\1"))

nrow(aa_umap)
head(aa_umap)
#Third_base[!Third_base$ID %in% aa_umap$ID,]
aa_umap[!aa_umap$ID %in% Third_base$ID,]

aa2=merge(aa_umap[,c(1,5,2)],base,by.x="ID",by.y="ID")


head(aa2)

base2 = aa2 %>% group_by(Species,Class)  %>% dplyr::summarise(Total=sum(BaseNum))
aa2=merge(aa2,base2,by=c("Species","Class"))
aa2$CG="CG"
aa2[aa2$base %in% c("T","A")]$CG="AT"
head(aa2)
base2 = aa2 %>% group_by(Species,Class,CG)  %>% dplyr::summarise(CG_Frac=sum(BaseNum/Total))
head(base2)
length(unique(base2$Species))

base2$Class=gsub(" ","",base2$Class)
base2[base2$Class=="human",]$Class="othermammals"

base2$Class=factor(base2$Class,levels = c("Liliopsida_Plant", "eudicotyledons_Plant","otherplants",
                                          "fish", "reptile", "bird",
                                          "othermammals", "bat",
                                          "RNAvirus", "DNAvirus"))


data=filter(base2,CG=="CG")
totalCG=data
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
    group1 <- data$CG_Frac[data$Class == classes[i]]
    group2 <- data$CG_Frac[data$Class == classes[j]]
    test_result <- wilcox.test(group1, group2)
    results <- rbind(results, data.frame(
      Group1 = as.character(classes[i]),
      Group2 = as.character(classes[j]),
      p_value = test_result$p.value
    ))
  }
}

# 创建一个对称的 p 值矩阵
p_matrix <- matrix(1, nrow = length(classes), ncol = length(classes))
rownames(p_matrix) <- levels(classes)
colnames(p_matrix) <- levels(classes)
for (k in 1:nrow(results)) {
  p_matrix[results$Group1[k], results$Group2[k]] <- results$p_value[k]
  p_matrix[results$Group2[k], results$Group1[k]] <- results$p_value[k]
}

# 使用 multcompView 进行字母标记
letters <- multcompView :: multcompLetters(p_matrix, compare="<=", threshold = 0.05)$Letters
letters_df <- data.frame(Class = names(letters), Letter = letters)


# 计算每个组的中位数
medians <- tapply(data$CG_Frac, data$Class, median)

# 按中位数大小对组别进行排序（从大到小）
sorted_classes <- names(sort(medians, decreasing = TRUE))  # 这里添加了decreasing = TRUE

# 重新排序p矩阵的行和列，使其与中位数排序一致
p_matrix_sorted <- p_matrix[sorted_classes, sorted_classes]

# 使用排序后的矩阵生成字母标记，确保字母顺序反映中位数大小（从大到小）
letters <- multcompView::multcompLetters(p_matrix_sorted, compare="<=", threshold = 0.05)$Letters

# 保持原始类别名称与字母的对应关系
letters_df <- data.frame(
  Class = names(letters), 
  Letter = as.character(letters),
  stringsAsFactors = FALSE
)

# 绘制箱线图并标记字母
ggplot(data, aes(x = Class, y = CG_Frac,col=as.factor(Class))) +
  geom_boxplot() +#
  #geom_jitter(mapping=aes(),alpha=0.2)+
  geom_text(data = letters_df, aes(x = Class, y = max(data$CG_Frac), label = Letter),color="black", vjust = -1) +
  theme_classic()+
  scale_color_manual(values =c("#006400","#00FF00","#b8ebb0",
                               "#B09C85FF","black","#FFC0CB",
                               "#ADB6B6FF","#631879FF",
                               "#57C3F3","blue"))+
  scale_y_continuous(breaks = c(seq(0.2,0.8,0.1)),limits = c(0.2,0.8)) #+
#my_theme

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/CG.pdf', width= 8 , height= 5 , units='in')

mid=merge(filter(Third_base2,CG=="CG")[,-3],filter(totalCG,CG=="CG")[,-3],by=c("Species","Class"))
head(mid)
ggplot(mid,aes(x=CG_Frac,y= CG3_Frac,col=Class))+
  geom_point()+
  #facet_grid(Sample~.,scales = "free")+
  theme_classic()+
  scale_x_continuous(breaks = c(seq(0,1,0.1)),limits = c(0,1)) +
  scale_y_continuous(breaks = c(seq(0,1,0.1)),limits = c(0,1)) +
  #labs(x="The RSCU values of host codons",y="The RSCU values of viral codons ")+
  geom_smooth(method = "glm",col="red",se=F)+
  geom_abline(slope = 1,intercept = 0)+
  #ggrepel::geom_text_repel(aes(label = Codon),color="blue")+
  ggpubr::stat_cor(method = "spearman",col="red",size=5)+
  #my_theme+
  scale_color_manual(values =c("#006400","#00FF00","#b8ebb0",
                               "#B09C85FF","black","#FFC0CB",
                               "#ADB6B6FF","#631879FF",
                               "#57C3F3","blue"))


#################################################################

#host VS virus CAI
#figure 1C

################################################################

#viral CAI

Wij<-read.table("WijTable.txt",header = T,sep="\t")#filter stop codons, log10
Wij$Host= sub("_.*", "", Wij$Organism_name)
Wij=filter(Wij,class!="Virus")
head(Wij)

Virus_CUB<-fread("EnsemblHost.virus.CUB.txt",header = F,sep = "\t")
colnames(Virus_CUB)=c("Host_class","Host","Organism_name","species","Codon","Num")
Virus_CUB$class="Virus"
head(Virus_CUB)
length(unique(Virus_CUB$Organism_name))
#filter(filter(aa_umap,grepl("irus",Class)), !Species %in% Virus_CUB$Organism_name)
#sum(filter(Virus_CUB,grepl("GCF_009858895",Organism_name),Host=="Sus")$Num)
#sum(filter(Virus_CUB,grepl("GCF_009858895",Organism_name),Host=="Sus",Codon %in% c("TAA","TAG","TGA"))$Num)
#sum(filter(Virus_CUB,grepl("GCF_009858895",Organism_name),Host=="Homo",Codon %in% c("TAA","TAG","TGA"))$Num)

CUB_virus<-filter(Virus_CUB,!Codon %in% c("TGA","TAG","TAA"),Organism_name %in% CUB$Species)
head(CUB_virus)
length(unique(filter(CUB_virus,Organism_name %in% aa_umap$Species)$Organism_name))

CAI<-data.frame()
head(CUB)
#filter(Host_Virus_CAI,grepl("GCF_000859885",Organism_name))

for (i in unique(CUB_virus$Host)) {
  CUB2=filter(Wij,Host==i)
  CUB_virus3=filter(CUB_virus,Host==i)
  for (j in unique(CUB2$Organism_name)) {
    CUB3<-filter(CUB2,Organism_name==j)
    for (virus in unique(CUB_virus3$Organism_name)) {
      CUB_virus2=filter(CUB_virus3,Organism_name==virus)
      codon_table<-merge(CUB_virus2,CUB3[,c("Codon","Wij")],by=c("Codon"))
      codon_table$CAI1<-(codon_table$Wij)*(codon_table$Num)
      CAI_mid<-codon_table %>% group_by(Organism_name) %>%  dplyr::summarise(CAI=10^(sum(CAI1)/sum(Num)),Host=j,CodonSum=sum(Num))
      CAI<-rbind(CAI_mid,CAI)
      print(paste(i,j,virus,sep="  "))
    }
  }
}
write.table(CAI,"Virus_CAI.txt",row.names = F,quote = F,sep = "\t")
length(unique(filter(CAI,Organism_name %in% aa_umap$Species)$Organism_name))
filter(aa_umap,grepl("virus",Class),!Species %in% CAI$Organism_name)$Species
length(unique(filter(aa_umap,grepl("virus",Class),!Species %in% CAI$Organism_name)$Species))
Virus_CAI=fread("Virus_CAI.txt")
head(Virus_CAI)
Virus_CUB<-fread("EnsemblHost.virus.CUB.txt",header = F,sep = "\t") 
colnames(Virus_CUB)=c("Host_class","Host","Organism_name","species","Codon","Num")
Virus_CUB=filter(Virus_CUB,!Codon %in% c("TAA","TAG","TGA"))
Virus_CUB2=unique(Virus_CUB[,c(3,4,5,6)]) %>% group_by(Organism_name) %>% dplyr::summarise(CodonSum=sum(Num))
Virus_CAI=merge(Virus_CAI[,1:3],Virus_CUB2,by=c("Organism_name"))
Virus_CAI[grepl("GCF_009858895.2",Organism_name)]$CodonSum=Virus_CAI[grepl("GCF_009858895.2",Organism_name)]$CodonSum/2
write.table(Virus_CAI,"Virus_CAI.txt",row.names = F,quote = F,sep = "\t")

#CAI=fread("Virus_CAI.txt")
#length(unique(CAI$Organism_name))

#host CAI
Wij<-read.table("WijTable.txt",header = T,sep="\t")[,c(2,3,7)]#filter stop codons, log10
head(Wij)
Wij=Wij[!grepl(":",Wij$Organism_name),]
CAI<-data.frame()
for (i in unique(Wij$Organism_name)) {
  CUB2=filter(Wij,Organism_name==i)
  print(i)
  print(paste0("LongestCDS_CUB/",i,".CUB.txt"))
  species=fread(paste0("LongestCDS_CUB/",i,".CUB.txt"),header=F,sep="\t")
  species$V1=toupper(species$V1)
  colnames(species)=c("Codon","Num","ID")
  codon_table<-merge(CUB2,species,by=c("Codon"))
  codon_table$CAI1<-(codon_table$Wij)*(codon_table$Num)
  CAI_mid<-codon_table %>% group_by(ID) %>%  dplyr::summarise(CAI=10^(sum(CAI1)/sum(Num)),Host=i,CodonSum=sum(Num))
  CAI<-rbind(CAI_mid,CAI)
  
}

head(CAI)
unique(CAI$Host)
CAI=fread("Host_CAI.txt") 

CAI_mid=filter(CAI,CodonSum>=100) %>% group_by(Host)  %>%  dplyr::summarise(median=median(CAI),low=quantile(CAI,0.025),high=quantile(CAI,0.975))
write.table(CAI_mid,"Host_CAI_median.txt",row.names = F,quote = F,sep = "\t")

{
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



species=fread("EnsemblHost.virus.CUB.species.txt",sep="\t",header = F)
colnames(species)[1:2]=c("Organism_name","Species")
head(species) 

mid=merge(Host_Virus_CAI,species,by=c("Organism_name","Species"))
mid$Host_Class <- factor(mid$Host_Class,
                         levels = c("other mammals","human","bat",
                                    "Liliopsida_Plant", "eudicotyledons_Plant","other plants",
                                    "fish", "reptile", "bird"))
ggplot(data = mid, 
       mapping = aes(x = median, y = CAI, color = Host_Class)) +
  
  # 先画底层点
  geom_point(data = filter(mid, Host_Class %in% c("other mammals","Liliopsida_Plant", "eudicotyledons_Plant","other plants", "reptile")),size=0.7) +
  geom_point(data = filter(mid, Host_Class %in% c("fish", "bird")),size=0.7) +  #  bat 放后面置顶
  
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
  geom_point(data = filter(mid, Host_Class == "bat"),size=0.7) +
  geom_point(data = filter(mid, Host_Class == "human"),size=0.7) +
  
  # 对角线
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
  
  #scale_y_continuous(breaks = seq(0.2,0.8,0.1), limit=c(0.27,0.8)) +
  #scale_x_continuous(breaks = seq(0.4,0.8,0.1), limit=c(0.49,0.8)) +
  ggpubr::stat_cor(method = "spearman",col="red",size=5,cor.coef.name = "rho")+
  # 颜色完全不变
  scale_color_manual(values = class_cols)+
  labs(x = "Host CAI median", y = "Viral CAI") +
  #labs(x = "", y = "") +
  theme_classic(base_line_size = 1.5, base_size = 20) +
  theme(
    #axis.text = element_blank(),
    #axis.ticks = element_blank(),
    legend.position = "none"
  )
length(unique(mid$Organism_name))




mid=fread("max_CAI_Virus_host2.txt")
#mid=fread("max_CAI_Virus_host2.cov.txt")
head(mid)
filter(mid,CodonSum>10000,Virus_Class=="RNA virus")
# 严格对应：9组 ↔ 9个颜色
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


#################################################################

#figure 1D

################################################################


mid=fread("max_CAI_Virus_host2.txt")
head(mid)
test=unique(mid[,c(2,4,6)]) %>% group_by(Host_Class,Host_genus) %>% dplyr::summarise(count=n()) %>% filter(count==1)
test=filter(mid,Host_genus %in% test$Host_genus) %>% group_by(Host_Class,Host_genus,Host) %>% dplyr::summarise(count=n())

result <- filter(test, !grepl("plants",Host_Class)) %>%
  group_by(Host_Class) %>%               # 按第1列分组
  slice_max(count, n = 1) %>%       # 每组挑 col4 最大的那一行
  pull(Host)  #Host_Class,                     # 只要第3列的值

violin_plot=filter(mid,Host %in% result)[,c("Host_Class","Host","CAI")]
violin_plot$class="virus"
colnames(violin_plot)
Host_CAI=fread("Host_CAI.txt")
Host_CAI=filter(Host_CAI,Host %in% result )
Host_CAI=merge(filter(Host_CAI,Host %in% result ),unique(mid[,c(2,4,6)]),by=c("Host"))
head(Host_CAI)
Host_CAI=Host_CAI[,c("Host_Class","Host","CAI")]
Host_CAI$class="host"
violin_plot=rbind(violin_plot,Host_CAI)
head(violin_plot)
violin_plot %>% group_by(Host_Class,class) %>% dplyr::summarise(count=n())
violin_plot$Host_Class=factor(violin_plot$Host_Class,levels = c( "human","bat","other mammals","bird", 
                                                                 "reptile","fish", #      "plants", 
                                                                 "Liliopsida_Plant", "eudicotyledons_Plant"#,"plants","other plants",
))

ggplot(violin_plot,
       aes(x = Host_Class, y = CAI, fill = class)) +
  geom_boxplot(outliers = F, alpha = .7)+
  scale_fill_manual(values = c("#E41A1C","#377EB8")) +
  theme_classic()
ggsave('Host_Viral_medianCAI2.pdf', width= 8 , height= 4 , units='in')

install.packages("gghalves")
library(gghalves)  
ggplot(violin_plot,
       aes(x = Host_Class, y = CAI, fill = class)) +
  
  # RNA 左半边
  geom_half_violin(
    data = . %>% filter(class == "host"),
    side = "l", alpha = .7, color = NA
  ) +
  # geom_half_boxplot(                     # ← 半个 box
  #   data = . %>% filter(class == "host"),
  #   side = "l", width = .2, outlier.shape = NA, color = "black"
  # ) +
  
  # DNA 右半边
  geom_half_violin(
    data = . %>% filter(class == "virus"),
    side = "r", alpha = .7, color = NA
  ) +
  # geom_half_boxplot(                     # ← 半个 box
  #   data = . %>% filter(class == "virus"),
  #   side = "r", width = .2, outlier.shape = NA, color = "black"
  # ) +
  
  scale_fill_manual(values = c("#E41A1C","#377EB8")) +
  #labs(x = "Host class", y = "CAI", fill = "Virus class") +
  theme_classic()
#theme(axis.text.x = element_text(angle = 45,hjust = 1,color="black"))

# ggplot(violin_plot,
#        aes(x =  class, y = CAI, fill = class)) +
#   
#   geom_boxplot( ) +
#   scale_fill_manual(values = c("#E41A1C","#377EB8")) +
#   #labs(x = "Host class", y = "CAI", fill = "Virus class") +
#   theme_classic()+
#   facet_grid(~Host_Class)+
#   ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = list(c("host","virus")),label = 'p.signif')
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Host_Viral_medianCAI.pdf', width= 14 , height= 8 , units='in')

library(ggplot2)
library(dplyr)
library(gghalves)   # 需安装：install.packages("gghalves")

# 假设您的数据框名为 violin_plot，包含三列：Host_Class, CAI, class
# class 取值为 "host" 和 "virus"
violin_plot$Host_Class = factor(violin_plot$Host_Class,levels = rev(levels(violin_plot$Host_Class)))

# 1. 计算各分组的中位数
median_df <- violin_plot %>%
  group_by(Host_Class, class) %>%
  summarise(median_CAI = median(CAI, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    x_pos = as.numeric(factor(Host_Class, levels = unique(Host_Class))),
    xmin = ifelse(class == "host", x_pos - 0.4, x_pos),
    xmax = ifelse(class == "host", x_pos, x_pos + 0.4)
  ) 
# 2. 绘图
ggplot(violin_plot, aes(x = Host_Class, y = CAI, fill = class)) +
  # 左半（host）
  geom_half_violin(
    data = filter(violin_plot, class == "host"),
    side = "l", alpha = 0.7, color = NA
  ) +
  # 右半（virus）
  geom_half_violin(
    data = filter(violin_plot, class == "virus"),
    side = "r", alpha = 0.7, color = NA
  ) +
  # 中位数线段
  geom_segment(
    data = median_df,
    aes(x = xmin, xend = xmax, y = median_CAI, yend = median_CAI),
    color = "black", size = 1.2, lineend = "round"
  ) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  theme_classic() +
  ylim(0.2,1)+
  labs(x = "Host class", y = "CAI", fill = "Virus class")+
  coord_flip()  
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Host_Viral_medianCAI3.pdf', width= 8 , height= 14 , units='in')
