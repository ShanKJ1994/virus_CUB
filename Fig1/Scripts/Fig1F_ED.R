{
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  se <- function(x) sd(x)/sqrt(length(x))
  
}

#/lustre/user/lulab/shankj/miniconda3/envs/monocle_env/bin/R
{
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  library(parallel)        # 已随 R 自带
  library(data.table)
}

#################################################################

#RSCU的欧几里得距离

################################################################

RSCUij=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/WijTable.txt")
head(RSCUij)
RSCUij$Codon=paste0(RSCUij$AA,":",RSCUij$Codon)
CUB<-reshape2::dcast(RSCUij,class+Organism_name~Codon,value.var = "RSCUij")
nrow(CUB)
CUB=na.omit(CUB)
colnames(CUB)[2]<-"Species"
head(CUB)

aa_umap=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/umap_plot2.txt")
head(aa_umap)
aa2=merge(aa_umap[,1:2],CUB[,-1],by="Species")
aa2$Class=as.character(aa2$Class)



# 把数据矩阵化（只保留数值）
mat <- as.matrix(aa2[, -c("Species", "Class")])
rownames(mat) <- aa2$Species

# 多线程计算
n.cores <- max(1, 30)   # 留 2 个核给系统
pair.idx <- combn(nrow(mat), 2, simplify = FALSE)

dist.dt <- mclapply(pair.idx, function(idx) {
  i <- idx[1]; j <- idx[2]
  d <- sqrt(sum((mat[i, ] - mat[j, ])^2))
  data.table(species1 = rownames(mat)[i],
             species2 = rownames(mat)[j],
             Euclidean_Distance = d)
}, mc.cores = n.cores) |> rbindlist()

fwrite(dist.dt, "/gpfs2/shankj/CodonUsage2/Ensembl/20250703/euclidean_distance.txt",
       sep = "\t", quote = FALSE)

#RSCUij=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/WijTable.txt")
#head(RSCUij)
#RSCUij$Codon=paste0(RSCUij$AA,":",RSCUij$Codon)
#CUB<-reshape2::dcast(RSCUij,class+Organism_name~Codon,value.var = "RSCUij")
#nrow(CUB)
#CUB=na.omit(CUB)
#colnames(CUB)[2]<-"Species"
#head(CUB)

#aa_umap=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/umap_plot2.txt")
#head(aa_umap)
#aa2=merge(aa_umap[,1:2],CUB[,-1],by="Species")
#aa2$Class=as.character(aa2$Class)

## 计算每两行欧几里得距离的函数
#calculate_euclidean_distance <- function(df) {
  #df=data.frame(df)
  ## 获取数据行数
  #n_rows <- nrow(df)
  
  ## 数值列索引（从第3列到最后一列）
  #value_cols_idx <- 3:ncol(df)

  #result_df=data.frame()
  
  ## 对每两行计算欧几里得距离
  #for (i in seq(1, n_rows-1, by = 1)) {
    #for (j in seq(i+1, n_rows, by = 1)) {
        ## 提取两行的注释信息
        #species1_class <- df[i, 2]
        #species1_name <- df[i, 1]
        #species2_class <- df[j, 2]
        #species2_name <- df[j, 1]
        
        ## 提取两行的所有RSCU值
        #row1_values <- as.numeric(df[i, value_cols_idx])
        #row2_values <- as.numeric(df[j, value_cols_idx])
        
        ## 计算欧几里得距离：√[Σ(x₁i - x₂i)²]
        #distance <- sqrt(sum((row1_values - row2_values)^2))
        
        #cat(paste0(species1_name,"\tvs\t",species2_name,"\n"))
        ## 添加到结果数据框
        #result_df=rbind(result_df,
        #data.frame(
          ##species1_class=species1_class, 
          #species1_name=species1_name, 
          ##species2_class=species2_class, 
          #species2_name=species2_name, 
          #Euclidean_Distance=distance
        #)
        #)
    #}
  #}
  
  #return(result_df)
#}

## 调用函数计算欧几里得距离
#result <- calculate_euclidean_distance(aa2)

#write.table(result,"/gpfs2/shankj/CodonUsage2/Ensembl/20250703/euclidean_distance.txt",row.names = F,quote = F,sep = "\t")






ED=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/euclidean_distance.txt")
ED$ID1 <- stringr::word(ED$species1, -1, sep = ":")
ED$ID2 <- stringr::word(ED$species2, -1, sep = ":")
head(ED)
##################

#Genus

#################

Host_virus=fread("/gpfs2/shankj/CodonUsage2/Ensembl/VirusHost_allhost.txt",header=F)
Host_virus$V2=tolower(Host_virus$V2)
#Mammal Genus from NCBI https://www.ncbi.nlm.nih.gov/taxonomy/?term=txid40674[Organism:noexp]
Genus=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/Mammalia_cleaned.txt",header=F,sep = "\t")
Mammalia <- tolower(sapply(strsplit(Genus$V1, " "), `[`, 1))
Host_virus=filter(Host_virus,V2 %in% Mammalia)

#Host_virus=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/EnsemblHost.virus.CUB.txt",header=F)
#Host_virus=unique(Host_virus[,c(1,2,3)])
Host_virus$ID <- stringr::word(Host_virus$V3, -1, sep = ":")

head(Host_virus)


class=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/umap_plot2.txt")
virus=filter(class, grepl("RNA|DNA",Class))
virus$Class="Virus"
virus$ID <- stringr::word(virus$Species, -1, sep = ":")

virus_ED=filter(ED,ID1 %in% virus$ID, ID2 %in% virus$ID,ID1 %in%  Host_virus$ID, ID2 %in%  Host_virus$ID)
nrow(virus_ED)
virus_ED$barcode=paste0(virus_ED$species1,"_vs_",virus_ED$species2)
virus_ED$class="diff"


data=data.frame()
for (host in unique(Host_virus$V2)) {
   genus=filter(virus_ED, ID1 %in% filter(Host_virus, V2==host)$ID , ID2 %in% filter(Host_virus, V2==host)$ID)
   genus$genus=host
   genus$class="Same"
   genus$barcode=paste0(genus$species1,"_vs_",genus$species2)
   data=rbind(data,genus)
}

virus_ED[virus_ED$barcode %in% data$barcode,]$class="same"
head(virus_ED)

virus_ED %>% group_by(class) %>% dplyr::summarise(virus_pair_count=n(),Mean=mean(Euclidean_Distance),sd=sd(Euclidean_Distance),Median=median(Euclidean_Distance),
                                                  Q0.025=quantile(Euclidean_Distance,0.025),Q0.975=quantile(Euclidean_Distance,0.975),
                                                  Q0.25=quantile(Euclidean_Distance,0.25),Q0.75=quantile(Euclidean_Distance,0.75))


Homo_same=filter(virus_ED, ID1 %in% filter(Host_virus, V2=="homo")$ID , ID2 %in% filter(Host_virus, V2=="homo")$ID)

Homo_diff=rbind(filter(virus_ED, class=="diff",ID1 %in% filter(Host_virus, V2=="homo")$ID  , !ID2 %in% filter(Host_virus, V2=="homo")$ID),
                filter(virus_ED, class=="diff",!ID1 %in% filter(Host_virus, V2=="homo")$ID  , ID2 %in% filter(Host_virus, V2=="homo")$ID))

Homo=rbind(Homo_same,Homo_diff)


Homo %>% group_by(class) %>% dplyr::summarise(virus_pair_count=n(),Mean=mean(Euclidean_Distance),sd=sd(Euclidean_Distance),Median=median(Euclidean_Distance),
                                               Q0.025=quantile(Euclidean_Distance,0.025),Q0.975=quantile(Euclidean_Distance,0.975),
                                               Q0.25=quantile(Euclidean_Distance,0.25),Q0.75=quantile(Euclidean_Distance,0.75))

#ggplot(Homo,aes(x=class,y=Euclidean_Distance)) +
  #geom_violin(mapping=aes(fill=class))+
  #geom_boxplot(width=0.3,col="black")+
  ##geom_jitter(alpha=0.01)+
  ##facet_wrap(~type,scales = "free")+
  ##my_theme+
  #theme_classic()
#ggsave('/gpfs2/shankj/CodonUsage2/Ensembl/20250703/Homo_Euclidean_Distance.pdf', width= 3 , height= 5 , units='in')
write.table(Homo,"/gpfs2/shankj/CodonUsage2/Ensembl/20250703/Homo_Euclidean_Distance.txt",row.names = F,quote = F,sep = "\t")


library(gghalves)             
ggplot(Homo,
       aes(x = Type, y = Euclidean_Distance, fill = class)) +
  
  # RNA 左半边
  geom_half_violin(
    data = . %>% filter(class == "same"),
    side = "l", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "same"),
    side = "l", width = .2, outlier.shape = NA, color = "black"
  ) +
  
  # DNA 右半边
  geom_half_violin(
    data = . %>% filter(class == "diff"),
    side = "r", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "diff"),
    side = "r", width = .2, outlier.shape = NA, color = "black"
  ) +
  scale_y_continuous(breaks = c(seq(2.5,6.5,1)),limits = c(1,7)) +
  scale_fill_manual(values = c( "#968175","#AB3282")) +
  #labs(x = "Host class", y = "CAI", fill = "Virus class") +
  theme_classic()+
  coord_cartesian(ylim = c(2.5,6.5))





nonHomo_same=filter(virus_ED, class=="same",!ID1 %in% filter(Host_virus, V2=="homo")$ID , !ID2 %in% filter(Host_virus, V2=="homo")$ID)
nonHomo_diff=rbind(filter(virus_ED, class=="diff",!ID1 %in% filter(Host_virus, V2=="homo")$ID , !ID2 %in% filter(Host_virus, V2=="homo")$ID))

nonHomo=rbind(nonHomo_same,nonHomo_diff) %>% filter(!barcode %in% Homo$barcode)


nonHomo %>% group_by(class) %>% dplyr::summarise(virus_pair_count=n(),Mean=mean(Euclidean_Distance),sd=sd(Euclidean_Distance),Median=median(Euclidean_Distance),
                                               Q0.025=quantile(Euclidean_Distance,0.025),Q0.975=quantile(Euclidean_Distance,0.975),
                                               Q0.25=quantile(Euclidean_Distance,0.25),Q0.75=quantile(Euclidean_Distance,0.75))


wilcox.test(filter(nonHomo ,class=="same")$Euclidean_Distance,filter(nonHomo ,class=="diff")$Euclidean_Distance)
wilcox.test(filter(Homo ,class=="same")$Euclidean_Distance,filter(Homo ,class=="diff")$Euclidean_Distance)

write.table(nonHomo,"/gpfs2/shankj/CodonUsage2/Ensembl/20250703/nonHomoMammal_Euclidean_Distance.txt",row.names = F,quote = F,sep = "\t")

ggplot(nonHomo,
       aes(x = Type, y = Euclidean_Distance, fill = class)) +
  
  # RNA 左半边
  geom_half_violin(
    data = . %>% filter(class == "same"),
    side = "l", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "same"),
    side = "l", width = .2, outlier.shape = NA, color = "black"
  ) +
  
  # DNA 右半边
  geom_half_violin(
    data = . %>% filter(class == "diff"),
    side = "r", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "diff"),
    side = "r", width = .2, outlier.shape = NA, color = "black"
  ) +
  scale_y_continuous(breaks = c(seq(2.5,6.5,1)),limits = c(1,7)) +
  scale_fill_manual(values = c( "#968175","#AB3282")) +
  #labs(x = "Host class", y = "CAI", fill = "Virus class") +
  theme_classic()+
  coord_cartesian(ylim = c(2.5,6.5))


##################

#Species

#################
Host_virus=fread("/gpfs2/shankj/CodonUsage2/Ensembl/VirusHost_allhost_matchedSpecies.txt",header=F)
Host_virus$V2=tolower(Host_virus$V2)
#Mammal Genus from NCBI https://www.ncbi.nlm.nih.gov/taxonomy/?term=txid40674[Organism:noexp]
Genus=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/Mammalia_cleaned.txt",header=F,sep = "\t")
Mammalia <- sapply(strsplit(Genus$V1, " "), function(x) {
  tolower(paste(head(x, 2), collapse = "_"))
})
Host_virus=filter(Host_virus,V1 %in% Mammalia)

#Host_virus=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/EnsemblHost.virus.CUB.txt",header=F)
#Host_virus=unique(Host_virus[,c(1,2,3)])
Host_virus$ID <- stringr::word(Host_virus$V3, -1, sep = ":")

#head(Host_virus)




class=fread("/gpfs2/shankj/CodonUsage2/Ensembl/20250703/umap_plot2.txt")
virus=filter(class, grepl("RNA|DNA",Class))
virus$Class="Virus"
virus$ID <- stringr::word(virus$Species, -1, sep = ":")

virus_ED=filter(ED,ID1 %in% virus$ID, ID2 %in% virus$ID,ID1 %in%  Host_virus$ID, ID2 %in%  Host_virus$ID)
nrow(virus_ED)
virus_ED$barcode=paste0(virus_ED$species1,"_vs_",virus_ED$species2)
virus_ED$class="diff"


data=data.frame()
for (host in unique(Host_virus$V2)) {
   genus=filter(virus_ED, ID1 %in% filter(Host_virus, V2==host)$ID , ID2 %in% filter(Host_virus, V2==host)$ID)
   genus$genus=host
   genus$class="Same"
   genus$barcode=paste0(genus$species1,"_vs_",genus$species2)
   data=rbind(data,genus)
}

virus_ED[virus_ED$barcode %in% data$barcode,]$class="same"
head(virus_ED)

virus_ED %>% group_by(class) %>% dplyr::summarise(virus_pair_count=n(),Mean=mean(Euclidean_Distance),sd=sd(Euclidean_Distance),Median=median(Euclidean_Distance),
                                                  Q0.025=quantile(Euclidean_Distance,0.025),Q0.975=quantile(Euclidean_Distance,0.975),
                                                  Q0.25=quantile(Euclidean_Distance,0.25),Q0.75=quantile(Euclidean_Distance,0.75))


Homo_same=filter(virus_ED, ID1 %in% filter(Host_virus, V1=="homo_sapiens")$ID , ID2 %in% filter(Host_virus, V1=="homo_sapiens")$ID)

Homo_diff=rbind(filter(virus_ED, class=="diff",ID1 %in% filter(Host_virus, V1=="homo_sapiens")$ID  , !ID2 %in% filter(Host_virus, V1=="homo_sapiens")$ID),
                filter(virus_ED, class=="diff",!ID1 %in% filter(Host_virus, V1=="homo_sapiens")$ID  , ID2 %in% filter(Host_virus, V1=="homo_sapiens")$ID))

Homo=rbind(Homo_same,Homo_diff)


Homo %>% group_by(class) %>% dplyr::summarise(virus_pair_count=n(),Mean=mean(Euclidean_Distance),sd=sd(Euclidean_Distance),Median=median(Euclidean_Distance),
                                               Q0.025=quantile(Euclidean_Distance,0.025),Q0.975=quantile(Euclidean_Distance,0.975),
                                               Q0.25=quantile(Euclidean_Distance,0.25),Q0.75=quantile(Euclidean_Distance,0.75))

#ggplot(Homo,aes(x=class,y=Euclidean_Distance)) +
  #geom_violin(mapping=aes(fill=class))+
  #geom_boxplot(width=0.3,col="black")+
  ##geom_jitter(alpha=0.01)+
  ##facet_wrap(~type,scales = "free")+
  ##my_theme+
  #theme_classic()
#ggsave('/gpfs2/shankj/CodonUsage2/Ensembl/20250703/Homo_Euclidean_Distance.pdf', width= 3 , height= 5 , units='in')
write.table(Homo,"/gpfs2/shankj/CodonUsage2/Ensembl/20250703/Homo_Euclidean_Distance.txt",row.names = F,quote = F,sep = "\t")


library(gghalves)             
ggplot(Homo,
       aes(x = Type, y = Euclidean_Distance, fill = class)) +
  
  # RNA 左半边
  geom_half_violin(
    data = . %>% filter(class == "same"),
    side = "l", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "same"),
    side = "l", width = .2, outlier.shape = NA, color = "black"
  ) +
  
  # DNA 右半边
  geom_half_violin(
    data = . %>% filter(class == "diff"),
    side = "r", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "diff"),
    side = "r", width = .2, outlier.shape = NA, color = "black"
  ) +
  scale_y_continuous(breaks = c(seq(2.5,6.5,1)),limits = c(1,7)) +
  scale_fill_manual(values = c( "#968175","#AB3282")) +
  #labs(x = "Host class", y = "CAI", fill = "Virus class") +
  theme_classic()+
  coord_cartesian(ylim = c(2.5,6.5))





nonHomo_same=filter(virus_ED, class=="same",!ID1 %in% filter(Host_virus, V1=="homo_sapiens")$ID , !ID2 %in% filter(Host_virus, V1=="homo_sapiens")$ID)
nonHomo_diff=rbind(filter(virus_ED, class=="diff",!ID1 %in% filter(Host_virus, V1=="homo_sapiens")$ID , !ID2 %in% filter(Host_virus, V1=="homo_sapiens")$ID))

nonHomo=rbind(nonHomo_same,nonHomo_diff) %>% filter(!barcode %in% Homo$barcode)


nonHomo %>% group_by(class) %>% dplyr::summarise(virus_pair_count=n(),Mean=mean(Euclidean_Distance),sd=sd(Euclidean_Distance),Median=median(Euclidean_Distance),
                                               Q0.025=quantile(Euclidean_Distance,0.025),Q0.975=quantile(Euclidean_Distance,0.975),
                                               Q0.25=quantile(Euclidean_Distance,0.25),Q0.75=quantile(Euclidean_Distance,0.75))


wilcox.test(filter(nonHomo ,class=="same")$Euclidean_Distance,filter(nonHomo ,class=="diff")$Euclidean_Distance)
wilcox.test(filter(Homo ,class=="same")$Euclidean_Distance,filter(Homo ,class=="diff")$Euclidean_Distance)

write.table(nonHomo,"/gpfs2/shankj/CodonUsage2/Ensembl/20250703/nonHomoMammal_Euclidean_Distance.txt",row.names = F,quote = F,sep = "\t")

ggplot(nonHomo,
       aes(x = Type, y = Euclidean_Distance, fill = class)) +
  
  # RNA 左半边
  geom_half_violin(
    data = . %>% filter(class == "same"),
    side = "l", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "same"),
    side = "l", width = .2, outlier.shape = NA, color = "black"
  ) +
  
  # DNA 右半边
  geom_half_violin(
    data = . %>% filter(class == "diff"),
    side = "r", alpha = .4, color = NA
  ) +
  geom_half_boxplot(                     # ← 半个 box
    data = . %>% filter(class == "diff"),
    side = "r", width = .2, outlier.shape = NA, color = "black"
  ) +
  scale_y_continuous(breaks = c(seq(2.5,6.5,1)),limits = c(1,7)) +
  scale_fill_manual(values = c( "#968175","#AB3282")) +
  #labs(x = "Host class", y = "CAI", fill = "Virus class") +
  theme_classic()+
  coord_cartesian(ylim = c(2.5,6.5))