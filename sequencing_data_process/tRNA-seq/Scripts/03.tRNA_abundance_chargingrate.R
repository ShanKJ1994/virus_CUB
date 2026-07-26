se <- function(x) sd(x)/sqrt(length(x))

# 定义函数来实现反向互补
reverse_complement <- function(sequence) {
  # 定义碱基互补映射
  complement_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  # 将输入序列转换为单个碱基的向量
  bases <- strsplit(sequence, "")[[1]]
  
  # 获取互补碱基
  complement_bases <- complement_map[bases]
  
  # 反转互补碱基向量
  reverse_complement_bases <- rev(complement_bases)
  
  # 将反转后的互补碱基向量组合成字符串
  reverse_complement_sequence <- paste(reverse_complement_bases, collapse = "")
  
  return(reverse_complement_sequence)
}

# 定义处理向量的函数
reverse_complement_vector <- function(seq_vector) {
  result <- sapply(seq_vector, reverse_complement)
  return(result)
}

# 示例向量
sequence_vector <- c("ATCG", "GCAT", "CCGG")

# 调用函数进行反向互补操作
result_vector <- reverse_complement_vector(sequence_vector)

# 输出结果
print("原始序列向量:")
print(sequence_vector)
print("反向互补序列向量:")
print(result_vector)

############################################################

#tRNA charge rate

############################################################

CCAprops=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/CCAprops.csv")
CCAprops=filter(CCAprops,!grepl("eColi",gene),!grepl("mito",gene)) %>% 
  separate(gene,into = c('AA','anti-codon','Isodecoder'),sep = '-') %>% filter(!AA %in% c("tRX","Sup"),!`anti-codon` %in% c("ACT"))
CCAprops$Codon=reverse_complement_vector(CCAprops$`anti-codon`)
head(CCAprops)
unique(CCAprops$Codon)
unique(CCAprops$`anti-codon`)
length(unique(CCAprops$`anti-codon`))
nrow(unique(CCAprops[,c(1,2)]))
#unique(filter(CCAprops, !`anti-codon` %in% FSJ$`anti-codon`)[,c(1,2)])
#length(FSJ$`anti-codon`)

CCAprops=CCAprops %>% group_by(AA,`anti-codon`,Codon,end,sample) %>% dplyr::summarise(count=sum(count)) 
CCAprops$sample=gsub("/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj/clean_","",CCAprops$sample)
CCAprops$sample=gsub("_R1.fq.gz.unpaired_uniq.bam","",CCAprops$sample)
CCAprops$Charge=CCAprops$end
CCAprops[CCAprops$Charge!="CA",]$Charge="No"
CCAprops[CCAprops$Charge=="CA",]$Charge="Yes"
CCAprops=CCAprops %>% group_by(`anti-codon`,AA,Charge,sample) %>% dplyr::summarise(count=sum(count)) 

mid=CCAprops %>% group_by(`anti-codon`,AA,sample) %>% dplyr::summarise(sum=sum(count)) 
head(mid)
CCAprops=merge(CCAprops,mid,by=c("sample","AA","anti-codon"))
CCAprops$frac=CCAprops$count/CCAprops$sum
CCAprops$Sample=gsub("_\\d+","",CCAprops$sample)
CCAprops=filter(CCAprops,Charge=="Yes")
head(CCAprops)

CCAprops$Gene=paste0(CCAprops$AA,"-",CCAprops$`anti-codon`)
mid=CCAprops %>% group_by(Sample,Gene) %>% dplyr::summarise(mean=mean(frac),se=se(frac))
mid=reshape2::dcast(mid,Gene~Sample,value.var = "mean")
rownames(mid)=mid$Gene
mid=mid[,-1]
gene_order_sorted
my.colors <- colorRampPalette(c("grey","white", "red"))(50)
gene_order_sorted=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/gene_order_sorted.csv")

plt=pheatmap::pheatmap(t(mid)[,gene_order_sorted], angle_col = 90,               # 你的数值矩阵
                       color = my.colors,cluster_rows = F,cluster_cols = F,  # 自定义渐变
                       breaks = seq(0.4,0.9, 
                                    length.out = 50))  # 保证对称
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/tRNA_charge.pdf', plt, width= 8 , height= 3 , units='in')


codon_table<-fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/UMAP/codon_table.txt",header=F)
colnames(codon_table)=c("Codon","AA3","AA","AAname")
mid1=filter(CCAprops,Charge=="Yes")
head(mid1)
mid1$timepoint=0
mid1[grepl("X1",mid1$sample),]$timepoint=1
mid1[grepl("X2",mid1$sample),]$timepoint=2
mid1[grepl("X4",mid1$sample),]$timepoint=3
#mid1=filter(mid1,grepl("X",mid1$sample))
head(mid1)
#mid1=merge(unique(codon_table[,2:3]),mid1,by.x="AA3",by.y="AA")
name="Glu"
mid1$AA=paste0(mid1$AA,"-",mid1$`anti-codon`)
mid1$replicate=gsub("Ctrl_","",mid1$sample)
mid1$replicate=gsub("X\\d+_","",mid1$replicate)
mid1$replicate=as.numeric(mid1$replicate)
str(mid1)
df=data.frame(AA=unique(mid1[,c(2)]),pvalue=100)
str(new)
colnames(mid1)
#name="C-GCA"
for (name in unique(mid1$AA)){
  new=filter(mid1,AA==name)[,c(9,10,5,6)]
  tryCatch({
    #拟合负二项回归模型
    model <- MASS::glm.nb(
      count ~ timepoint + offset(log(sum)),
      data = new,
      control = glm.control(maxit = 500)
    )
    # 使用 drop1 进行模型简化并获取 p 值
    drop1_result <- drop1(model, test = "Chisq")
    p_value <- drop1_result["timepoint", "Pr(>Chi)"]
    
    # 将 p 值赋给原始表格的最后一列
    df[df$AA == name, "pvalue"] <- p_value
  },error =function(e){
    #捕获错误并打印警告，但继续执行其他操作
    warning(sprintf("Error processing %s:%s",name, e$message))
  })
}
df$FDR=p.adjust(df$pvalue,method = "fdr")
filter(df,FDR<0.05)
df$sig="n.s."
df[df$FDR<0.05,]$sig="*"
df[df$FDR<0.01,]$sig="**"
df[df$FDR<0.001,]$sig="***"
mid1 = mid1 %>% group_by(AA,Charge,Sample) %>% dplyr::summarise(Count=mean(count),Count_se=se(count),
                                                                Frac=mean(frac),Frac_se=se(frac)) 

# 创建初始向量
vec <- unique(mid1$Sample)
# 组合所有可能的向量
comb <- combn(vec, 2)
# 将所有向量形成一个列表
my_pair <- lapply(seq_len(ncol(comb)), function(i) comb[,i])
ggplot(mid1,aes(x=Sample,y=Frac,col=Sample)) +#,grepl("\\w",external_gene_name)
  #geom_violin(mid2,mapping=aes(x=class,y=log10(mean),fill = class))+
  geom_boxplot(width=0.5,col="black")+
  geom_point()+
  #facet_wrap(~type,scales = "free")+
  #my_theme+
  theme_classic()+
  scale_y_continuous(breaks = c(seq(0.4,0.9,0.1)),limits = c(0.4,0.9)) +
  geom_line(aes(group=AA), linetype="dashed", col="grey",size=0.05)+
  geom_line(filter(mid1, AA %in% filter(df,FDR<0.05)$AA),mapping=aes(group=AA), linetype="dashed", col="red",size=0.3)+
  labs(y="AA Charge rate")+
  ggrepel::geom_text_repel(filter(mid1,Sample=="Ctrl", AA %in% filter(df,FDR<0.05)$AA),mapping=aes(x="Ctrl",y=Frac,label = AA),color="#FF4040")+
  #ggpubr::stat_compare_means(mapping = aes(group = Sample),method = 'wilcox.test', comparisons = my_pair,paired = T,label = 'p.signif' ) +
  scale_color_manual(values=c("#B7B5B6", "#B7D0EC", "#89B3D9", "#408FBF"))
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Charge_tRNA_point.pdf', width= 8 , height= 5 , units='in')

# df[df$sig=="n.s.",]$sig=""
# ggplot(merge(mid1,df,by=c("AA")),
#        aes(x=AA,y=Frac,fill=Sample)) +
#   geom_bar(position="dodge", stat="identity", colour="black")+
#   scale_y_continuous(breaks = c(seq(0,0.9,0.1)),limits = c(0,0.9)) +
#   #my_theme+
#   scale_fill_manual(values=c("#B7B5B6", "#B7D0EC", "#89B3D9", "#408FBF")) +
#   theme_classic()+
#   geom_text(aes(x=AA,y=0.9,label = sig))+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1,color="black"))+
#   labs(title = paste0("tRNA charged rate"),y="charged rate")#+
# # geom_errorbar(aes(ymin=(Frac-Frac_se),ymax=(Frac+Frac_se)),width=.2,col='red',size=.2,position=position_dodge(.5))
# ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Charge_tRNA_20250709.pdf', width= 10 , height= 5 , units='in')
# 
########################

#total tRNA charge rate

########################

mid1=filter(CCAprops,Charge=="Yes")
head(mid1)
mid1=mid1 %>% group_by(sample,Sample) %>% dplyr::summarise(count=sum(count),sum=sum(sum))
mid1$timepoint=0
mid1[grepl("X1",mid1$sample),]$timepoint=1
mid1[grepl("X2",mid1$sample),]$timepoint=2
mid1[grepl("X4",mid1$sample),]$timepoint=3
#mid1=filter(mid1,grepl("X",mid1$sample))
mid1

mid1$replicate=gsub("Ctrl_","",mid1$sample)
mid1$replicate=gsub("X\\d+_","",mid1$replicate)
mid1$replicate=as.numeric(mid1$replicate)
str(mid1)
df=data.frame(pvalue=100)
str(new)
#for (name in unique(mid1$AA)){
new=mid1[,c(4,5,2,3)]
tryCatch({
  #拟合负二项回归模型
  model <- MASS::glm.nb(
    count ~ timepoint + offset(log(sum)),
    data = new,
    control = glm.control(maxit = 500)
  )
  # 使用 drop1 进行模型简化并获取 p 值
  drop1_result <- drop1(model, test = "Chisq")
  p_value <- drop1_result["timepoint", "Pr(>Chi)"]
  
  # 将 p 值赋给原始表格的最后一列
  df[, "pvalue"] <- p_value
},error =function(e){
  #捕获错误并打印警告，但继续执行其他操作
  warning(sprintf("Error processing %s:%s",name, e$message))
})
#}

mid = mid1 %>% group_by(Sample) %>% dplyr::summarise(Frac=mean(count/sum),Frac_se=se(count/sum)) 

ggplot(mid,aes(x=Sample,y=Frac,col=Sample)) +
  geom_point(size=5)+
  geom_line(group=1)+
  #facet_wrap(~type,scales = "free")+
  #my_theme+
  theme_classic()+
  scale_y_continuous(breaks = c(seq(0.65,0.8,0.05)),limits = c(0.65,0.8)) +
  labs(y="tRNA Charge rate")+
  scale_color_manual(values=c("#B7B5B6", "#B7D0EC", "#89B3D9", "#408FBF"))+
  geom_errorbar(aes(ymin = Frac - Frac_se, ymax = Frac + Frac_se), 
                width = 0.2, position = position_dodge(0.1) )

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/ChargeTotal_point_20250709.pdf', width= 6 , height= 5 , units='in')


############################################################

#AA charge rate

############################################################
mid=CCAprops %>% group_by(AA,Charge,Sample,sample) %>% dplyr::summarise(Count=sum(count),Sum=sum(sum)) 
mid=mid %>% group_by(Sample,AA) %>% dplyr::summarise(mean=mean(Count/Sum),se=se(Count/Sum))
mid=reshape2::dcast(mid,AA~Sample,value.var = "mean")
rownames(mid)=mid$AA
mid=mid[,-1]
my.colors <- colorRampPalette(c("grey","white", "red"))(50)
plt=pheatmap::pheatmap(t(mid), angle_col = 90,               # 你的数值矩阵
                       color = my.colors,cluster_rows = F,cluster_cols = F,  # 自定义渐变
                       breaks = seq(0.4, 0.9, 
                                    length.out = 50))  # 保证对称
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/AA_charge.pdf', plt, width= 6 , height= 3 , units='in')



mid1=CCAprops %>% group_by(AA,Charge,sample) %>% dplyr::summarise(Count=sum(count),Sum=sum(sum)) 
head(mid1)
#mid1$sample=mid1$Sample
mid1$timepoint=0
mid1[grepl("X1",mid1$sample),]$timepoint=1
mid1[grepl("X2",mid1$sample),]$timepoint=2
mid1[grepl("X4",mid1$sample),]$timepoint=3
#mid1=filter(mid1,grepl("X",mid1$sample))
head(mid1)
name="Glu"
mid1$replicate=gsub("Ctrl_","",mid1$sample)
mid1$replicate=gsub("X\\d+_","",mid1$replicate)
mid1$replicate=as.numeric(mid1$replicate)
str(mid1)
df=data.frame(AA=unique(mid1[,c(1)]),pvalue=100)
str(new)
for (name in unique(mid1$AA)){
  new=filter(mid1,AA==name)#[,c(8,9,4,5)]
  tryCatch({
    #拟合负二项回归模型
    model <- MASS::glm.nb(
      Count ~ timepoint + offset(log(Sum)),
      data = new,
      control = glm.control(maxit = 500)
    )
    # 使用 drop1 进行模型简化并获取 p 值
    drop1_result <- drop1(model, test = "Chisq")
    p_value <- drop1_result["timepoint", "Pr(>Chi)"]
    
    # 将 p 值赋给原始表格的最后一列
    df[df$AA == name, "pvalue"] <- p_value
  },error =function(e){
    #捕获错误并打印警告，但继续执行其他操作
    warning(sprintf("Error processing %s:%s",name, e$message))
  })
}
df$FDR=p.adjust(df$pvalue,method = "fdr")
filter(df,FDR<0.05)


mid=CCAprops %>% group_by(AA,Charge,Sample,sample) %>% dplyr::summarise(Count=sum(count),Sum=sum(sum)) 
mid1=mid %>% group_by(Sample,AA) %>% dplyr::summarise(mean=mean(Count/Sum),se=se(Count/Sum))
# 创建初始向量
vec <- unique(mid1$Sample)
# 组合所有可能的向量
comb <- combn(vec, 2)
# 将所有向量形成一个列表
my_pair <- lapply(seq_len(ncol(comb)), function(i) comb[,i])
ggplot(mid1,aes(x=Sample,y=mean,col=Sample)) +#,grepl("\\w",external_gene_name)
  #geom_violin(mid2,mapping=aes(x=class,y=log10(mean),fill = class))+
  geom_boxplot(width=0.5,col="black")+
  geom_point()+
  #facet_wrap(~type,scales = "free")+
  #my_theme+
  theme_classic()+
  scale_y_continuous(breaks = c(seq(0.4,0.9,0.1)),limits = c(0.4,0.9)) +
  geom_line(aes(group=AA), linetype="dashed", col="grey",size=0.05)+
  geom_line(filter(mid1, AA %in% filter(df,FDR<0.05)$AA),mapping=aes(group=AA), linetype="dashed", col="red",size=0.3)+
  labs(y="AA Charge rate")+
  ggrepel::geom_text_repel(filter(mid1,Sample=="Ctrl", AA %in% filter(df,FDR<0.05)$AA),mapping=aes(x="Ctrl",y=mean,label = AA),color="#FF4040")+
  #ggpubr::stat_compare_means(mapping = aes(group = Sample),method = 'wilcox.test', comparisons = my_pair,paired = T,label = 'p.signif' ) +
  scale_color_manual(values=c("#B7B5B6", "#B7D0EC", "#89B3D9", "#408FBF"))
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/ChargeAA_point.pdf', width= 8 , height= 5 , units='in')



##################################

#wobble

##################################

wobble=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/wobble.csv")
head(wobble)
Deseq2_tRNA=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/cyto_Anticodon_counts_DESEqNormalized.csv")
Deseq2_tRNA=na.omit(filter(Deseq2_tRNA,grepl("Homo_sapiens",Anticodon),
                           !grepl("Sup-TTA",Anticodon),!grepl("iMet",Anticodon),!grepl("SeC",Anticodon)))  %>%
  separate(Anticodon,into = c('Name','AA','anticodon'),sep = '-')
head(Deseq2_tRNA)
colnames(Deseq2_tRNA)=gsub("_R1.fq.gz.unpaired_uniq","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)=gsub("clean_","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)

Deseq2_tRNA=Deseq2_tRNA[,c(2:11)]
head(Deseq2_tRNA)

head(wobble)
Deseq2_tRNA=merge(wobble,Deseq2_tRNA,by.x="anticodon",by.y="anticodon")
Deseq2_tRNA2=tidyr:: gather(Deseq2_tRNA, Sample, value, Ctrl_1:X4_2)
Deseq2_tRNA2$Sample=sub("X","XS",Deseq2_tRNA2$Sample)

Deseq2_tRNA2$tAI=Deseq2_tRNA2$value*(1-Deseq2_tRNA2$s)
tAI=Deseq2_tRNA2 %>% group_by(Codon,Sample) %>% dplyr::summarise(W=sum(tAI))
write.csv(tAI,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/tRNA_Expression_W.csv",row.names = F,quote = F)

#tAI number of charged tRNA
CCAprops=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/CCAprops.csv")
CCAprops=filter(CCAprops,!grepl("eColi",gene),!grepl("mito",gene)) %>%
  separate(gene,into = c('AA','anti-codon','Isodecoder'),sep = '-') %>% filter(!AA %in% c("tRX","Sup"),
                                                                               !grepl("iMet",AA),!grepl("SeC",AA))
CCAprops$Codon=reverse_complement_vector(CCAprops$`anti-codon`)
head(CCAprops)
unique(CCAprops$Codon)
unique(CCAprops$`anti-codon`)
unique(gsub("-\\d+","",CCAprops$gene))
CCAprops$sample=gsub("/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj/clean_","",CCAprops$sample)
CCAprops$sample=gsub("_R1.fq.gz.unpaired_uniq.bam","",CCAprops$sample)
CCAprops$Charge=CCAprops$end
CCAprops[CCAprops$Charge!="CA",]$Charge="No"
CCAprops[CCAprops$Charge=="CA",]$Charge="Yes"
CCAprops=CCAprops %>% group_by(Codon,AA,`anti-codon`,Charge,sample) %>% dplyr::summarise(count=sum(count))

mid1=filter(CCAprops,Charge=="Yes",
            !grepl("iMet",AA),!grepl("SeC",AA))

mid2=CCAprops %>% group_by(Codon,AA,`anti-codon`,sample) %>% dplyr::summarise(total=sum(count))
mid=merge(mid1,mid2,by=c("Codon","AA","anti-codon","sample"))
mid

unique(Deseq2_tRNA$anticodon)
unique(CCAprops$'anti-codon')

Deseq2_tRNA=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/cyto_Anticodon_counts_DESEqNormalized.csv")
Deseq2_tRNA=na.omit(filter(Deseq2_tRNA,grepl("Homo_sapiens",Anticodon),
                           !grepl("Sup-TTA",Anticodon),!grepl("iMet",Anticodon),!grepl("SeC",Anticodon)))  %>%
  separate(Anticodon,into = c('Name','AA','anticodon'),sep = '-')
head(Deseq2_tRNA)
colnames(Deseq2_tRNA)=gsub("_R1.fq.gz.unpaired_uniq","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)=gsub("clean_","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)

Deseq2_tRNA=Deseq2_tRNA[,c(2:11)]
Deseq2_tRNA2=tidyr:: gather(Deseq2_tRNA, Sample, value, Ctrl_1:X4_2)
head(Deseq2_tRNA2)
head(mid[,c(2:4,6,7)])
Deseq2_tRNA2=merge(Deseq2_tRNA2,mid[,c(2:4,6,7)],by.x=c("Sample","AA","anticodon"),by.y=c("sample","AA","anti-codon"))
Deseq2_tRNA2=merge(data.frame(wobble),data.frame(Deseq2_tRNA2),by.x="anticodon",by.y="anticodon")

head(Deseq2_tRNA2)

Deseq2_tRNA2$tAI=Deseq2_tRNA2$value*(Deseq2_tRNA2$count/Deseq2_tRNA2$total)*(1-Deseq2_tRNA2$s)
tAI=Deseq2_tRNA2 %>% group_by(Codon,AA,Sample) %>% dplyr::summarise(W_charged=sum(tAI))
tAI$Sample=sub("X","XS",tAI$Sample)

total_tAI=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/tRNA_Expression_W.csv")
head(tAI)
head(total_tAI)

tAI=merge(tAI,total_tAI,by=c("Codon","Sample"))
write.csv(tAI,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/tRNA_Deseq2_Expression_W.csv",row.names = F,quote = F)

############################################################

#tAI in the transcriptome

############################################################

tAI=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/tRNA_Deseq2_Expression_W.csv")
colnames(tAI)[2]="Sample"
head(tAI)
tAI$sample=gsub("_\\d+","",tAI$Sample)
tAI$charge=tAI$W_charged/tAI$W


mid <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/CUB_Tr_RPF/AsiteStall/HEK293FT_RSCU.txt")
# 排序：先 AA3 → 再 RSCU 降序（你要求的排序）
mid_sorted <- mid %>%
  arrange(
    AA3,                  # 第一排序：氨基酸
    desc(RSCU)            # 第三排序：RSCU 从大到小
  )
#mid_sorted$order=paste0(mid_sorted$AA3,":",mid_sorted$Codon)
mid_sorted$order=paste0(mid_sorted$Codon)

mid=tAI %>% group_by(sample,Codon,AA) %>% dplyr::summarise(charge_mean=mean(charge),charge_se=se(charge))
str(mid)
#mid$Codon=paste0(mid$AA,":",mid$Codon)
mid$Codon=paste0(mid$Codon)
mid=reshape2::dcast(mid,Codon~sample,value.var = "charge_mean")
rownames(mid)=mid$Codon
mid=mid[,-1]
my.colors <- colorRampPalette(c("grey","white","red"))(50)
plt=pheatmap::pheatmap(t(mid)[,mid_sorted$order], angle_col = 90,               # 你的数值矩阵
                       color = my.colors,cluster_rows = F,cluster_cols = F,  # 自定义渐变
                       breaks = seq(0.4, 0.9, 
                                    length.out = 50))  # 保证对称
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Codon_charge.pdf', plt, width= 18 , height= 3 , units='in')
# 
# mid=tAI %>% group_by(sample,Codon,AA) %>% dplyr::summarise(charge_mean=mean(charge),charge_se=se(charge))
# ggplot(mid,mapping=aes(x=paste0(AA,":",Codon),y=charge_mean,fill= sample)) +
#   geom_bar(position="dodge", stat="identity", colour="black")+
#   #my_theme+
#   scale_y_continuous(breaks = c(seq(0,0.9,0.1)),limits = c(0,0.9)) +
#   labs(x="",y="charged rate")+
#   #scale_alpha_discrete(range = c(1, 0.1)) +
#   scale_fill_manual(values=c("#B7B5B6", "#B7D0EC", "#89B3D9", "#408FBF")) 
# ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/ChargeCodon.pdf', width= 30 , height= 5 , units='in')

head(tAI)
mid1=tAI
head(mid1)
mid1$timepoint=0
mid1[grepl("XS1",mid1$Sample),]$timepoint=1
mid1[grepl("XS2",mid1$Sample),]$timepoint=2
mid1[grepl("XS4",mid1$Sample),]$timepoint=3
#mid1=filter(mid1,grepl("X",mid1$sample))
head(mid1)
name="AAA"
mid1$replicate=gsub("Ctrl_","",mid1$Sample)
mid1$replicate=gsub("XS\\d+_","",mid1$replicate)
mid1$replicate=as.numeric(mid1$replicate)
colnames(mid1)
mid1$Codon=paste0(mid1$AA,":",mid1$Codon)
df=data.frame(Codon=unique(mid1[,c(1)]),pvalue=100)
str(new)
for (name in unique(mid1$Codon)){
  new=filter(mid1,Codon==name)[,c(8,9,4,5)]
  tryCatch({
    #拟合负二项回归模型
    model <- MASS::glm.nb(
      W_charged ~ timepoint + offset(log(W)),
      data = new,
      control = glm.control(maxit = 500)
    )
    # 使用 drop1 进行模型简化并获取 p 值
    drop1_result <- drop1(model, test = "Chisq")
    p_value <- drop1_result["timepoint", "Pr(>Chi)"]
    
    # 将 p 值赋给原始表格的最后一列
    df[df$Codon == name, "pvalue"] <- p_value
  },error =function(e){
    #捕获错误并打印警告，但继续执行其他操作
    warning(sprintf("Error processing %s:%s",name, e$message))
  })
}
df$FDR=p.adjust(df$pvalue,method = "fdr")
df=filter(df,FDR<0.05)[order(filter(df,FDR<0.05)$Codon,decreasing = F),]


mid=tAI %>% group_by(sample,Codon,AA) %>% dplyr::summarise(charge_mean=mean(charge),charge_se=se(charge))
# 创建初始向量
vec <- unique(mid$sample)
# 组合所有可能的向量
comb <- combn(vec, 2)
# 将所有向量形成一个列表
my_pair <- lapply(seq_len(ncol(comb)), function(i) comb[,i])
mid$AA=paste0(mid$AA,":",mid$Codon)
ggplot(mid,aes(x=sample,y=charge_mean,col=sample)) +#,grepl("\\w",external_gene_name)
  #geom_violin(mid2,mapping=aes(x=class,y=log10(mean),fill = class))+
  geom_boxplot(width=0.5,col="black")+
  geom_point()+
  #facet_wrap(~type,scales = "free")+
  #my_theme+
  theme_classic()+
  scale_y_continuous(breaks = c(seq(0.4,0.9,0.1)),limits = c(0.4,0.9)) +
  geom_line(mid,mapping=aes(group=AA), linetype="dashed", col="grey",size=0.02)+
  geom_line(filter(mid, AA %in% filter(df,FDR<0.05)$Codon),mapping=aes(group=AA), linetype="dashed", col="red",size=0.3)+
  labs(y="codon Charge rate")+
  ggrepel::geom_text_repel(filter(mid,sample=="Ctrl", AA %in% filter(df,FDR<0.05)$Codon)[seq(1,19,2),],mapping=aes(x="Ctrl",y=charge_mean,label = AA),color="#FF4040")+
  ggrepel::geom_text_repel(filter(mid,sample=="XS4", AA %in% filter(df,FDR<0.05)$Codon)[seq(2,18,2),],mapping=aes(x="XS4",y=charge_mean,label = AA),color="#FF4040")+
  scale_color_manual(values=c("#B7B5B6", "#B7D0EC", "#89B3D9", "#408FBF"))#+
#ggpubr::stat_compare_means(mapping = aes(group = Sample),method = 'wilcox.test', comparisons = my_pair,paired = T,label = 'p.signif' )
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Charge_codon_point.pdf', width= 8 , height= 5 , units='in')

unique(filter(mid, AA %in% filter(df,FDR<0.05)$Codon)$AA)


#########################

#生物学重复的相关性

#########################

tRNA=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/cyto_Isodecoder_counts_DESEqNormalized.csv")
head(tRNA)
tRNA=filter(tRNA,!grepl("coli",isodecoder))
# 计算每行的平均值，忽略 NA
tRNA$row_means <- rowMeans(tRNA[,-c("isodecoder","size")], na.rm = TRUE)
tRNA=filter(tRNA,row_means>=100)

cor.test(log10(tRNA$clean_Ctrl_1_R1.fq.gz.unpaired_uniq),log10(tRNA$clean_Ctrl_2_R1.fq.gz.unpaired_uniq),method = "p")
ggplot(tRNA,aes(x=log10(clean_Ctrl_1_R1.fq.gz.unpaired_uniq),y=log10(clean_Ctrl_2_R1.fq.gz.unpaired_uniq))) +
  geom_point()+
  my_theme+
  geom_abline(slope = 1,intercept = 0,col="red")+
  scale_x_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  scale_y_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",`anti-codon`)),color="blue")+#paste0(AA,":",Codon)
  ggpubr::stat_cor(method = "p", col="red")+
  labs(x="tRNA expression in Ctrl replication 1",y="tRNA expression in Ctrl replication 2")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Ctrl_rep_tRNA.pdf', width= 5 , height= 5 , units='in')

cor.test(log10(tRNA$clean_X1_1_R1.fq.gz.unpaired_uniq),log10(tRNA$clean_X1_2_R1.fq.gz.unpaired_uniq),method = "p")
ggplot(tRNA,aes(x=log10(clean_X1_1_R1.fq.gz.unpaired_uniq),y=log10(clean_X1_2_R1.fq.gz.unpaired_uniq))) +
  geom_point()+
  my_theme+
  geom_abline(slope = 1,intercept = 0,col="red")+
  scale_x_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  scale_y_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",`anti-codon`)),color="blue")+#paste0(AA,":",Codon)
  ggpubr::stat_cor(method = "p", col="red")+
  labs(x="tRNA expression in XS1 replication 1",y="tRNA expression in XS1 replication 2")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/XS1_rep_tRNA.pdf', width= 5 , height= 5 , units='in')

#cor.test(log10(tRNA$clean_X4_1_R1.fq.gz.unpaired_uniq),log10(tRNA$clean_X4_2_R1.fq.gz.unpaired_uniq),method = "s")
cor.test(log10(tRNA$clean_X2_1_R1.fq.gz.unpaired_uniq),log10(tRNA$clean_X2_2_R1.fq.gz.unpaired_uniq),method = "p")
ggplot(tRNA,aes(x=log10(clean_X2_1_R1.fq.gz.unpaired_uniq),y=log10(clean_X2_2_R1.fq.gz.unpaired_uniq))) +
  geom_point()+
  my_theme+
  geom_abline(slope = 1,intercept = 0,col="red")+
  scale_x_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  scale_y_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",`anti-codon`)),color="blue")+#paste0(AA,":",Codon)
  ggpubr::stat_cor(method = "p", col="red")+
  labs(x="tRNA expression in XS2 replication 1",y="tRNA expression in XS2 replication 2")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/XS2_rep_tRNA.pdf', width= 5 , height= 5 , units='in')

cor.test(log10(tRNA$clean_X4_1_R1.fq.gz.unpaired_uniq),log10(tRNA$clean_X4_2_R1.fq.gz.unpaired_uniq),method = "p")
ggplot(tRNA,aes(x=log10(clean_X4_1_R1.fq.gz.unpaired_uniq),y=log10(clean_X4_2_R1.fq.gz.unpaired_uniq))) +
  geom_point()+
  my_theme+
  geom_abline(slope = 1,intercept = 0,col="red")+
  scale_x_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  scale_y_continuous(breaks = c(seq(1,7,1)),limits = c(1,7)) +
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",`anti-codon`)),color="blue")+#paste0(AA,":",Codon)
  ggpubr::stat_cor(method = "p", col="red")+
  labs(x="tRNA expression in XS4 replication 1",y="tRNA expression in XS4 replication 2")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/XS4_rep_tRNA.pdf', width= 5 , height= 5 , units='in')

head(tRNA[,2:3])
head(rowMeans(tRNA[,2:3], na.rm = TRUE))

ggplot(tRNA,aes(x=log10(clean_X4_1_R1.fq.gz.unpaired_uniq))) +
  geom_density()+
  my_theme

#########################

#tRNA abundance

#########################

tRNA_num=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/X4vsCtrl_diffexpr-results.csv")
#tRNA_num=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/X1vsX4_diffexpr-results.csv")
tRNA_num$Gene=gsub("Homo_sapiens_tRNA-","",tRNA_num$Gene)
tRNA_num=na.omit(filter(tRNA_num,!grepl("eColi",Gene),!grepl("Sup-TTA",Gene),!grepl("Ser-ACT",Gene),baseMean>100))
colnames(tRNA_num)=gsub("_R1.fq.gz.unpaired_uniq","",colnames(tRNA_num))
colnames(tRNA_num)=gsub("clean_X","clean_XS",colnames(tRNA_num))
colnames(tRNA_num)=gsub("clean_","",colnames(tRNA_num))
head(tRNA_num)

 tRNA_num=tRNA_num[,-c(1,3:8,17)]
# wilcox.test(log(tRNA_num$XS1_1+tRNA_num$XS1_2),log(tRNA_num$XS2_1+tRNA_num$XS2_2),paired = T)
# wilcox.test((tRNA_num$XS1_1+tRNA_num$XS1_2),(tRNA_num$XS2_1+tRNA_num$XS2_2),paired = T)
# 
 tRNA_num=gather(tRNA_num,sample,Exp,Ctrl_1:XS4_2)
 tRNA_num$Sample=gsub("_\\d","",tRNA_num$sample)
 tRNA_num2=tRNA_num %>% group_by(Gene,Sample) %>% dplyr::summarise(mean=mean(Exp),se=se(Exp))
 head(tRNA_num2)

 # 1. 拆分 AA + anticodon，计算平均丰度，排序
 tRNA_sorted <- tRNA_num2 %>%
   # 拆分 Gene 列为 AA 和 anticodon
   separate(Gene, into = c("AA", "anticodon"), sep = "-", remove = FALSE) %>%
   
   # 计算4个样本的平均丰度
   mutate(mean_abun = rowMeans(across(c(Ctrl, XS1, XS2, XS4)), na.rm = TRUE)) %>%
   
   # 先按 AA 排序，AA 内部再按 mean_abun 从大到小排序
   arrange(AA, desc(mean_abun))
 
 # 查看最终结果（完美排序）
 tRNA_sorted %>% select(AA, anticodon, Gene, mean_abun, Ctrl, XS1, XS2, XS4)
 gene_order_sorted <- tRNA_sorted$Gene
 gene_order_sorted
 write.csv(
   data.frame(Gene = gene_order_sorted),  # 转成数据框
   file = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/gene_order_sorted.csv",        # 文件名
   row.names = FALSE                      # 不导出行号
 )
# # 创建初始向量
# vec <- unique(tRNA_num2$Sample)
# # 组合所有可能的向量
# comb <- combn(vec, 2)
# # 将所有向量形成一个列表
# my_pair <- lapply(seq_len(ncol(comb)), function(i) comb[,i])
# ggplot(tRNA_num2,aes(x=Gene,y=log10(mean),fill=Sample)) +#,grepl("\\w",external_gene_name)
#   geom_bar(stat = "identity", color = "black", position = position_dodge()) +
#   #geom_point()+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1,color="black"))+
#   labs(y="tRNA abundance (log10)")+
#   scale_fill_manual(values=c("#B7B5B6", "#B7D0EC", "#89B3D9", "#408FBF"))
# ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/tRNA_abundace.pdf', width= 8 , height= 3 , units='in')

 tRNA_num2=reshape2::dcast(tRNA_num2,Gene~Sample,value.var = "mean")
 rownames(tRNA_num2)=tRNA_num2$Gene
 
 my.colors <- colorRampPalette(c("grey","white", "red"))(50)
 plt=pheatmap::pheatmap(log10(tRNA_num2[gene_order_sorted,-1]),                # 你的数值矩阵
                    color = my.colors,cluster_rows = F,cluster_cols = F,  # 自定义渐变
                    breaks = seq(min(log10(tRNA_num2[,-1])), max(log10(tRNA_num2[,-1])), 
                    length.out = 50))  # 保证对称
 ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/tRNA_abundace.pdf', plt, width= 3 , height= 8 , units='in')

 my.colors <- colorRampPalette(c("blue", "white", "red"))(11)
 
 tRNA_num=tRNA_num2[,-1] 
 tRNA_num=log10(tRNA_num/tRNA_num$Ctrl)
 min(tRNA_num)
 max(tRNA_num)
 plt=pheatmap::pheatmap(tRNA_num[,-1],                # 你的数值矩阵
                    color = my.colors,
                    cluster_rows = F,cluster_cols = F,  # 自定义渐变
                    breaks = seq(-0.5, 0.5,0.1), 
                    length.out = 11)  # 保证对称
 plt
 ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/tRNA_abundace.pdf', plt, width= 3 , height= 8 , units='in')

 
 ############################################################
 
 #tRNA modification
 
 ############################################################
 # 自定义颜色梯度：
 # - 从-0.2（深蓝）到0（白色）：生成100个蓝色渐变
 # - 从0（白色）到0.6（深红）：生成100个红色渐变
 blue_values <- seq(-0.6, 0, length.out = 100)
 red_values <- seq(0, 0.6, length.out = 100)
 
 # 生成颜色向量（蓝色渐变 + 红色渐变，去除重复的白色）
 blue_colors <- colorRampPalette(c("#00008B", "white"))(100)  # 深蓝到白
 red_colors <- colorRampPalette(c("white", "#8B0000"))(100)   # 白到深红
 custom_colors <- c(blue_colors[-100], red_colors)  # 合并，去除重复的白色
 Keep=10
 
 XS1=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/X1vsCtrl_logOR.csv")
 XS1_heatmap=as.matrix(XS1[,-1])
 rownames(XS1_heatmap)=XS1$V1
 
 # 
 pheatmap::pheatmap((XS1_heatmap),
                    cluster_cols = F,
                    cluster_rows = F,
                    color = custom_colors,
                    breaks = c(blue_values, red_values[-1]),  # 颜色分割点（与颜色对应
                    angle_col = 90)
 
 XS1_pos= gather(XS1,Pos,OR,`0`:`76`)
 head(XS1_pos)
 XS1_pos=filter(XS1_pos,OR!=0) %>% group_by(Pos) %>% dplyr::summarise(count=n()) %>% filter(count>=Keep)
 
 colSums(XS1[,-1])
 XS1=XS1[,c("V1","58")]
 nrow(filter(XS1,`58`>0))
 nrow(filter(XS1,`58`<0))
 colnames(XS1)[2]="XS1vsCtrl"
 
 
 XS2=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/X2vsCtrl_logOR.csv")
 XS2_heatmap=as.matrix(XS2[,-1])
 rownames(XS2_heatmap)=XS2$V1
 pheatmap::pheatmap((XS2_heatmap),
                    cluster_cols = F,
                    cluster_rows = F,
                    color = custom_colors,
                    breaks = c(blue_values, red_values[-1]),  # 颜色分割点（与颜色对应
                    angle_col = 90)
 
 head(XS2)
 XS2_pos= gather(XS2,Pos,OR,`0`:`76`)
 head(XS2_pos)
 XS2_pos=filter(XS2_pos,OR!=0) %>% group_by(Pos) %>% dplyr::summarise(count=n()) %>% filter(count>=Keep)
 
 colSums(XS2[,-1])
 nrow(filter(XS2,`58`>0))
 nrow(filter(XS2,`58`<0))
 XS2=XS2[,c("V1","58")]
 colnames(XS2)[2]="XS2vsCtrl"
 
 
 XS4=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/X4vsCtrl_logOR.csv")
 # XS4_heatmap=as.matrix(XS4[,-1])
 # rownames(XS4_heatmap)=XS4$V1
 pheatmap::pheatmap(XS4_heatmap,
                    cluster_cols = F,
                    cluster_rows = F,
                    #color = custom_color,
                    #breaks = c(blue_values, red_values[-1]),  # 颜色分割点（与颜色对应
                    angle_col = 90)
 head(XS4)
 
 nrow(filter(XS4,`58`>0))
 nrow(filter(XS4,`58`<0))
 # rownames(XS4)=XS4$V1
 colSums(XS4[,-1])
 XS4_pos= gather(XS4,Pos,OR,`0`:`76`)
 head(XS4_pos)
 XS4_pos=filter(XS4_pos,OR!=0) %>% group_by(Pos) %>% dplyr::summarise(count=n()) %>% filter(count>=Keep)
 
 
 XS4=XS4[,c("V1","58")]
 colnames(XS4)[2]="XS4vsCtrl"
 colSums(XS4[,-1])
 
 XS1_XS2=merge(XS1,XS2,by="V1")
 XS1_XS2_XS4=merge(XS1_XS2,XS4,by="V1")
 colnames(XS1_XS2_XS4)[1]="Isodecoder"
 XS1_XS2_XS4$Isodecoder=gsub("-multi","",XS1_XS2_XS4$Isodecoder)
 XS1_XS2_XS4=data.frame(XS1_XS2_XS4)
 rownames(XS1_XS2_XS4)=XS1_XS2_XS4$Isodecoder
 
 pheatmap::pheatmap(t(XS1_XS2_XS4[,-1]),
                    cluster_cols = F,
                    cluster_rows = F,
                    color = custom_colors,
                    breaks = c(blue_values, red_values[-1]),  # 颜色分割点（与颜色对应
                    angle_col = 90)
 
 keep3=unique(c(XS1_pos$Pos,XS2_pos$Pos,XS4_pos$Pos))
 sample=c(paste0(keep3,"_XS1"),paste0(keep3,"_XS2"),paste0(keep3,"_XS4"))
 sample=sample[order(sample)]
 {
   XS1=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/X1vsCtrl_logOR.csv")
   colnames(XS1)[1]="isodecoder"
   #colnames(XS1)[-1]=paste0("XS1_",colnames(XS1)[-1])
   colnames(XS1)[-1]=paste0(as.character(colnames(XS1)[-1]),"_XS1")
   
   XS2=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/X2vsCtrl_logOR.csv")
   colnames(XS2)[1]="isodecoder"
   #colnames(XS2)[-1]=paste0("XS2_",colnames(XS2)[-1])
   colnames(XS2)[-1]=paste0(as.character(colnames(XS2)[-1]),"_XS2")
   #XS2$Sample="XS2"
   
   XS4=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/X4vsCtrl_logOR.csv")
   colnames(XS4)[1]="isodecoder"
   #colnames(XS4)[-1]=paste0("XS4_",colnames(XS4)[-1])
   colnames(XS4)[-1]=paste0(as.character(colnames(XS4)[-1]),"_XS4")
   #XS4$Sample="XS4"
   
   head(XS4)
   XS1_XS2_XS4=merge(XS1,XS2,by=c("isodecoder"))
   XS1_XS2_XS4=merge(XS1_XS2_XS4,XS4,by=c("isodecoder"))
   XS1_XS2_XS4$isodecoder=gsub("-multi","",XS1_XS2_XS4$isodecoder)
   (XS1_XS2_XS4[1:4,1:4])
   #rownames(XS1_XS2_XS4)=XS1_XS2_XS4$Isodecoder
   
   plt=as.data.frame(XS1_XS2_XS4[,-1])
   rownames(plt)=XS1_XS2_XS4$isodecoder
   # (plt[1:4,1:4])
   # pheatmap::pheatmap(t(plt[,c("9_XS1","9_XS2","9_XS4",
   #                             "20_XS1","20_XS2","20_XS4",
   #                             "32_XS1","32_XS2","32_XS4",
   #                             "37_XS1","37_XS2","37_XS4",
   #                             "58_XS1","58_XS2","58_XS4")]),
   #                    cluster_cols = F,
   #                    cluster_rows = F,
   #                    color = custom_colors,
   #                    breaks = c(blue_values, red_values[-1]),  # 颜色分割点（与颜色对应
   #                    angle_col = 90#,
   #                    #filename = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/tRNA_mod.pdf"
   #                    )
   
   pheatmap::pheatmap(t(plt[,sample]),
                      cluster_cols = F,
                      cluster_rows = F,
                      color = custom_colors,
                      breaks = c(blue_values, red_values[-1]),  # 颜色分割点（与颜色对应
                      angle_col = 90#,
                      #filename = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/tRNA_mod.pdf"
   ) 
   pheatmap::pheatmap(plt,
                      cluster_cols = F,
                      cluster_rows = F,
                      color = custom_colors,
                      breaks = c(blue_values, red_values[-1]),  # 颜色分割点（与颜色对应
                      angle_col = 90)
   
 }
 min(plt)
 max(plt)
 
 ###############
 
 #OR
 
 ##################################
mid=tAI %>% group_by(sample,Codon,AA) %>% dplyr::summarise(W_charge_mean=mean(W_charged),W_mean=mean(W))
filter(mid,Codon=="TTA")

OR=data.frame()
for (i in unique(mid$Codon)) {
 dat=filter(mid,Codon==i)
# 1. 把 Ctrl 单独拿出来
ctrl_tbl <- dat %>% 
  #ungroup() %>%  
  filter(sample == "Ctrl") %>% 
  select(Codon,AA, Ctrl_W = W_charge_mean, Ctrl_Total = W_mean)
colnames(ctrl_tbl)[1]="Sample"
colnames(ctrl_tbl)[3]="AA1"

# 2. 非 Ctrl 样本
test_tbl <- dat %>% 
  filter(sample != "Ctrl")

# 3. 合并、建 2×2 表、跑 Fisher
res <- test_tbl %>% 
  left_join(ctrl_tbl, by = "Codon") %>% 
  mutate(
    a = W_charge_mean,                  # 实验组“成功”
    b = W_mean - W_charge_mean,         # 实验组“失败”
    c = Ctrl_W,                         # Ctrl“成功”
    d = Ctrl_Total - Ctrl_W             # Ctrl“失败”
  ) %>% 
  rowwise() %>% 
  summarise(
    sample = sample,
    Codon  = Codon,
    AA  = AA,
    OR     = fisher.test(matrix(round(c(a,b,c,d)), nrow = 2))$estimate,
    p      = fisher.test(matrix(round(c(a,b,c,d)), nrow = 2))$p.value,
    .groups = "drop"
  ) 
OR=rbind(res,OR)
}
OR

library(dplyr)

OR <- OR %>%
  group_by(sample) %>%                # 每个样本内部校正
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    sig   = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  ) %>%
  ungroup()

unique(OR$sample)
OR$Codon=paste0(OR$AA,":",OR$Codon)
OR=OR[order(OR$Codon,decreasing = F),]
# 2. 颜色：蓝-白-红，以 1 为中心
my_col <- colorRampPalette(c("steelblue","white","firebrick"))(100)

mat_or  <- OR %>%
  select(Codon, sample, OR) %>%        # 只拿三列
  pivot_wider(names_from = sample, values_from = OR) %>%
  column_to_rownames("Codon") %>%      # 行名
  as.matrix()

mat_sig <-OR %>%
  select(Codon, sample, sig) %>%
  pivot_wider(names_from = sample, values_from = sig) %>%
  column_to_rownames("Codon") %>%
  as.matrix()

pheatmap::pheatmap(
  mat_or,
  color           = my_col,
  breaks          = seq(min(mat_or), max(mat_or), length.out = 101),
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = mat_sig,
  number_color    = "black",
  fontsize_number = 8,
  legend_labels   = c("low OR", "high OR"),
  main            = "OR vs Ctrl (Fisher test)"
)
 



##################################

#wobble

##################################

wobble=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/wobble.csv")
head(wobble)
Deseq2_tRNA=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/cyto_Isodecoder_counts_DESEqNormalized.csv")
Deseq2_tRNA=na.omit(filter(Deseq2_tRNA,grepl("Homo_sapiens", isodecoder),
                           !grepl("Sup-TTA", isodecoder),!grepl("iMet", isodecoder),!grepl("SeC", isodecoder)))  %>%
  separate( isodecoder,into = c('Name','AA','isodecoder','ID'),sep = '-')
head(Deseq2_tRNA)
colnames(Deseq2_tRNA)=gsub("_R1.fq.gz.unpaired_uniq","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)=gsub("clean_","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)

Deseq2_tRNA=Deseq2_tRNA[,c(2:12)]
head(Deseq2_tRNA)

head(wobble)
Deseq2_tRNA=merge(wobble,Deseq2_tRNA,by.x="anticodon",by.y="isodecoder", allow.cartesian = T)

Deseq2_tRNA2=tidyr:: gather(Deseq2_tRNA, Sample, value, Ctrl_1:X4_2)
Deseq2_tRNA2$Sample=sub("X","XS",Deseq2_tRNA2$Sample)

Deseq2_tRNA2$tAI=Deseq2_tRNA2$value*(1-Deseq2_tRNA2$s)
tAI=Deseq2_tRNA2 %>% group_by(anticodon,ID, Codon,Sample) %>% dplyr::summarise(W=sum(tAI))
write.csv(tAI,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/tRNA_Expression_W.csv",row.names = F,quote = F)

#tAI number of charged tRNA
CCAprops=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/CCAprops.csv")
unique(CCAprops$gene)
CCAprops=filter(CCAprops,!grepl("eColi",gene),!grepl("mito",gene)) %>%
  separate(gene,into = c('AA','Isodecoder','ID'),sep = '-') %>% filter(!AA %in% c("tRX","Sup"),
                                                                               !grepl("iMet",AA),!grepl("SeC",AA))
CCAprops$Codon=reverse_complement_vector(CCAprops$`Isodecoder`)
head(CCAprops)
unique(CCAprops$Codon)
unique(CCAprops$`Isodecoder`)
unique(gsub("-\\d+","",CCAprops$gene))
CCAprops$sample=gsub("/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj/clean_","",CCAprops$sample)
CCAprops$sample=gsub("_R1.fq.gz.unpaired_uniq.bam","",CCAprops$sample)
CCAprops$Charge=CCAprops$end
CCAprops[CCAprops$Charge!="CA",]$Charge="No"
CCAprops[CCAprops$Charge=="CA",]$Charge="Yes"
CCAprops=CCAprops %>% group_by(Codon,AA,Isodecoder,ID,Charge,sample) %>% dplyr::summarise(count=sum(count))

mid1=filter(CCAprops,Charge=="Yes",
            !grepl("iMet",AA),!grepl("SeC",AA))

mid2=CCAprops %>% group_by(Codon,AA,Isodecoder,ID,sample) %>% dplyr::summarise(total=sum(count))
mid=merge(mid1,mid2,by=c("Codon","AA","Isodecoder","ID","sample"))
mid


Deseq2_tRNA=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/cyto_Isodecoder_counts_DESEqNormalized.csv")
unique(Deseq2_tRNA$isodecoder)
Deseq2_tRNA=na.omit(filter(Deseq2_tRNA,grepl("Homo_sapiens", isodecoder),
                           !grepl("Sup-TTA", isodecoder)))  %>%
  separate( isodecoder,into = c('Name','AA','isodecoder','ID'),sep = '-')
head(Deseq2_tRNA)
colnames(Deseq2_tRNA)=gsub("_R1.fq.gz.unpaired_uniq","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)=gsub("clean_","",colnames(Deseq2_tRNA))
colnames(Deseq2_tRNA)

Deseq2_tRNA=Deseq2_tRNA[,c(2:12)]
head(Deseq2_tRNA)
Deseq2_tRNA2=tidyr:: gather(Deseq2_tRNA, Sample, value, Ctrl_1:X4_2)
head(Deseq2_tRNA2)
Deseq2_tRNA2=merge(data.frame(wobble),data.frame(Deseq2_tRNA2),by.x="anticodon",by.y="isodecoder")
head(Deseq2_tRNA2)

Deseq2_tRNA2$tAI=Deseq2_tRNA2$value*(1-Deseq2_tRNA2$s)
tAI=Deseq2_tRNA2 %>% group_by(Codon,anticodon,ID,AA,Sample) %>% dplyr::summarise(W=sum(tAI))
tAI$isodecoder=paste(tAI$AA,tAI$anticodon,tAI$ID,sep = "-")
unique(tAI$isodecoder)
head(tAI)
tAI$isodecoder <- sub("/.*", "", tAI$isodecoder)

test=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/20250709/modPos_totalMisincProp_skj.csv")
test=filter(test,canon_pos==58,!grepl("eColi",isodecoder),!grepl("mito",isodecoder))[,c(1,3,6)]
test$Sample=gsub("/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj/clean_","",test$bam)
test$Sample=gsub("_R1.fq.gz.unpaired_uniq.bam","",test$Sample)
test$isodecoder=gsub("-multi","",test$isodecoder)
test$isodecoder=gsub("tRX-","",test$isodecoder)
test$isodecoder <- sub("/.*", "", test$isodecoder)
unique(test$isodecoder)
head(test)

test2=merge(tAI,test,by=c("isodecoder","Sample"),all = T)
head(test2)
# 找出所有带 NA 的行
na_rows=test2[!complete.cases(test2), ] %>% data.frame()
test2[is.na(test2)]=0


mid=test2 %>% group_by(AA,Codon,Sample) %>% dplyr::summarise(W_total=sum(W),W_m1A59=sum(W*misinc_proportion))
mid=filter(mid,Codon!=0,Sample!=0,AA!=0)
mid$condition=gsub("_1","",mid$Sample)
mid$condition=gsub("_2","",mid$condition)


mid=mid %>% group_by(AA,Codon,condition) %>% dplyr::summarise(m1A59_frac=mean(W_m1A59/W_total))
mid=filter(mid,!AA %in% c("iMet"))

Final=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/tRNA/Codon_Final.csv")
head(Final)
Final$condition=sub("XS","X",Final$condition)
mid=merge(mid,Final[,c("Codon","RSCU", "condition")],by=c("Codon","condition"))
plot_list <- list()
temp=unique(mid$AA)
i=1
for(i in 1:length(temp)){
  
  # 当前 AA
  current_AA <- temp[i]
  
  # 筛选当前 AA 数据
  subdata <- mid %>% filter(AA == current_AA)
  
  # 按 Ctrl 的 m1A59_frac 从大到小排序 Codon
  ctrl_sub <- subdata %>% filter(condition == "Ctrl") %>% arrange(desc(RSCU))#arrange(desc(m1A59_frac))
  codon_order <- ctrl_sub$Codon
  subdata$Codon <- factor(subdata$Codon, levels = codon_order)
  
  # 画图
  plot_list[[i]] <- ggplot(subdata, aes(x = Codon, y = m1A59_frac, color = condition, group = condition)) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = c("#575959","#E08251","#9D71A6","#4088B4")) + 
    labs(y = "m1A59 level", x = "Codon", title = current_AA) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# 4. 把所有图拼在一起输出
gridExtra::grid.arrange(grobs = plot_list)
library(dplyr)
# 按 condition 分组，分别做 Spearman 相关
cor_result <- data.frame()

for(cond in unique(mid$condition)){
  
  # 筛选当前组
  sub <- mid[mid$condition == cond, ]
  
  # 计算 Spearman 相关
  test <- cor.test(sub$RSCU, sub$m1A59_frac, method = "spearman")
  test2 <- cor.test(sub$RSCU, sub$m1A59_frac, method = "p")
  
  # 保存结果
  cor_result <- rbind(cor_result, data.frame(
    condition = cond,
    n=nrow(sub),
    rho = test$estimate,
    p_spearman = test$p.value,
    R = test2$estimate,
    p_pearson = test2$p.value
  ))
}

# 查看结果
cor_result

cor.test(mid$m1A59_frac,mid$RSCU,method = "s")
cor.test(mid$m1A59_frac,mid$RSCU,method = "p")
