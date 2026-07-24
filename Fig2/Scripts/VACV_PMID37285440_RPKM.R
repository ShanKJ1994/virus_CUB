##############################################

#Sum

##############################################

Our=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/featureCount/featureCounts_VACV_HeLa229_SciAdv.txt")[,-c(2:5)]
Our$Geneid=gsub("gene:","",Our$Geneid)
head(Our)

mRNA=colnames(Our)[grepl("mRNA",colnames(Our))]
Ribo=colnames(Our)[grepl("Ribo",colnames(Our))]
TE=gsub("Ribo","TE",Ribo)

#Genes with >= 5 reads
virus=Our[!grepl("ENS",Our$Geneid),]
Host=Our[grepl("ENS",Our$Geneid),]
min_values <- apply(data.frame(Host)[, -c(1:2)], 1, function(row) min(row))
condition <- min_values >= 0
Host <- Host[condition, ]
Our=rbind(Host,virus)

#RPKM
norm=colSums(Our[, -c(1:2)])
# norm=data.frame(norm)
# norm$sample=rownames(norm)
# norm$Type=gsub("_\\w+","",norm$sample)
# norm[grepl("Ribo",norm$sample),]$Type="Ribo"
# norm[grepl("mRNA",norm$sample),]$Type="mRNA"
# norm=norm %>% group_by(sample,Type) %>% dplyr::summarise(count=sum(norm))
# min_norm=norm %>% group_by(Type) %>% dplyr::summarise(count_min=min(count))
# norm=merge(norm,min_norm,by=c("Type"))
# norm$size=round(norm$count/norm$count_min,4)
# norm

rpkm_values <- sweep(Our[, -c(1:2)], 1, Our$Length, "/")  # 每个基因的Counts除以基因长度（kb）
rpkm_values <- sweep(rpkm_values, 2,norm, "/")  # 然后除以样本的总Read计数（百万）
rpkm_values <- rpkm_values*1000000000
rpkm_values$Geneid=Our$Geneid
rpkm_values$Length=Our$Length
head(rpkm_values)
mRNA_RPKM=rpkm_values
head(mRNA_RPKM)
virus=mRNA_RPKM[!grepl("ENS",mRNA_RPKM$Geneid),]
head(virus)

#Genes with mRNA > 0 and Ribo > 0
min_values <- apply(data.frame(mRNA_RPKM)[, mRNA], 1, function(row) min(row))
condition <- min_values > 0
mRNA_RPKM <- mRNA_RPKM[condition, ]
min_values <- apply(data.frame(mRNA_RPKM)[, Ribo], 1, function(row) min(row))
condition <- min_values > 0
Host <- mRNA_RPKM[condition, ]
Host=Host[grepl("ENS",Host$Geneid),]

rpkm_values=rbind(Host,virus)
head(rpkm_values)

TE_ratio=rpkm_values[,Ribo]/rpkm_values[,mRNA]
colnames(TE_ratio)=gsub("Ribo","TE",colnames(TE_ratio))
head(TE_ratio)

rpkm_values=cbind(TE_ratio,rpkm_values)
rpkm_values=rpkm_values[,c("Geneid","Length",mRNA,Ribo,TE)]
write.csv(rpkm_values,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/featureCount/VACV_HeLa229_normSum.csv",row.names = F,quote = F)

tail(rpkm_values)

##############################################

#Deseq2

##############################################

Our=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/featureCount/featureCounts_VACV_HeLa229_SciAdv.txt")[,-c(2:5)]
Our$Geneid=gsub("gene:","",Our$Geneid)
head(Our)

mRNA=colnames(Our)[grepl("mRNA",colnames(Our))]
Ribo=colnames(Our)[grepl("Ribo",colnames(Our))]
TE=gsub("Ribo","TE",Ribo)

#Genes with >= 5 reads
virus=Our[!grepl("ENS",Our$Geneid),]
Host=Our[grepl("ENS",Our$Geneid),]
min_values <- apply(data.frame(Host)[, -c(1:2)], 1, function(row) min(row))
condition <- min_values >= 0
Host <- Host[condition, ]
Our=rbind(Host,virus)


#mRNA
FP=as.matrix(data.frame(Our)[grepl("ENS",Our$Geneid), mRNA])
nrow(FP)
head(FP)

condition=as.matrix(data.frame(condition=factor(mRNA)))
rownames(condition)<-colnames(FP)

library("DESeq2")
#使用count matrix(cts)和样本信息(coldata)表格, 构建DESeqDataSet(dds)对象；
dds <- DESeqDataSetFromMatrix(countData = floor(FP),
                              colData = condition,
                              design = ~ condition)
# 计算大小因子
dds <- estimateSizeFactors(dds)
dds$sizeFactor

norm=dds$sizeFactor
Our=data.frame(Our)
rpkm_values <- sweep(Our[, mRNA], 1, Our$Length, "/")  # 每个基因的Counts除以基因长度（kb）
rpkm_values <- sweep(rpkm_values, 2,norm, "/")  # 然后除以样本的总Read计数（百万）
rpkm_values$Geneid=Our$Geneid
rpkm_values$Length=Our$Length
head(rpkm_values)
mRNA_RPKM=rpkm_values
head(mRNA_RPKM)


#Ribo
FP=as.matrix(data.frame(Our)[grepl("ENS",Our$Geneid), Ribo])
nrow(FP)
head(FP)

condition=as.matrix(data.frame(condition=factor(Ribo)))
rownames(condition)<-colnames(FP)

#使用count matrix(cts)和样本信息(coldata)表格, 构建DESeqDataSet(dds)对象；
dds <- DESeqDataSetFromMatrix(countData = floor(FP),
                              colData = condition,
                              design = ~ condition)
# 计算大小因子
dds <- estimateSizeFactors(dds)
dds$sizeFactor

norm=dds$sizeFactor
Our=data.frame(Our)
rpkm_values <- sweep(Our[, Ribo], 1, Our$Length, "/")  # 每个基因的Counts除以基因长度（kb）
rpkm_values <- sweep(rpkm_values, 2,norm, "/")  # 然后除以样本的总Read计数（百万）
rpkm_values$Geneid=Our$Geneid
rpkm_values$Length=Our$Length
head(rpkm_values)
Ribo_RPKM=rpkm_values
head(Ribo_RPKM)

rpkm_values=merge(mRNA_RPKM,Ribo_RPKM,by=c("Geneid","Length"))
TE_ratio=rpkm_values[,Ribo]/rpkm_values[,mRNA]
colnames(TE_ratio)=gsub("Ribo","TE",colnames(TE_ratio))
head(TE_ratio)

rpkm_values=cbind(TE_ratio,rpkm_values)
rpkm_values=rpkm_values[,c("Geneid","Length",mRNA,Ribo,TE)]

write.csv(rpkm_values,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/featureCount/VACV_HeLa229_Deseq2.csv",row.names = F,quote = F)

##############################################

#Host Sum

##############################################

Our=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/featureCount/featureCounts_VACV_HeLa229_SciAdv.txt")[,-c(2:5)]
Our$Geneid=gsub("gene:","",Our$Geneid)
Our=Our[grepl("ENS",Our$Geneid),]
head(Our)

mRNA=colnames(Our)[grepl("mRNA",colnames(Our))]
Ribo=colnames(Our)[grepl("Ribo",colnames(Our))]
TE=gsub("Ribo","TE",Ribo)

#Genes with >= 5 reads
Host=Our[grepl("ENS",Our$Geneid),]
min_values <- apply(data.frame(Host)[, -c(1:2)], 1, function(row) min(row))
condition <- min_values >= 0
Host <- Host[condition, ]
Our=Host

#RPKM
norm=colSums(Our[, -c(1:2)])
rpkm_values <- sweep(Our[, -c(1:2)], 1, Our$Length, "/")  # 每个基因的Counts除以基因长度（kb）
rpkm_values <- sweep(rpkm_values, 2,norm, "/")  # 然后除以样本的总Read计数（百万）
rpkm_values <- rpkm_values*1000000000
rpkm_values$Geneid=Our$Geneid
rpkm_values$Length=Our$Length
head(rpkm_values)
mRNA_RPKM=rpkm_values
head(mRNA_RPKM)

#Genes with mRNA > 0 and Ribo > 0
min_values <- apply(data.frame(mRNA_RPKM)[, mRNA], 1, function(row) min(row))
condition <- min_values > 0
mRNA_RPKM <- mRNA_RPKM[condition, ]
min_values <- apply(data.frame(mRNA_RPKM)[, Ribo], 1, function(row) min(row))
condition <- min_values > 0
Host <- mRNA_RPKM[condition, ]
Host=Host[grepl("ENS",Host$Geneid),]

rpkm_values=Host
head(rpkm_values)

TE_ratio=rpkm_values[,Ribo]/rpkm_values[,mRNA]
colnames(TE_ratio)=gsub("Ribo","TE",colnames(TE_ratio))
head(TE_ratio)

rpkm_values=cbind(TE_ratio,rpkm_values)
rpkm_values=rpkm_values[,c("Geneid","Length",mRNA,Ribo,TE)]
write.csv(rpkm_values,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/featureCount/VACV_HeLa229_HostnormSum.csv",row.names = F,quote = F)

tail(rpkm_values)
