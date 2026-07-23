temp <- list.files(path = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/Frame012/",pattern = "*.count")

i=1
Host_Virus=data.frame()
for(i in c(3:5,31:33,37:38,40:43)){
name=temp[i]
name=gsub("\\.count","",name)
file=read.table(paste0("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/Frame012/",temp[i]),sep='\t',header = T)

file_Tr_mean=file %>% group_by(Tr_ID) %>% dplyr::summarise(TrMean=mean(Ribo/mRNA))
file=merge(file,file_Tr_mean,by=c("Tr_ID"))
Psite=filter(file,!grepl("ENS",Tr_ID)) %>% group_by(Psite) %>% dplyr::summarise(Psite_Mean=mean((Ribo/mRNA)/TrMean),Psite_Median=median((Ribo/mRNA)/TrMean),Virus_Psite_Codon_num=n())
Asite=filter(file,!grepl("ENS",Tr_ID)) %>% group_by(Asite) %>% dplyr::summarise(Asite_Mean=mean((Ribo/mRNA)/TrMean),Asite_Median=median((Ribo/mRNA)/TrMean),Virus_Asite_Codon_num=n())
Occupancy=merge(Psite,Asite,by.x="Psite",by.y="Asite")
colnames(Occupancy)[1]="Codon"
Occupancy1= tidyr::gather(Occupancy, Type, virus, Psite_Mean:Asite_Median)

Psite=filter(file,grepl("ENS",Tr_ID)) %>% group_by(Psite) %>% dplyr::summarise(Psite_Mean=mean((Ribo/mRNA)/TrMean),Psite_Median=median((Ribo/mRNA)/TrMean),Host_Psite_Codon_num=n())
Asite=filter(file,grepl("ENS",Tr_ID)) %>% group_by(Asite) %>% dplyr::summarise(Asite_Mean=mean((Ribo/mRNA)/TrMean),Asite_Median=median((Ribo/mRNA)/TrMean),Host_Asite_Codon_num=n())
Occupancy=merge(Psite,Asite,by.x="Psite",by.y="Asite")
colnames(Occupancy)[1]="Codon"
Occupancy2= tidyr::gather(Occupancy, Type, host, Psite_Mean:Asite_Median)
Occupancy2$Sample=name
Occupancy2$Tr_num=length(unique(file$Tr_ID))
Host_Virus=rbind(merge(Occupancy1,Occupancy2,by=c("Codon","Type")),Host_Virus)
}
codon_table<-fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/UMAP/codon_table.txt",header=F)
colnames(codon_table)=c("Codon","AA3","AA","AAname")
Host_Virus=merge(Host_Virus,codon_table[,c(1,3)],by=c("Codon"))
head(Host_Virus)
ggplot(Host_Virus,aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  facet_grid(Type~Sample,scales = "free")+
  theme_classic()+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")

ggplot(filter(Host_Virus,grepl("zzHCMV",Sample)),aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  facet_grid(Type~Sample,scales = "free")+
  my_theme+
  labs(y="Virus P-site Occupancy",x="Host P-site Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman",size=5,col="red")

temp
i=8

Host_RNAVirus=data.frame()
for(i in c(1,6:8,10:14,17:19,21:29,45:46)){
  name=temp[i]
  name=gsub("\\.count","",name)
  file=read.table(paste0("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Ribo/Virus/Frame012/",temp[i]),sep='\t',header = T)
  file_Tr_mean=file %>% group_by(Tr_ID) %>% dplyr::summarise(TrMean=mean(Ribo/mRNA),TrMean_Ribo=mean(Ribo))
  file=merge(file,file_Tr_mean,by=c("Tr_ID"))
  Psite=filter(file,!grepl("ENS",Tr_ID)) %>% group_by(Psite) %>% dplyr::summarise(Psite_Mean=mean((Ribo)/TrMean_Ribo),Psite_Median=median((Ribo)/TrMean_Ribo),Virus_Psite_Codon_num=n())
  Asite=filter(file,!grepl("ENS",Tr_ID)) %>% group_by(Asite) %>% dplyr::summarise(Asite_Mean=mean((Ribo)/TrMean_Ribo),Asite_Median=median((Ribo)/TrMean_Ribo),Virus_Asite_Codon_num=n())
  Occupancy=merge(Psite,Asite,by.x="Psite",by.y="Asite")
  colnames(Occupancy)[1]="Codon"
  Occupancy1= tidyr::gather(Occupancy, Type, virus, Psite_Mean:Asite_Median)
  
  Psite=filter(file,grepl("ENS",Tr_ID)) %>% group_by(Psite) %>% dplyr::summarise(Psite_Mean=mean((Ribo/mRNA)/TrMean),Psite_Median=median((Ribo/mRNA)/TrMean),Host_Psite_Codon_num=n())
  Asite=filter(file,grepl("ENS",Tr_ID)) %>% group_by(Asite) %>% dplyr::summarise(Asite_Mean=mean((Ribo/mRNA)/TrMean),Asite_Median=median((Ribo/mRNA)/TrMean),Host_Asite_Codon_num=n())
  Occupancy=merge(Psite,Asite,by.x="Psite",by.y="Asite")
  colnames(Occupancy)[1]="Codon"
  Occupancy2= tidyr::gather(Occupancy, Type, host, Psite_Mean:Asite_Median)
  Occupancy2$Sample=name
  Occupancy2$Tr_num=length(unique(file$Tr_ID))
  
  Host_RNAVirus=rbind(merge(Occupancy1,Occupancy2,by=c("Codon","Type")),Host_RNAVirus)
}
Host_RNAVirus=merge(Host_RNAVirus,codon_table[,c(1,3)],by=c("Codon"))
Host_RNAVirus$Sample=gsub("PMID_34433827_","",Host_RNAVirus$Sample)
Host_RNAVirus

ggplot(Host_RNAVirus,aes(x=log10(host+0.001),y=log10(virus+0.001)))+
  geom_point()+
  facet_grid(Type~Sample,scales = "free")+
  theme_classic()+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")

write.table(rbind(Host_Virus,Host_RNAVirus),"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Virus_VS_Host_Occupancy3.txt",row.names = F,quote = F,sep = "\t")

head(Host_RNAVirus)
head(Host_Virus)

Host_Virus=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Virus_VS_Host_Occupancy3.txt")
mid=filter(Host_Virus,Tr_num>1000,virus>0) %>% group_by(Sample) %>% dplyr::summarise(count=n()) %>% filter(count==61*4)
unique(filter(Host_Virus,Tr_num>1000,virus>0)[,c(7,8)])
unique(Host_Virus$Sample)
filter(Host_Virus,grepl("VACV",Sample))
ggplot(filter(Host_Virus,Sample %in% mid$Sample),aes(x=log10(host+0.01),y=log10(virus+0.01)))+
  geom_point()+
  facet_grid(Type~Sample,scales = "free")+
  #my_theme+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "s")+
  theme_classic()



ggplot(filter(Host_Virus,grepl("HCMV\\.05hpi",Sample),Type=="Asite_Median"),aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  facet_grid(~Sample,scales = "free")+
  my_theme+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")+
  scale_x_continuous(breaks = round(seq(-0.3,0.3,0.1),2),limits = c(-0.3,0.3))+
  scale_y_continuous(breaks = round(seq(-0.3,0.3,0.1),2),limits = c(-0.3,0.3))

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/HCMV.5hpi.pdf', width= 3 , height= 3 , units='in')

head(Host_Virus)
ggplot(filter(Host_Virus,grepl("VACV.8",Sample),!grepl("zVACV.8",Sample),Type=="Asite_Median"),aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  facet_grid(~Sample,scales = "free")+
  my_theme+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")+
  scale_x_continuous(breaks = round(seq(-0.4,0.4,0.1),2),limits = c(-0.4,0.4))+
  scale_y_continuous(breaks = round(seq(-0.4,0.4,0.1),2),limits = c(-0.4,0.4))

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/VACV.8hpi.pdf', width= 5 , height= 5 , units='in')

unique(Host_Virus$Sample)

unique(mid$Sample)
ggplot(filter(Host_Virus,Sample %in% c("zzzVero.SCV2_05hpi"),Type=="Asite_Median"),aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  my_theme+
  facet_wrap(~Sample,scales = "free")+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")+
  scale_x_continuous(breaks = seq(-0.3,0.3,0.1),limits = c(-0.3,0.3))+
  scale_y_continuous(breaks = seq(-0.4,0.4,0.1),limits = c(-0.4,0.4))

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Vero.SCV2_05hpi_Psite_correlation.pdf', width= 6 , height= 6 , units='in')

ggplot(filter(Host_Virus,Sample %in% c("EToV.1hpi"),Type=="Asite_Median"),aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  my_theme+
  facet_wrap(~Sample,scales = "free")+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")+
  scale_y_continuous(breaks = seq(-0.5,0.3,0.1),limits = c(-0.5,0.3))+
  scale_x_continuous(breaks = seq(-0.3,0.2,0.1),limits = c(-0.3,0.2))
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/EToV.1hpi_Psite_correlation.pdf', width= 6 , height= 6 , units='in')


ggplot(filter(Host_Virus,Sample %in% c("MHV.5hpi"),Type=="Asite_Median"),aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  my_theme+
  facet_wrap(~Sample,scales = "free")+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")+
  scale_x_continuous(breaks = seq(-0.2,0.2,0.1),limits = c(-0.2,0.2))+
  scale_y_continuous(breaks = seq(-0.2,0.3,0.1),limits = c(-0.2,0.3))
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/MHV.5hpi_Psite_correlation.pdf', width= 6 , height= 6 , units='in')


ggplot(filter(Host_Virus,Sample %in% c("IAV.2hpi"),Type=="Asite_Median"),aes(x=log10(host),y=log10(virus)))+
  geom_point()+
  my_theme+
  facet_wrap(~Sample,scales = "free")+
  #labs(y="Occupancy")+
  geom_smooth(method = "glm",col="red",se=F)+
  #ggrepel::geom_text_repel(aes(label = paste0(AA,":",Codon)),color="blue")+
  ggpubr::stat_cor(method = "spearman")+
  scale_x_continuous(breaks = seq(-0.3,0.2,0.1),limits = c(-0.3,0.2))+
  scale_y_continuous(breaks = seq(-0.3,0.3,0.1),limits = c(-0.3,0.3))
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/IAV.2hpi_Psite_correlation.pdf', width= 7 , height= 6 , units='in')



test=filter(Host_Virus,Sample %in% mid$Sample) %>% group_by(Sample,Type) %>% dplyr::summarise(rho=cor.test(host,virus,method = "s")$estimate,
                                                                                              P=cor.test(host,virus,method = "s")$p.value)

test=filter(test,Type=="Asite_Median",!grepl("SCV2_24hpi",Sample))
test$Sample=gsub("z","0",test$Sample)
test$Sample=gsub("VACV\\.","0VACV\\.",test$Sample)
test$Sample=gsub("000Vero\\.SCV2","SCV2\\.Vero",test$Sample)
test$Sample=gsub("VACV_Early","VACV_0Early",test$Sample)

write.table(test,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Virus_VS_Host_Occupancy2.txt",row.names = F,quote = F,sep = "\t")
test=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/Virus_VS_Host_Occupancy2.txt")
ggplot(test,aes(x=Sample,y=rho)) +
  geom_bar(position="dodge", stat="identity", colour="black")+
  #facet_wrap( ~ Study,scales = "free",nrow = 4)+
  my_theme2
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/CodonUsageBias2/Ensembl/Plot/RNA_virus_Psite_correlation2.pdf', width= 15 , height= 10 , units='in')
