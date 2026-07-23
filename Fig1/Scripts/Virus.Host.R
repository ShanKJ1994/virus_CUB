library(data.table)
# 读取CSV文件
data <- fread("/gpfs2/shankj/CodonUsage2/Ensembl/VirusCDS/sequences.csv",header=T)
data2 <- fread("/gpfs2/shankj/CodonUsage2/Ensembl/VirusCDS/SCV2.csv",header=T)

# 提取Organism_Name, Species, Host, Molecule_type四列数据
selected_data <- subset(data, select = c("Organism_Name", "Species", "Host", "Molecule_type"))
selected_data2 <- subset(data2, select = c("Organism_Name", "Species", "Host", "Molecule_type"))

# 去除重复行
unique_data <- unique(rbind(selected_data,selected_data2))
write.table(unique_data,"/gpfs2/shankj/CodonUsage2/Ensembl/Virus.Host.txt",row.names = F,quote = F,sep = "\t")
