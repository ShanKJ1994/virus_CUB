#####################################################################

#Host CUB from ggsearch

#####################################################################

cd /gpfs2/shankj/CodonUsage2/Ensembl/CDS/
ls |while read id;do(
perl /gpfs2/shankj/CodonUsage2/bin/LongestCDS.pl $id > /gpfs2/shankj/CodonUsage2/Ensembl/LongestCDS/Animal/${id%.cds.all.fa.gz}.longest.txt
);done


cd /gpfs2/shankj/CodonUsage2/Ensembl/Ensembl_plants/
ls /gpfs2/shankj/CodonUsage2/Ensembl/Ensembl_plants/*/*/*.gz |while read id;do(
perl /gpfs2/shankj/CodonUsage2/bin/LongestCDS.pl $id > /gpfs2/shankj/CodonUsage2/Ensembl/LongestCDS/Plant/`basename ${id%.cds.all.fa.gz}.longest.txt`
);done


sh /gpfs2/shankj/CodonUsage2/Test/ggsearch/ggsearch_test.sh
sh /gpfs2/shankj/CodonUsage2/Ensembl/NCBI/NCBI_ribosomalprotein.sh

cd /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein
cat ./animals/*.txt |sed 's/.ribomalCDS.fas//g'  > /gpfs2/shankj/CodonUsage2/Ensembl/Animal.ribosomalprotein.CUB.txt &
cat ./plants/*.txt |sed 's/.ribomalCDS.fas//g'  > /gpfs2/shankj/CodonUsage2/Ensembl/Plant.ribosomalprotein.CUB.txt &

#####################################################################

#Host longest CDS CUB

#####################################################################

ls /gpfs2/shankj/CodonUsage2/Ensembl/LongestCDS/*/* |while read id;do(
perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage_ForEachSeq.pl $id > /gpfs2/shankj/CodonUsage2/Ensembl/LongestCDS_CUB/`basename ${id%.longest.txt}.CUB.txt` 
);done


#####################################################################

#virus CUB

#####################################################################

#Accession,Organism_Name,Submitters,Organization,Org_location,Release_Date,Isolate,Species,Molecule_type,Length,Nuc_Completeness,Geo_Location,USA,Host,Isolation_Source,Collection_Date
~/miniconda3/envs/monocle_env/bin/Rcript /gpfs2/shankj/CodonUsage2/Ensembl/Virus.Host.R

cd /gpfs2/shankj/CodonUsage2/Ensembl/
#cat /gpfs2/shankj/CodonUsage2/Ensembl/VirusCDS/sequences.fasta >> /gpfs2/shankj/CodonUsage2/Ensembl/VirusCDS/TotalRefseq.fasta
perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage_ForEachVirus2.pl /gpfs2/shankj/CodonUsage2/Ensembl/VirusCDS/TotalRefseq.fasta > /gpfs2/shankj/CodonUsage2/Ensembl/Virus.CDS.CUB.txt
perl /gpfs2/shankj/CodonUsage2/bin/VirusHostPre.pl  /gpfs2/shankj/CodonUsage2/Ensembl/Virus.CDS.CUB.txt > /gpfs2/shankj/CodonUsage2/Ensembl/VirusCDS/AllVirus_CUB.Final.txt

#only kept viruses whose host exists in Ensembl and NCBI, although similar results
perl /gpfs2/shankj/CodonUsage2/bin/VirusHost.pl > /gpfs2/shankj/CodonUsage2/Ensembl/EnsemblHost.virus.CUB.txt
perl /gpfs2/shankj/CodonUsage2/bin/VirusHost_matchedSpecies.pl  |awk -F '\t' '{print $3"\t"$1}'|sort -u > /gpfs2/shankj/CodonUsage2/Ensembl/EnsemblHost.virus.CUB.species.txt

#all host not only our hosts download from Ensemblc
perl /gpfs2/shankj/CodonUsage2/Ensembl/VirusHost_allhost.pl|awk -F '\t' '{print $1"\t"$2"\t"$3}' |sort -u > /gpfs2/shankj/CodonUsage2/Ensembl/VirusHost_allhost.txt
perl /gpfs2/shankj/CodonUsage2/Ensembl/VirusHost_matchedSpecies_allhost.pl|awk -F '\t' '{print $1"\t"$2"\t"$3}' |sort -u > /gpfs2/shankj/CodonUsage2/Ensembl/VirusHost_allhost_matchedSpecies.txt

