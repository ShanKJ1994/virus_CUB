#data from ensembl
cd /gpfs2/shankj/CodonUsage2/Test/
wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.46.pep.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz

gunzip *.fa.gz

#easier for checking
ls *.fa|while read id;do(
sed -i 's/ /__/g' $id &
);done

#blast index
cd /gpfs2/shankj/CodonUsage2/Test/
makeblastdb -in Mus_musculus.GRCm39.pep.all.fa -dbtype prot -out Mus_musculus_db &
makeblastdb -in Homo_sapiens.GRCh38.pep.all.fa -dbtype prot -out Homo_sapiens_db &
makeblastdb -in Drosophila_melanogaster.BDGP6.46.pep.all.fa -dbtype prot -out Drosophila_melanogaster_db &
makeblastdb -in Arabidopsis_thaliana.TAIR10.pep.all.fa -dbtype prot -out Arabidopsis_thaliana_db &

#obtain the reviewed ribosomal proteins from Peixiang's list
ls *.csv|while read id;do(sed -i 's/\r//g' $id);done
perl /gpfs2/shankj/CodonUsage2/Test/ExtractTargetSeq.pl /gpfs2/shankj/CodonUsage2/Test/ribosomal_protein_Ath_cytosolic_with_canonical_transcript.csv Arabidopsis_thaliana.TAIR10.pep.all.fa  > ribo.index.fas
perl /gpfs2/shankj/CodonUsage2/Test/ExtractTargetSeq.pl /gpfs2/shankj/CodonUsage2/Test/ribosomal_protein_Hs_id_with_canonical_transcript.csv  Homo_sapiens.GRCh38.pep.all.fa   >> ribo.index.fas
perl /gpfs2/shankj/CodonUsage2/Test/ExtractTargetSeq.pl /gpfs2/shankj/CodonUsage2/Test/ribosomal_protein_Dm_with_canonical_transcript.csv Drosophila_melanogaster.BDGP6.46.pep.all.fa  >> ribo.index.fas

#longest peptide
perl /gpfs2/shankj/CodonUsage2/Test/longeastPep.pl Mus_musculus.GRCm39.pep.all.fa > Mus_musculus.GRCm39.longest.fasta

cd /gpfs2/shankj/CodonUsage2/Test/
#blast
#blastp -query ribo.index.fas \
   #-db Mus_musculus_db   \
   #-outfmt 6 -evalue 1e-5 -num_threads 3 > potentialRP.outfmt6 
  #bioawk -c fastx '{print $name"\t"length($seq)}' ${id}ta > ./${id}.length

#ggsearch36
cd /gpfs2/shankj/CodonUsage2/Test/
~/miniconda3/envs/test20240311/bin/ggsearch36  ribo.index.fas Mus_musculus.GRCm39.pep.all.fa  -d 10 -E 0.01  -T 5 -p -m 8  > /gpfs2/shankj/CodonUsage2/Test/potentialRP.ggsearch36


#check
awk -F ',' '{print $2}' /gpfs2/shankj/CodonUsage2/Test/Mm_RP_missed_hit.csv |sed 's/\"//g'|sed 's/\r//g' |while read id;do(echo $id; grep $id /gpfs2/shankj/CodonUsage2/Test/potentialRP.outfmt6|awk '{print $3,$4,$8,$9}');done
awk -F ',' '{print $2}' /gpfs2/shankj/CodonUsage2/Test/Mm_RP_wrong_hit.csv |sed 's/\"//g'|sed 's/\r//g' |while read id;do(echo $id; grep $id /gpfs2/shankj/CodonUsage2/Test/Mus_musculus.GRCm39.pep.all.fa);done



#evalue< 1e-5 Identity >= 0.8 
#$mid[1]!~/pseudo/i && $mid[10]<0.00001 && $mid[2]/$mid[3]>=0.8 && ($mid[9]-$mid[8]+1)/$mid[3]>=0.8
perl /gpfs2/shankj/CodonUsage2/Test/ExtractPotentialRP.pl /gpfs2/shankj/CodonUsage2/Test/potentialRP.ggsearch36 |sort -u > potentialRP.txt
perl /gpfs2/shankj/CodonUsage2/Test/ExtractPotentialRP_fasta.pl Mus_musculus.GRCm39.pep.all.fa  potentialRP.txt > Mus_musculus.potentialRP.fasta


#blast again
#blastp -query Mus_musculus.potentialRP.fasta \
   #-db Homo_sapiens_db  \
   #-outfmt 6 -evalue 1e-5 -num_threads 3 > Mus_musculus.potentialRP.blastp.outfmt6 

#blastp -query Mus_musculus.potentialRP.fasta \
   #-db Drosophila_melanogaster_db   \
   #-outfmt 6 -evalue 1e-5 -num_threads 3 >> Mus_musculus.potentialRP.blastp.outfmt6 

#blastp -query Mus_musculus.potentialRP.fasta \
   #-db Arabidopsis_thaliana_db   \
   #-outfmt 6 -evalue 1e-5 -num_threads 3 >> Mus_musculus.potentialRP.blastp.outfmt6 

#ggsearch36  again
~/miniconda3/envs/test20240311/bin/ggsearch36  Mus_musculus.potentialRP.fasta Homo_sapiens.GRCh38.pep.all.fa -d 1 -E 0.01 -p  -T 5  -m 8   > Mus_musculus.potentialRP.blastp.outfmt6
~/miniconda3/envs/test20240311/bin/ggsearch36  Mus_musculus.potentialRP.fasta Drosophila_melanogaster.BDGP6.46.pep.all.fa -d 1 -p -E 0.01  -T 5  -m 8   >> Mus_musculus.potentialRP.blastp.outfmt6 
~/miniconda3/envs/test20240311/bin/ggsearch36  Mus_musculus.potentialRP.fasta Arabidopsis_thaliana.TAIR10.pep.all.fa -d 1 -p -E 0.01  -T 5  -m 8   >> Mus_musculus.potentialRP.blastp.outfmt6 

#Extract Ribosomal Protein fasta
cat *.csv |grep -v merge |grep -v missed |grep -v wrong > merge.csv
perl /gpfs2/shankj/CodonUsage2/Test/ExtractRP.pl Mus_musculus.GRCm39.pep.all.fa /gpfs2/shankj/CodonUsage2/Test/ribo.index.fas  Mus_musculus.potentialRP.blastp.outfmt6  > Mus_musculus.RP.final2.fasta
awk -F ',' '{print $2}' /gpfs2/shankj/CodonUsage2/Test/Mm_RP_missed_hit.csv |sed 's/\"//g'|sed 's/\r//g' |while read id;do(echo $id; grep $id Mus_musculus.RP.final2.fasta);done
