##data from ensembl
#cd /gpfs2/shankj/CodonUsage2/Test/ggsearch/
#wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.46.pep.all.fa.gz
#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz



##obtain the reviewed ribosomal proteins from Peixiang's list
#cd /gpfs2/shankj/CodonUsage2/Test/ggsearch/
#ls *.csv|while read id;do(sed -i 's/\r//g' $id);done
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractTargetSeq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomal_protein_Ath_cytosolic_with_canonical_transcript.csv Arabidopsis_thaliana.TAIR10.pep.all.fa  > /gpfs2/shankj/CodonUsage2/Test/ggsearch/Arabidopsis_thaliana.TAIR10.ribomalAA.fas
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractTargetSeq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomal_protein_Hs_id_with_canonical_transcript.csv  Homo_sapiens.GRCh38.pep.all.fa   > /gpfs2/shankj/CodonUsage2/Test/ggsearch/Homo_sapiens.GRCh38.ribomalAA.fas
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractTargetSeq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomal_protein_Dm_with_canonical_transcript.csv Drosophila_melanogaster.BDGP6.46.pep.all.fa  > /gpfs2/shankj/CodonUsage2/Test/ggsearch/Drosophila_melanogaster.BDGP6.46.ribomalAA.fas
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractTargetSeq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomal_protein_Mm_with_canonical_transcript.csv Mus_musculus.GRCm39.pep.all.fa  > /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mus_musculus.GRCm39.ribomalAA.fas
#cat /gpfs2/shankj/CodonUsage2/Test/ggsearch/Arabidopsis_thaliana.TAIR10.ribomalAA.fas /gpfs2/shankj/CodonUsage2/Test/ggsearch/Homo_sapiens.GRCh38.ribomalAA.fas /gpfs2/shankj/CodonUsage2/Test/ggsearch/Drosophila_melanogaster.BDGP6.46.ribomalAA.fas /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mus_musculus.GRCm39.ribomalAA.fas > /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribo.index.fas

#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomalCDS_seq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/Homo_sapiens.GRCh38.ribomalAA.fas /gpfs2/shankj/CodonUsage2/Ensembl/CDS/Homo_sapiens.GRCh38.cds.all.fa.gz > /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Homo_sapiens.GRCh38.ribomalCDS.fas
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomalCDS_seq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mus_musculus.GRCm39.ribomalAA.fas /gpfs2/shankj/CodonUsage2/Ensembl/CDS/Mus_musculus.GRCm39.cds.all.fa.gz > /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Mus_musculus.GRCm39.ribomalCDS.fas
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomalCDS_seq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/Drosophila_melanogaster.BDGP6.46.ribomalAA.fas /gpfs2/shankj/CodonUsage2/Ensembl/CDS/Drosophila_melanogaster.BDGP6.46.cds.all.fa.gz > /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Drosophila_melanogaster.BDGP6.46.ribomalCDS.fas
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomalCDS_seq.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/Arabidopsis_thaliana.TAIR10.ribomalAA.fas /gpfs2/shankj/CodonUsage2/Ensembl/Ensembl_plants/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz > /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/Arabidopsis_thaliana.TAIR10.ribomalCDS.fas

#perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/Arabidopsis_thaliana.TAIR10.ribomalCDS.fas >  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/Arabidopsis_thaliana.TAIR10.ribomalCDS.CUB.txt
#perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Homo_sapiens.GRCh38.ribomalCDS.fas >  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Homo_sapiens.GRCh38.ribomalCDS.CUB.txt
#perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Mus_musculus.GRCm39..ribomalCDS.fas >  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Mus_musculus.GRCm39.ribomalCDS.CUB.txt
#perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Drosophila_melanogaster.BDGP6.46.ribomalCDS.fas >  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/Drosophila_melanogaster.BDGP6.46.ribomalCDS.CUB.txt



cd /gpfs2/shankj/CodonUsage2/Ensembl/EnsemblProtein/animals/
ls |grep -v drosophila_melanogaster |grep -v mus_musculus |grep -v homo_sapiens  |while read id;do(
sh /gpfs2/shankj/CodonUsage2/Test/ggsearch/AnimalRiboCDS.sh $id
);done
cd  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/
ls *.ribomalCDS.fas |while read id;do(
perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl $id >  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/animals/${id%.ribomalCDS.fas}.ribomalCDS.CUB.txt
);done


cd /gpfs2/shankj/CodonUsage2/Ensembl/EnsemblProtein/plants/
ls |grep -v arabidopsis_thaliana |while read id;do(
sh /gpfs2/shankj/CodonUsage2/Test/ggsearch/PlantRiboCDS.sh $id
);done
cd  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/
ls *.ribomalCDS.fas |while read id;do(
perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl $id >  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/${id%.ribomalCDS.fas}.ribomalCDS.CUB.txt
);done


#{
##GeneTrID_seq
#gunzip *.fa.gz

##easier for checking
#ls *.fa|while read id;do(
#sed -i 's/ /__/g' $id 
#);done

#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/GeneTrID_seq.pl Mus_musculus.GRCm39.pep.all.fa > Mus_musculus.GRCm39.fas

##ggsearch36
#cd /gpfs2/shankj/CodonUsage2/Test/ggsearch/
#~/miniconda3/envs/test20240311/bin/ggsearch36  /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribo.index.fas Mus_musculus.GRCm39.fas -d 10 -E 0.01  -T 5 -p -m 8 -M 20 > ./potentialRP.ggsearch36

##evalue< 0.01 
##Identity >= 0.8 
##kept isoform's highest Identity for identical gene
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractPotentialRP.pl ./potentialRP.ggsearch36 |sort -u > potentialRP.txt
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractPotentialRP_fasta.pl Mus_musculus.GRCm39.pep.all.fa  potentialRP.txt > Mus_musculus.potentialRP.fasta
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomalCDS_seq.pl ./Mus_musculus.potentialRP.fasta /gpfs2/shankj/CodonUsage2/Ensembl/CDS/Mus_musculus.GRCm39.cds.all.fa.gz > Mus_musculus.GRCm39.ribomalCDS.fas
#perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl Mus_musculus.GRCm39.ribomalCDS.fas >  Mus_musculus.GRCm39.ribomalCDS.CUB.txt

#}


##sed -i 's/\r//g' /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mm_RP_10_wrong_hits.csv
##perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomalCDS_seq_nolike.pl /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mm_RP_10_wrong_hits.csv /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mus_musculus.potentialRP.fasta /gpfs2/shankj/CodonUsage2/Ensembl/CDS/Mus_musculus.GRCm39.cds.all.fa.gz > Mus_musculus.GRCm39.ribomalCDS.nolike.fas
##perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl Mus_musculus.GRCm39.ribomalCDS.nolike.fas > Mus_musculus.GRCm39.ribomalCDS.nolike.CUB.txt





















##ggsearch36  again
#~/miniconda3/envs/test20240311/bin/ggsearch36  Mus_musculus.potentialRP.fasta Homo_sapiens.GRCh38.pep.all.fa  -d 20  -E 0.01 -p  -T 5  -m 8 -M 20  > Mus_musculus.potentialRP.blastp.outfmt6
#~/miniconda3/envs/test20240311/bin/ggsearch36  Mus_musculus.potentialRP.fasta Drosophila_melanogaster.BDGP6.46.pep.all.fa -d 20 -p   -E 0.01 -T 5  -m 8  -M 20 >> Mus_musculus.potentialRP.blastp.outfmt6 
#~/miniconda3/envs/test20240311/bin/ggsearch36  Mus_musculus.potentialRP.fasta Arabidopsis_thaliana.TAIR10.pep.all.fa -d 20 -p   -T 5  -m 8  -E 0.01 -M 20 >> Mus_musculus.potentialRP.blastp.outfmt6 


##Extract Ribosomal Protein fasta
#ls *.fa|grep -v Mus_musculus|while read id;do(grep '>' $id);done |grep -v -i 'ribosomal__protein'|grep -i 'transcript_biotype:protein_coding' > unribosomal.protein.txt
#perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractRP.pl Mus_musculus.GRCm39.pep.all.fa unribosomal.protein.txt  Mus_musculus.potentialRP.blastp.outfmt6  > Mus_musculus.RP.final2.fasta

##check
#awk -F ',' '{print $2}' /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mm_RP_missed_hit.csv |sed 's/\"//g'|sed 's/\r//g' |while read id;do(echo $id; grep $id /gpfs2/shankj/CodonUsage2/Test/ggsearch/potentialRP.outfmt6|awk '{print $3,$4,$8,$9}');done
#awk -F ',' '{print $2}' /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mm_RP_wrong_hit.csv |sed 's/\"//g'|sed 's/\r//g' |while read id;do(echo $id; grep $id /gpfs2/shankj/CodonUsage2/Test/ggsearch/Mus_musculus.GRCm39.pep.all.fa);done



#~/miniconda3/envs/test20240311/bin/ggsearch36  Mus_musculus.potentialRP.fasta ribo.index.fas -d 20  -E 0.01 -p  -T 5  -m 8   > test.blastp.outfmt6
