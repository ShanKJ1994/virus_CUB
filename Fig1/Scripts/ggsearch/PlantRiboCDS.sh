cd /gpfs2/shankj/CodonUsage2/Ensembl/EnsemblProtein/plants/$1/pep/
#GeneTrID_seq
gunzip *.fa.gz

#easier for checking
ls *.pep.all.fa|while read id;do(
sed -i 's/ /__/g' $id 

perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/GeneTrID_seq.pl $id > ${id%.pep.all.fa}.fas

#ggsearch36
~/miniconda3/envs/test20240311/bin/ggsearch36  /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribo.index.fas ./${id%.pep.all.fa}.fas   -d 10 -E 0.01 -T 5 -p -m 8 -M 20 > ./potentialRP.ggsearch36
#evalue< 0.01 
#Identity >= 0.8 
#kept isoform's highest Identity for identical gene
perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractPotentialRP.pl ./potentialRP.ggsearch36 |sort -u > potentialRP.txt
perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ExtractPotentialRP_fasta.pl $id  potentialRP.txt > ${id%.pep.all.fa}.potentialRP.fasta
perl /gpfs2/shankj/CodonUsage2/Test/ggsearch/ribosomalCDS_seq.pl ${id%.pep.all.fa}.potentialRP.fasta /gpfs2/shankj/CodonUsage2/Ensembl/Ensembl_plants/$1/cds/${id%.pep.all.fa}.cds.all.fa.gz > /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/${id%.pep.all.fa}.ribomalCDS.fas
#perl /gpfs2/shankj/CodonUsage2/bin/CodonUsage.pl /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/${id%.pep.all.fa}.ribomalCDS.fas >  /gpfs2/shankj/CodonUsage2/Ensembl/ribosomalProtein/plants/${id%.pep.all.fa}.ribomalCDS.CUB.txt

);done

