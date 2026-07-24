#! /bin/bash

#align to genome

output_dir=/lustre/user/lulab/gaopx/projects/virus_host_translation/RiboSeq/VACV_HeLaS3
virus_genome_dir=/lustre/user/lulab/gaopx/projects/virus_host_translation/RiboSeq/VACV_HeLaS3/virus_genome/index

for i in /gpfs2/shankj/riboSeq/VACV_PMID28003488/*/*mRNA*_clean_tRNA.fq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id="${sample_id%_clean_tRNA.fq.gz}-host"
	out_prefix="${output_dir}/alignment_host/${sample_id}_"
	STAR --runThreadN 32 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
	 --genomeDir /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/index --readFilesCommand zcat \
	 --readFilesIn $i --outFileNamePrefix $out_prefix \
	 --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

for i in $(ls /gpfs2/shankj/riboSeq/VACV_PMID28003488/*/*_clean_tRNA.fq.gz | grep -v mRNA)
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id="${sample_id%_clean_tRNA.fq.gz}-host"
  out_prefix="${output_dir}/alignment_host/${sample_id}_"
  STAR --runThreadN 32 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
   --genomeDir /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/index \
   --readFilesCommand zcat \
   --alignEndsType EndToEnd \
   --readFilesIn $i --outFileNamePrefix $out_prefix \
   --outSAMtype BAM SortedByCoordinate --outSAMattributes All \
   --quantMode TranscriptomeSAM
done

STAR --genomeLoad Remove --genomeDir /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/index \
  --outFileNamePrefix remove_genome_output
rm remove_genome_output*

for i in ${output_dir}/alignment_host/*_Aligned.sortedByCoord.out.bam
do
  samtools index $i
done


# align to virus genome
for i in /gpfs2/shankj/riboSeq/VACV_PMID28003488/*/*mRNA*_clean_tRNA.fq.gz
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id="${sample_id%_clean_tRNA.fq.gz}-virus"
	out_prefix="${output_dir}/alignment_virus/${sample_id}_"
	STAR --runThreadN 32 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
	 --genomeDir ${virus_genome_dir} --readFilesCommand zcat \
	 --readFilesIn $i --outFileNamePrefix $out_prefix \
	 --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done

for i in $(ls /gpfs2/shankj/riboSeq/VACV_PMID28003488/*/*_clean_tRNA.fq.gz | grep -v mRNA)
do
  [[ -e "$i" ]] || break
  sample_id=${i##*/}
  sample_id="${sample_id%_clean_tRNA.fq.gz}-virus"
  out_prefix="${output_dir}/alignment_virus/${sample_id}_"
  STAR --runThreadN 32 --genomeLoad LoadAndKeep --limitBAMsortRAM 64000000000 \
   --genomeDir ${virus_genome_dir} \
   --readFilesCommand zcat \
   --alignEndsType EndToEnd \
   --readFilesIn $i --outFileNamePrefix $out_prefix \
   --outSAMtype BAM SortedByCoordinate --outSAMattributes All \
   --quantMode TranscriptomeSAM
done

STAR --genomeLoad Remove --genomeDir ${virus_genome_dir} \
  --outFileNamePrefix remove_genome_output
rm remove_genome_output*

for i in ${output_dir}/alignment_virus/*_Aligned.sortedByCoord.out.bam
do
  samtools index $i
done

mkdir /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/pic
~/miniconda3/envs/Rstudio/bin/Rscript /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/get_psite_offset.R

ln -s /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_host/uninfect_mRNA-host_Aligned.sortedByCoord.out.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/merge/uninfect_mRNA_Aligned.sortedByCoord.out.bam
ln -s /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_host/uninfect_Ribo-host_Aligned.sortedByCoord.out.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/merge/uninfect_Ribo_Aligned.sortedByCoord.out.bam

mkdir /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/merge
cd /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/merge
samtools merge -o mRNA_VACV_Early.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_host/VACV_Early_mRNA-host_Aligned.sortedByCoord.out.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_virus/VACV_Early_mRNA-virus_Aligned.sortedByCoord.out.bam 
samtools merge -o mRNA_VACV_8hpi.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_host/VACV_mRNA-host_Aligned.sortedByCoord.out.bam  /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_virus/VACV_mRNA-virus_Aligned.sortedByCoord.out.bam 
samtools merge -o Ribo_VACV_Early.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_host/VACV_Early_Ribo-host_Aligned.sortedByCoord.out.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_virus/VACV_Early_Ribo-virus_Aligned.sortedByCoord.out.bam 
samtools merge -o Ribo_VACV_8hpi.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_host/VACV_Ribo-host_Aligned.sortedByCoord.out.bam /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/alignment_virus/VACV_Ribo-virus_Aligned.sortedByCoord.out.bam 

for i in /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/merge/*.bam
do
  samtools index $i
done




python /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/merged_bam_disambiguation_script.py /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID28003488/VACV.fasta  /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/merge    /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation         &


cat /lustre/user/lulab/gaopx/projects/virus_host_translation/RiboSeq/VACV_HeLaS3/virus_genome/VACV.gtf /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/gtf/Homo_sapiens.GRCh38.104.gtf  > /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/VACV.Human.gtf


nohup featureCounts -T 16 -s 1 -O -M -t CDS -g gene_id --fraction --fracOverlap 0.5 \
  -a /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/VACV.Human.gtf \
  -o featureCounts_VACV_HeLa229_SciAdv_S1.txt \
  /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation/*.bam \
  > featureCounts_VACV_HeLa229_SciAdv_S1.log 2>&1 &

nohup featureCounts -T 16 -s 2 -O -M -t CDS -g gene_id --fraction --fracOverlap 0.5 \
  -a /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/VACV.Human.gtf \
  -o featureCounts_VACV_HeLa229_SciAdv_S2.txt \
  /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation/*.bam \
  > featureCounts_VACV_HeLa229_SciAdv_S2.log 2>&1 &


cd /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/featureCount/
ls *.csv|while read id;do(
perl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/featureCount/TrID_GeneExpression.pl longest_cds.csv $id > ${id%.csv}_TrID.csv
);done

for i in /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation/*.bam
do
  samtools index $i
done

rm -r /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/merge/

#mRNA  
mkdir /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/Coverage
cd /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/
cd /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation/
ls *.disambiguated.bam|grep mRNA |while read id;do(
python /lustre/user/lulab/gaopx/tools/bamCoverage_modified/bamCoverage_fractionByNH/bamCoverage_fractionByNH.py --maxMultipleAlignSite 10 -b $id -o ../Coverage/${id%.disambiguated.bam}_fw.bw --filterRNAstrand forward -bs 1 -p 1  & 
python /lustre/user/lulab/gaopx/tools/bamCoverage_modified/bamCoverage_fractionByNH/bamCoverage_fractionByNH.py --maxMultipleAlignSite 10 -b $id -o ../Coverage/${id%.disambiguated.bam}_rc.bw --filterRNAstrand reverse -bs 1 -p 1  &
);done

cd /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/Coverage
python /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/UniqueMapped/bam_coverage_NH10.py  -u 32 -l 26 -v /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/uninfect_Ribo.txt  /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation/uninfect_Ribo.disambiguated.bam uninfect_Ribo &
python /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/UniqueMapped/bam_coverage_NH10.py  -u 32 -l 26 -v /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/VACV_Early_Ribo.txt  /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation/Ribo_VACV_Early.disambiguated.bam VACV_Early_Ribo &
python /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/UniqueMapped/bam_coverage_NH10.py  -u 32 -l 26 -v /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/VACV_Ribo.txt  /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/VACV_PMID37285440_SciAdv/disambiguation/Ribo_VACV_8hpi.disambiguated.bam VACV_8hpi_Ribo &

ls *_fw.bw|grep -v mRNA|while read id;do(
   python /gpfs2/zhangh/backup/mysrc/bigwig_coverage_per_nt.py $id ${id%fw.bw}rc.bw /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/VACV.Human.bed > ${id%fw.bw}transcript.count &
);done

ls *_fw.bw|grep mRNA|while read id;do(
   python /gpfs2/zhangh/backup/mysrc/bigwig_coverage_per_nt.py ${id%fw.bw}rc.bw $id /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/VACV.Human.bed > ${id%fw.bw}transcript.count &
);done


#/gpfs2/shankj/software/ORF_RATER/gtfToGenePred /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/gtf/Homo_sapiens.GRCh38.104.gtf  stdout | /gpfs2/shankj/software/ORF_RATER/genePredToBed stdin /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/Human.bed
gffread -x /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/Human.CDS.fasta  -g  /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa /lustre/user/lulab/gaopx/genome/GRCh38/ensemble/gtf/Homo_sapiens.GRCh38.104.gtf 
perl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/Final.transcript.pos.Psite.Asite.pl ./annot_longest_cds.csv /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/Human.CDS.fasta  > /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/Human.CDS.codon.pos.txt

#/gpfs2/shankj/software/ORF_RATER/gtfToGenePred /lustre/user/lulab/gaopx/projects/virus_host_translation/RiboSeq/VACV_HeLaS3/virus_genome/VACV.gtf stdout | /gpfs2/shankj/software/ORF_RATER/genePredToBed stdin /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV.bed
gffread -x /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV.CDS.fasta  -g  /lustre/user/lulab/gaopx/projects/virus_host_translation/RiboSeq/VACV_HeLaS3/virus_genome/VACV.fasta /lustre/user/lulab/gaopx/projects/virus_host_translation/RiboSeq/VACV_HeLaS3/virus_genome/VACV.gtf
perl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/Final.transcript.pos.Psite.Asite.pl /gpfs2/shankj/riboSeq/VACV_PMID28003488/VACV_Human/VACV_Human_gene_table.csv /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV.CDS.fasta    > /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV.CDS.codon.pos.txt
grep -v ENS  /gpfs2/shankj/riboSeq/VACV_PMID28003488/VACV_Human/VACV_Human_gene_table.csv  > /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV_gene_table.csv

#For gmRC raw data
perl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/CDS_read_per_codon_Virus.pl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/Human.CDS.codon.pos.txt,/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV.CDS.codon.pos.txt uninfect_mRNA_transcript.count uninfect_Ribo_transcript.count &
perl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/CDS_read_per_codon_Virus.pl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/Human.CDS.codon.pos.txt,/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV.CDS.codon.pos.txt mRNA_VACV_Early_transcript.count VACV_Early_Ribo_transcript.count &
perl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus/CDS_read_per_codon_Virus.pl /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/Human.CDS.codon.pos.txt,/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/Virus_Annotation_Expression/TrPosCodon/VACV.CDS.codon.pos.txt mRNA_VACV_8hpi_transcript.count VACV_8hpi_Ribo_transcript.count &
