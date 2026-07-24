#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
##############################################################################

#read in ID file
my $file2=shift;
my %ID;
open(IN,$file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	if ($line=~/\>/) {
		$line=~s/\>//;
		chomp $line;
		my @mid=split /__/,$line;
		foreach my $tr_id(@mid) {
			if ($tr_id=~/transcript:/) {
				$tr_id=~s/transcript://;
				$ID{$tr_id}=1;
			}
		}
		
	}
}
close IN;

##############################################################################

#ENSMUST00000124102.8 cds chromosome:GRCm39:16:43929986:43959381:-1 gene:ENSMUSG00000052459.14 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Atp6v1a description:ATPase, H+ transporting, lysosomal V1 subunit A [Source:MGI Symbol;Acc:MGI:1201780]

my $file=shift;#/gpfs2/shankj/CodonUsage2/Ensembl/CDS/Mus_musculus.GRCm39.cds.all.fa.gz
$/="\n>";
open(IN,"gunzip -c $file |" ) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
	$line=~s/\>//;
	my @seq=split /\n/,$line;
	my @gene=split /\s+/,$seq[0];
	my $seq;
	for (my $i=1;$i<=$#seq ;$i++) {
		$seq.=$seq[$i];
	}
	if (exists $ID{$gene[0]}) {
		print ">$seq[0]\n$seq\n";
		
	}
}
close IN;