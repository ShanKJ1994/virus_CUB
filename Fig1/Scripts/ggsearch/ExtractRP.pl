#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
##############################################################################

#ENSMUST00000124102.8 cds chromosome:GRCm39:16:43929986:43959381:-1 gene:ENSMUSG00000052459.14 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Atp6v1a description:ATPase, H+ transporting, lysosomal V1 subunit A [Source:MGI Symbol;Acc:MGI:1201780]

my $file=shift;#protein.fasta
$/="\n>";
my %cds;
open(IN,$file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
	$line=~s/\>//;
	my @seq=split /\n/,$line;
	$cds{$seq[0]}=$seq[1];
}
close IN;

##############################################################################

$/="\n";
#read in ID file
my $file2=shift;
my %ID;
open(IN,$file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	if ($line=~/\>/) {
		$line=~s/\>//;
		chomp $line;
		my @mid=split /__/,$line;
		$ID{$mid[0]}=1;
	}
}
close IN;

##############################################################################

#read in  blast ID file
$/="\n";
my $file3=shift;
my %RP;
my %no;
#0          1   2        3    4         5         6   7       8   9     10     11
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#read in blast file
open(IN,$file3) or die ("can not open $file3\n");
while (my $line=<IN>) {
	my @mid=split /\s+/,$line;
	if ($mid[0]!~/pseudo/i && $mid[10]<0.01 && $mid[2]>=70  &&
		$mid[0]!~/chloroplas/i && $mid[0]!~/plast/i && #$mid[0]=~/ribosomal protein/i &&
		$mid[0]!~/itochondrial/i && $mid[0]!~/chondrial/i && #$mid[0]!~/ubiquitin/i && 
		$mid[0]!~/retro/i &&
		$mid[0]=~/transcript_biotype:protein_coding/i) {
		my @gene=split /__/,$mid[0];
		my @gene1=split /__/,$mid[1];
		
		if (exists $ID{$gene1[0]}) {
			$no{$mid[0]}=1;
			print $line;
		}
		$RP{$mid[0]}=1;
	}
}
close IN;
#foreach my $ID(keys %RP) {
	#unless (exists $no{$ID}) {
		#print ">$ID\n$cds{$ID}\n";
	#}
#}
