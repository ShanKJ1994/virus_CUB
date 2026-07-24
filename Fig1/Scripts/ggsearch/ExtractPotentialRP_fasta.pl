#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
##############################################################################

#ENSMUST00000124102.8 cds chromosome:GRCm39:16:43929986:43959381:-1 gene:ENSMUSG00000052459.14 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Atp6v1a description:ATPase, H+ transporting, lysosomal V1 subunit A [Source:MGI Symbol;Acc:MGI:1201780]

my $file=shift;#protein.fasta
$/="\n>";
my %seq;
my %cds;
open(IN,$file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
	$line=~s/\>//;
	my @seq=split /\n/,$line;
	my $seq;
	for (my $i=1;$i<=$#seq ;$i++) {
		$seq.=$seq[$i];
	}
	my @gene=split /__/,$seq[0];
	my $ID;
	foreach my $name (@gene) {
		if ($name=~/transcript:/) {
			$ID=$name;
			$ID=~s/transcript://g;
		}
	}
	$seq{$ID}=$seq;
	$cds{$ID}=$seq[0];
}
close IN;


##############################################################################

#read in  blast ID file
$/="\n";
my $file2=shift;
open(IN2,$file2) or die ("can not open $file2\n");
while (my $line=<IN2>) {
	chomp $line;
	my @gene=split /:/,$line;
	if (exists $seq{$gene[1]} &&
		$cds{$gene[1]}!~/pseudo/i && 
		$cds{$gene[1]}!~/chloroplas/i && $cds{$gene[1]}!~/plast/i && #$cds{$gene[0]}=~/ribosomal protein/i &&
		$cds{$gene[1]}!~/itochondrial/i && $cds{$gene[1]}!~/chondrial/i && #$cds{$gene[0]}!~/ubiquitin/i && 
		$cds{$gene[1]}!~/retro/i &&
		$cds{$gene[1]}=~/transcript_biotype:protein_coding/i) {
		print ">$cds{$gene[1]}\n";
		print "$seq{$gene[1]}\n";
	}
	#my @gene=split /__/,$line;
	#foreach my $name(@gene) {
		#if ($name=~/gene:/) {
			#$name=~s/gene://;
			#if (exists $seq{$name}) {
				#print ">$line";
				#print "$seq{$name}\n";
			#}
		#}
	#}
}
close IN2;