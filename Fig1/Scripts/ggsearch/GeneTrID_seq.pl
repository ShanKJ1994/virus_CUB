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
	my @gene=split /__/,$seq[0];
	my $gene;
	my $tr_id;
	foreach my $name(@gene) {
		if ($name=~/gene:/) {
			$name=~s/gene://;
			$gene=$name;
		}
	}
	foreach my $name(@gene) {
		if ($name=~/transcript:/) {
			$name=~s/transcript://;
			$tr_id=$name;
		}
	}
	my $seq;
	for (my $i=1;$i<=$#seq ;$i++) {
		$seq.=$seq[$i];
	}
	print ">$gene:$tr_id\n$seq\n";
}
close IN;
