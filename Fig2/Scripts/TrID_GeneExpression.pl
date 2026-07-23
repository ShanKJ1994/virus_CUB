#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

#transcript_id,gene_id,cds_length,transcript_length
#ENST00000373020,ENSG00000000003,735,3768
#ENST00000373031,ENSG00000000005,951,1205
my $file=shift;#Tr_Gene
my %Tr_Gene;
open(IN1, $file) or die ("can not open $file\n");
while (my $line=<IN1>) {
	$line=~s/transcript://g;
	$line=~s/gene://g;
	my @gene=split /,/,$line;
	$Tr_Gene{$gene[1]}=$gene[0];
}
close IN1;



my $Express=shift;#Tr_Gene
open(IN1, $Express) or die ("can not open $Express\n");
while (my $line=<IN1>) {
	if ($line=~/Length/) {
		print "TrID,$line";
	}else{
		if ($line=~/^ENS/) {
			my @gene=split /,/,$line;
			print "$Tr_Gene{$gene[0]},$line";
		}else{
			my @gene=split /,/,$line;
			print "$gene[0],$line";
		}
	}
}
close IN1;
