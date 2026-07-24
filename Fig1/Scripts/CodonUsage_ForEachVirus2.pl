#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

##############################################################################

#read in viral CDS fas file
#0                 1          2           3      4         5        6           7                      
#Organism Name| Species| Molecule type| Host| Assembly| BioSample| BioProject| SRA Accession| GenBank Title
$/="\n>";
my $file=shift;
my %CDS;
my %test;
open(IN2,$file) or die ("can not open $file\n");
while (my $line=<IN2>) {
	$line=~s/\>//;
	my @mid=split /\n/,$line;
	my @virus=split /\|/,$mid[0];
	my $seq2;
	for (my $j=1;$j<=$#mid ;$j++) {
		$seq2.=$mid[$j];
	}                
	#unless ($line=~/coronavirus/i && $virus[8]=~/1a /i) {
	#if ($line=~/GCF_000866645/) {
		if (length($seq2) % 3 ==0 && $mid[1]=~/^ATG/ && $mid[0] !~/phage/i) {
			my @seq=split //,$seq2;
			for (my $i=0;$i<=$#seq ;$i+=3) {
				if ($seq[$i]=~/[ATCG]/i && $seq[$i+1]=~/[ATCG]/i && $seq[$i+2]=~/[ATCG]/i) {
					my $virus;
					my $name1="$virus[0]:$virus[4]";
					my $name2="$virus[0]:$virus[1]:$virus[2]:$virus[3]:$virus[4]";
					if (exists $test{$name1}) {
						$virus=$test{$name1};
					}else{
						$test{$name1}=$name2;
						$virus=$name2;
					}
					$CDS{"$seq[$i]$seq[$i+1]$seq[$i+2]"}{$virus}+=1;
				}
			}
		}
	#}
}
close IN2;

foreach my $codon(sort keys %CDS) {
	foreach  my $seq (sort keys %{$CDS{$codon}}) {
		print "$codon\t$CDS{$codon}{$seq}\t$seq\n";
	}
	
}



