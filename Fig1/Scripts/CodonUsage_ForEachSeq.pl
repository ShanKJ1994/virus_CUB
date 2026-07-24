#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

##############################################################################

#read in lognest CDS fasta file

$/="\n>";
my $file=shift;
my %CDS;
open(IN2,$file) or die ("can not open $file\n");
while (my $line=<IN2>) {
	chomp $line;
	$line=~s/\>//;
	my @mid=split /\n/,$line;
	my $seq;
	for (my $i=1;$i<=$#mid;$i++) {
		$seq.=$mid[$i];
	}
	if ( length($seq) % 3 ==0 ) {
		my @seq=split //,uc($seq);
		for (my $i=0;$i<=$#seq ;$i+=3) {
			if ($seq[$i]=~/[ATCG]/i && $seq[$i+1]=~/[ATCG]/i && $seq[$i+2]=~/[ATCG]/i) {
			$CDS{"$seq[$i]$seq[$i+1]$seq[$i+2]"}{$mid[0]}+=1;
			}
		}
	}
}
close IN2;

foreach my $codon(sort keys %CDS) {
	foreach  my $seq (sort keys %{$CDS{$codon}}) {
		print "$codon\t$CDS{$codon}{$seq}\t$seq\n";
	}
	
}

