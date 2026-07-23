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
	chomp  $line;
	my @mid=split /,/,$line;
	$mid[-1]=~s/\"//g;
	$ID{$mid[-1]}=1;
	#print $mid[-1]."\n";
}
close IN;

##############################################################################

#read in  fas file
$/="\n>";
my $file=shift;
open(IN2,$file) or die ("can not open $file\n");
while (my $line=<IN2>) {
	chomp $line;
	$line=~s/\>//;
	my @mid=split /\n/,$line;
	my @name=split /__|:/,$mid[0];
	my $count=0;
	foreach my $gene(@name) {
		#print "$gene\n";
		if (exists $ID{$gene}) {
			$count=1000;
			#print "$gene\n";
			last;
		}
	}
	if ($count == 1000) {
		print ">$mid[0]\n";
		my $seq2;
		for (my $j=1;$j<=$#mid ;$j++) {
			$seq2.=$mid[$j];
		}
		print "$seq2\n";
	}
}
close IN2;





