#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
##############################################################################
#                         0                 1                       2       3        4
#AAA     19      Abaca bunchy top virus :Abaca bunchy top virus:ssDNA(+):Musa sp.:GCF_000872625.1::PRJNA485481:
#AAA     2197    Abalone herpesvirus Victoria/AUS/2009 :Aurivirus haliotidmalaco1:dsDNA::GCF_000900375.1::PRJNA485481:
#AAA     450     Abalone shriveling syndrome-associated virus :Abalone shriveling syndrome-associated virus:dsDNA:Haliotis diversicolor aquatilis:GCF_000882555.1::PRJNA485481:
#AAA     17      Abelson murine leukemia virus :Abelson murine leukemia virus:ssRNA-RT::GCF_000848265.1::PRJNA485481:
#AAA     98      Abisko virus :Abisko virus:RNA:Epirrita autumnata:GCF_002270725.1::PRJNA485481:

my $CUB=shift;
my %name;
my %test;

open(IN2,$CUB) or die ("can not open $CUB\n");
while (my $line=<IN2>) {
	chomp $line;
	my @mid1=split /\t/,$line;
	my @mid=split /\s+:|:/,$mid1[-1];
	unless (exists $name{$mid[4]} && length($name{$mid[4]})<length($mid[-1])) {
		$name{$mid[4]}=$mid1[-1];
	}
	$test{"$mid[4]\t$mid1[0]"}+=$mid1[1];
}


foreach my $name(keys %test) {
	my @mid=split /\t/,$name;
	print "$mid[1]\t$test{$name}\t$name{$mid[0]}\n";
}

