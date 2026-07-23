#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;

##############################################################################

# Plant

##############################################################################

my $Plant="/gpfs2/shankj/CodonUsage2/Ensembl/Plant_class.txt";
my %host;
open(IN2,$Plant) or die ("can not open $Plant\n");
while (my $line=<IN2>) {
	chomp $line;
	my @mid=split /\t/,$line;
	my @shu=split /_/,$mid[0];
	$host{$shu[0]}=$mid[1];
}
close IN2;

##############################################################################

# Animal

##############################################################################

#Spiny chromis	Acanthochromis polyacanthus	fish	Acanthochromis_polyacanthus
#Eurasian sparrowhawk	Accipiter nisus	bird	Accipiter_nisus
#Giant panda	Ailuropoda melanoleuca	mammal	Ailuropoda_melanoleuca

my $Animal="/gpfs2/shankj/CodonUsage2/Ensembl/Animal_class.txt";
open(IN2,$Animal) or die ("can not open $Animal\n");
while (my $line=<IN2>) {
	chomp $line;
	my @mid=split /\t/,$line;
	my @shu=split /\s+/,$mid[1];
	$host{$shu[0]}=$mid[2];
}
close IN2;


##############################################################################

# virus-> host

##############################################################################

#Organism_Name   Species Host
#Frijoles virus  Frijoles phlebovirus    Lutzomyia
#Bird's-foot trefoil nucleorhabdovirus   Bird's-foot trefoil nucleorhabdovirus   Lotus corniculatus
#Erwinia phage Fifi44    Erwinia phage Fifi44    Erwinia amylovora
#Escherichia phage PO103-1       Escherichia phage PO103-1       Escherichia coli
#Escherichia phage MLP1  Escherichia phage MLP1  Escherichia coli CFT073
#Yersinia phage PYps3T   Yersinia phage PYps3T   Yersinia pseudotuberculosis
#Yersinia phage PYps23T  Yersinia phage PYps23T  Yersinia pseudotuberculosis


my $file="/gpfs2/shankj/CodonUsage2/Ensembl/Virus.Host.txt";
my %virus;
open(IN2,$file) or die ("can not open $file\n");
while (my $line=<IN2>) {
	chomp $line;
	my @mid=split /\t/,$line;
	if ($mid[0] =~ /Influenza B virus/) {
        $mid[0] = "Influenza B virus";
    }
    elsif ($mid[0] =~ /Influenza C virus/) {
        $mid[0] = "Influenza C virus";
    }
    elsif ($mid[0] =~ /Influenza A virus.*?\((H\d+N\d+)\)/) {
        $mid[0] = $1;
    }
	if ($mid[2]=~/[A-Z]/i) {
		my @shu=split /\s+/,$mid[2];
		if (exists $host{$shu[0]}) {
			$virus{"$mid[0]:$mid[1]"}{$shu[0]}=$host{$shu[0]};
		}
	}
}
close IN2;


##############################################################################

# virus CUB-> host

##############################################################################

#AAA     19      Abaca bunchy top virus :Abaca bunchy top virus:ssDNA(+):Musa sp.:GCF_000872625.1::PRJNA485481:
#AAA     2197    Abalone herpesvirus Victoria/AUS/2009 :Aurivirus haliotidmalaco1:dsDNA::GCF_000900375.1::PRJNA485481:
#AAA     450     Abalone shriveling syndrome-associated virus :Abalone shriveling syndrome-associated virus:dsDNA:Haliotis diversicolor aquatilis:GCF_000882555.1::PRJNA485481:
#AAA     17      Abelson murine leukemia virus :Abelson murine leukemia virus:ssRNA-RT::GCF_000848265.1::PRJNA485481:
#AAA     98      Abisko virus :Abisko virus:RNA:Epirrita autumnata:GCF_002270725.1::PRJNA485481:

my $CUB="/gpfs2/shankj/CodonUsage2/Ensembl/VirusCDS/AllVirus_CUB.Final.txt";
my %test;
open(IN2,$CUB) or die ("can not open $CUB\n");
while (my $line=<IN2>) {
	chomp $line;
	my @mid1=split /\t/,$line;
	my @mid=split /\s+:|:/,$mid1[-1];
	unless (exists $test{$mid[4]} && length($test{$mid[4]})>length($mid[-1])) {
		$test{$mid[4]}=$mid1[-1];
	}
}

open(IN2,$CUB) or die ("can not open $CUB\n");
while (my $line=<IN2>) {
	chomp $line;
	my @mid=split /\t/,$line;
	my @shu=split /\s+:|:/,$mid[-1];
	if ($shu[0] =~ /Influenza B virus/) {
        $shu[0] = "Influenza B virus";
    }
    elsif ($shu[0] =~ /Influenza C virus/) {
        $shu[0] = "Influenza C virus";
    }
    elsif ($shu[0] =~ /Influenza A virus.*?\((H\d+N\d+)\)/) {
        $shu[0] = $1;
    }
	if (exists $virus{"$shu[0]:$shu[1]"}) {
		foreach my $shu(keys %{$virus{"$shu[0]:$shu[1]"}}) {
			print $virus{"$shu[0]:$shu[1]"}{$shu};
			print "\t$shu\t$test{$shu[4]}\t$shu[1]\t$mid[0]\t$mid[1]\n";
		}
	}
	
}
close IN2;




