#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
##########################################################
#read in file which contain the cds pos -> codon in canonical transcript of one gene

#ENST00000322886 0       ATG GAG
#ENST00000322886 1       ATG GAG
#ENST00000322886 2       ATG GAG

my %pos;
my $file1=shift;
my @codonfile=split /,/,$file1;
foreach my $file(@codonfile) {
	open(IN, $file) or die ("can not open $file\n");
	while (my $line=<IN>) {
		chomp $line;
		my @mid=split /\t/,$line;
		$mid[2]=uc($mid[2]);
		$mid[3]=uc($mid[3]);
		if ($mid[2] ne "TGA" && $mid[2] ne "TAA" && $mid[2] ne "TAG" &&
			$mid[3] ne "TGA" && $mid[3] ne "TAA" && $mid[3] ne "TAG") {
			$mid[0]=~s/_transcript//g;
			$mid[0]=~s/transcript://g;
			#my @mid2=split /:/,$mid[0];
			$pos{$mid[0]}{$mid[1]}=("$mid[2]\t$mid[3]");
			#print "$mid2[1]\t$mid[1]\t$pos{$mid2[1]}{$mid[1]}\n";
		}
	}
	close IN;
}
##########################################################
#mRNA read in CDS pos read coverage

# ENST00000472194 4314    0.000000,0.000000,0.000000,0.000000,0.000000,

my $mRNA=shift;
my %mRNA;
open(IN, $mRNA) or die ("can not open $mRNA\n");
while (my $line=<IN>) {
	chomp $line;
	my @mid=split /\t/,$line;
	#$mid[0]=ENST00000472194
	my @array=split /,/,$mid[2];
	# 使用grep函数找出所有大于0的元素，并计算它们的数量
	my $count = scalar grep { $_ > 0 } @array;
	if ($count/$mid[1]>=0.5 && $mid[1]>=150*3) {#$count/$mid[1]>=0.5 && 
		$mid[0]=~s/_transcript//g;
		$mid[0]=~s/transcript://g;
		$mRNA{$mid[0]}=\@array;
		#print "$mid[0]\tmRNA\n";
	}
}
close IN;

##########################################################
#ribo read in CDS pos read coverage

# ENST00000472194 4314    0.000000,0.000000,0.000000,0.000000,0.000000,

my $Ribo1=shift;
my @Ribofile=split /,/,$Ribo1;
my %Ribo;
foreach my $Ribo(@Ribofile) {
	open(IN, $Ribo) or die ("can not open $Ribo\n");
	while (my $line=<IN>) {
		chomp $line;
		my @mid=split /\t/,$line;
		#$mid[0]=ENST00000472194
		my @array=split /,/,$mid[2];
		# 使用grep函数找出所有大于0的元素，并计算它们的数量
		#my $count = scalar grep { $_ > 0 } @array;
		if ($mid[1]>=150*3) {
			$mid[0]=~s/transcript://g;
			$mid[0]=~s/_transcript//g;
			if (exists $Ribo{$mid[0]}) {
				my @array2=@{$Ribo{$mid[0]}};
				my @result=();
				for (my $i=0;$i<=$#array;$i++) {
					$result[$i] = $array[$i] + $array2[$i];
				}
				$Ribo{$mid[0]}=\@result;
			}else{
				$Ribo{$mid[0]}=\@array;
			#print "$mid[0]\tRibo\n";
			}
		}
	}
	close IN;
}

$mRNA =~ s/mRNA_//g;
$mRNA =~ s/.disambiguated_transcript//g;
#open (OUT,">verbose.$mRNA")||die;
open (OUT0,">../frame/frame0.RiboCodon.$mRNA")||die;
open (OUT1,">../frame/frame012.RiboCodon.$mRNA")||die;
open (OUT2,">../frame/frame201.RiboCodon.$mRNA")||die;
#print OUT "Tr_ID\tPos\tPsite\tAsite\tmRNA\tRibo\n";
print OUT0 "Tr_ID\tPos\tPsite\tAsite\tmRNA\tRibo\n";
print OUT1 "Tr_ID\tPos\tPsite\tAsite\tmRNA\tRibo\n";
print OUT2 "Tr_ID\tPos\tPsite\tAsite\tmRNA\tRibo\n";

my %out012;
my %out201;
my %out0;
foreach my $tr(keys %Ribo) {
	if (exists $mRNA{$tr} && exists $pos{$tr}) {
		my @mRNA = @{$mRNA{$tr}}; 
		my @Ribo = @{$Ribo{$tr}}; 
		#for (my $i=0;$i<$#mRNA;$i++) {
			#my $mRNA=($mRNA[$i]);
			#my $Ribo=($Ribo[$i]);
			#if (exists $pos{$tr}{$i}) {
				#print OUT "$tr\t$i\t$pos{$tr}{$i}\t$mRNA\t$Ribo\n";
			#}
		#}
		#前15个codon和后5个codon去掉
		for (my $i=15*3;$i<$#mRNA-5*3;$i+=3) {
			if ($Ribo[$i]>0) {
				$out0{$tr}+=1/($#Ribo/3);
			}
			if (($Ribo[$i]+ $Ribo[$i+1] + $Ribo[$i+2])>0) {
				$out012{$tr}+=1/($#Ribo/3);
			}
			if (($Ribo[$i]+ $Ribo[$i+1] + $Ribo[$i-1])>0) {
				$out201{$tr}+=1/($#Ribo/3);
			}
		}
	}
}
#close OUT;

foreach my $tr(keys %Ribo) {
	if (exists $mRNA{$tr} && exists $pos{$tr}) {
		my @mRNA = @{$mRNA{$tr}}; 
		my @Ribo = @{$Ribo{$tr}}; 
		#前15个codon和后5个codon去掉
		for (my $i=15*3;$i<$#mRNA-5*3;$i+=3) {
			if ($mRNA[$i]>0 && exists $pos{$tr}{$i} && $out0{$tr}>=0.5) {
				print OUT0 "$tr\t$i\t$pos{$tr}{$i}\t$mRNA[$i]\t$Ribo[$i]\n";
			}
			if ($mRNA[$i]>0 && $mRNA[$i+1]>0 && $mRNA[$i+2]>0 && $out012{$tr}>=0.5) {
				if (exists $pos{$tr}{$i}) {
				my $mRNA=($mRNA[$i]+$mRNA[$i+1]+$mRNA[$i+2]);
				my $Ribo=($Ribo[$i]+$Ribo[$i+1]+$Ribo[$i+2]);
				print OUT1 "$tr\t$i\t$pos{$tr}{$i}\t$mRNA\t$Ribo\n";
				}
			}
			if ($mRNA[$i]>0 && $mRNA[$i+1]>0 && $mRNA[$i-1]>0 && $out201{$tr}>=0.5) {
				if (exists $pos{$tr}{$i}) {
					my $mRNA=($mRNA[$i]+$mRNA[$i+1]+$mRNA[$i-1]);
					my $Ribo=($Ribo[$i]+$Ribo[$i+1]+$Ribo[$i-1]);
					print OUT2 "$tr\t$i\t$pos{$tr}{$i}\t$mRNA\t$Ribo\n";
				}
			}
		}
	}
}
