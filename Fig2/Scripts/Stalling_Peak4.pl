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
my $file2=shift;#/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/00Final/UniqueMapped/Final.CDS.codon.pos.txt
my @file=split /,/,$file2;
foreach my $file (@file) {
open(IN, $file) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
	my @mid=split /\t/,$line;
	$mid[2]=uc($mid[2]);
	$mid[3]=uc($mid[3]);
	if ($mid[2] ne "TGA" && $mid[2] ne "TAA" && $mid[2] ne "TAG" &&
		$mid[3] ne "TGA" && $mid[3] ne "TAA" && $mid[3] ne "TAG") {
		$mid[0]=~s/_transcript//g;
		my $gene=$mid[0];
		if ($mid[0]=~/:/) {
			my @mid2=split /:/,$mid[0];
			$gene=$mid2[1];
		}
		
		$pos{$gene}{$mid[1]}=uc("$mid[2]\t$mid[3]");
		#print "$gene\t$mid[1]\t$pos{$gene}{$mid[1]}\n";
	}
}
close IN;
}
##########################################################
#mRNA read in CDS pos read coverage

# ENST00000472194 4314    0.000000,0.000000,0.000000,0.000000,0.000000,

my $mRNA_file=shift;
my @mRNA_file=split /,/,$mRNA_file;
my %mRNA;
foreach my $mRNA (@mRNA_file) {
	open(IN, $mRNA) or die ("can not open $mRNA\n");
	while (my $line=<IN>) {
		chomp $line;
		my @mid=split /\t/,$line;
		#$mid[0]=ENST00000472194
		my @array=split /,/,$mid[2];
		## 使用grep函数找出所有大于0的元素，并计算它们的数量
		#my $count = scalar grep { $_ > 0 } @array;
		if ($mid[1]>=150*3) {
			$mid[0]=~s/transcript://g;
			$mid[0]=~s/_transcript//g;
			if (exists $mRNA{$mid[0]}) {
				my @array2=@{$mRNA{$mid[0]}};
				my @result=();
				for (my $i=0;$i<=$#array;$i++) {
					$result[$i] = $array[$i] + $array2[$i];
				}
				# 使用grep函数找出所有大于等于5的元素，并计算它们的数量
				my $count = scalar grep { $_ >= 5 } @array;
				if ($count/$mid[1]>=0.5 && $mid[1]>=150*3) {#$count/$mid[1]>=0.5 && 
					$mid[0]=~s/_transcript//g;
					$mRNA{$mid[0]}=\@array;
					#print "$mid[0]\tmRNA\n";
				}
				$mRNA{$mid[0]}=\@result;
			}else{
				$mRNA{$mid[0]}=\@array;
			}
		}
		
	}
	close IN;
}

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

##########################################################

my $file_name=shift;
open (OUT0,">Stall4.frame0.$file_name.txt")||die;
open (OUT1,">Stall4.frame012.$file_name.txt")||die;
open (OUT2,">Stall4.frame201.$file_name.txt")||die;
print OUT0 "Tr_ID\tPos\tPsite\tAsite\tmRNA\tRibo\tmRNA_F\tRibo_F\tTr_RPF_Coverage\n";
print OUT1 "Tr_ID\tPos\tPsite\tAsite\tmRNA\tRibo\tmRNA_F\tRibo_F\tTr_RPF_Coverage\n";
print OUT2 "Tr_ID\tPos\tPsite\tAsite\tmRNA\tRibo\tmRNA_F\tRibo_F\tTr_RPF_Coverage\n";


my %out012;
my %out201;
my %out0;
foreach my $tr(keys %Ribo) {
	if (exists $mRNA{$tr} && exists $pos{$tr}) {
		my @mRNA = @{$mRNA{$tr}}; 
		my @Ribo = @{$Ribo{$tr}}; 
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
close OUT;

foreach my $tr(keys %Ribo) {
	if (exists $mRNA{$tr} && exists $pos{$tr}) {
		my @mRNA = @{$mRNA{$tr}}; 
		my @Ribo = @{$Ribo{$tr}}; 
		#前15个codon和后5个codon去掉
		#print "$tr\t@mRNA\t@Ribo\n";
		for (my $i=45*3;$i<=$#mRNA-35*3;$i+=3) {
			#print "$tr\t$i\n";
			if ($mRNA[$i]>=5 && exists $pos{$tr}{$i} && exists $out0{$tr}) {# && $out0{$tr}>=0.5
				my $mRNA=$mRNA[$i];
				my $Ribo=$Ribo[$i];
				my $mRNA_F=0;
				my $Ribo_F=0;
				for (my $j=$i-30*3;$j<=$i+30*3;$j+=3) {
					$mRNA_F+=$mRNA[$j];
					$Ribo_F+=$Ribo[$j];
				}
				print OUT0 "$tr\t$i\t$pos{$tr}{$i}\t$mRNA\t$Ribo\t$mRNA_F\t$Ribo_F\t$out0{$tr}\n";
			}
			if ($mRNA[$i]>=5 && $mRNA[$i+1]>=5 && $mRNA[$i+2]>=5 && exists $out012{$tr} ) {#&& $out012{$tr}>=0.5
				if (exists $pos{$tr}{$i}) {
				my $mRNA=($mRNA[$i]+$mRNA[$i+1]+$mRNA[$i+2]);
				my $Ribo=($Ribo[$i]+$Ribo[$i+1]+$Ribo[$i+2]);
				my $mRNA_F=0;
				my $Ribo_F=0;
					for (my $j=$i-30*3;$j<=$i+30*3;$j+=3) {
						$mRNA_F+=($mRNA[$j]+$mRNA[$j+1]+$mRNA[$j+2]);
						$Ribo_F+=($Ribo[$j]+$Ribo[$j+1]+$Ribo[$j+2]);
					}
				print OUT1 "$tr\t$i\t$pos{$tr}{$i}\t$mRNA\t$Ribo\t$mRNA_F\t$Ribo_F\t$out012{$tr}\n";
				}
			}
			if ($mRNA[$i]>=5 && $mRNA[$i+1]>=5 && $mRNA[$i+2]>=5 && exists $out201{$tr} ) {#&& $out201{$tr}>=0.5
				if (exists $pos{$tr}{$i}) {
					my $mRNA=($mRNA[$i]+$mRNA[$i+1]+$mRNA[$i+2]);
					my $Ribo=($Ribo[$i]+$Ribo[$i+1]+$Ribo[$i-1]);
					my $mRNA_F=0;
					my $Ribo_F=0;
					for (my $j=$i-30*3;$j<=$i+30*3;$j+=3) {
						$mRNA_F+=($mRNA[$j]+$mRNA[$j+1]+$mRNA[$j+2]);
						$Ribo_F+=($Ribo[$j]+$Ribo[$j+1]+$Ribo[$j-1]);
					}
					print OUT2 "$tr\t$i\t$pos{$tr}{$i}\t$mRNA\t$Ribo\t$mRNA_F\t$Ribo_F\t$out201{$tr}\n";
				}
			}
		}
	}
}
