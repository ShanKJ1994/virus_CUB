
#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;


#ENSMUST00000124102.8 cds chromosome:GRCm39:16:43929986:43959381:-1 gene:ENSMUSG00000052459.14 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:Atp6v1a description:ATPase, H+ transporting, lysosomal V1 subunit A [Source:MGI Symbol;Acc:MGI:1201780]

my $file=shift;#protein.fasta
$/="\n>";
my %seq;
my %cds;
open(IN,"gunzip -c $file |" ) or die ("can not open $file\n");
while (my $line=<IN>) {
	chomp $line;
	$line=~s/\>//;
	my @seq=split /\n/,$line;
	my @gene=split /\s+/,$seq[0];
	$gene[3]=~s/gene://;
	my $seq;
	for (my $i=1;$i<=$#seq ;$i++) {
		$seq.=$seq[$i];
	}
	unless (exists $seq{$gene[3]}) {
		$seq{$gene[3]}=$seq;
		$cds{$gene[3]}=$gene[0];
		
	}else{
		if (length($seq{$gene[3]})> length($seq)) {
			$seq{$gene[3]}=$seq;
			$cds{$gene[3]}=$gene[0];
		}
	}
}
close IN;

foreach my $gene(keys %seq) {
	print ">$gene:$cds{$gene}\n".$seq{$gene}."\n";
}