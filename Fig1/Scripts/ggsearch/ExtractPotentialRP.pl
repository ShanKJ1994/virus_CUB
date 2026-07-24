#!/usr/bin/perl -w
use strict;
use warnings;
#use lib 'D:/e-books';
#use BeginPerlBioinfo;
##############################################################################

#0          1   2        3    4         5         6   7       8   9     10     11
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#read in blast file
#AT3G52590.1__pep__chromosome:TAIR10:3:19505379:19506932:1__gene:AT3G52590__transcript:AT3G52590.1__gene_biotype:protein_coding__transcript_biotype:protein_coding__gene_symbol:UBQ1__description:ubiquitin__extension__protein__1__[Source:NCBI__gene__(formerly__Entrezgene);Acc:824425]	ENSMUSG00000090137.9:ENSMUST00000129909.8	92.19	128	10	0	1	128	1	128	1.3e-309	84.9
#AT3G52590.1__pep__chromosome:TAIR10:3:19505379:19506932:1__gene:AT3G52590__transcript:AT3G52590.1__gene_biotype:protein_coding__transcript_biotype:protein_coding__gene_symbol:UBQ1__description:ubiquitin__extension__protein__1__[Source:NCBI__gene__(formerly__Entrezgene);Acc:824425]	ENSMUSG00000090137.9:ENSMUST00000081940.11	92.19	128	10	0	1	128	1	128	1.3e-309	84.9
#AT3G52590.1__pep__chromosome:TAIR10:3:19505379:19506932:1__gene:AT3G52590__transcript:AT3G52590.1__gene_biotype:protein_coding__transcript_biotype:protein_coding__gene_symbol:UBQ1__description:ubiquitin__extension__protein__1__[Source:NCBI__gene__(formerly__Entrezgene);Acc:824425]	ENSMUSG00000090137.9:ENSMUST00000135446.8	92.19	128	10	0	1	128	1	128	1.3e-309	84.9

my $file2=shift;
my %seq;
my %cds;
open(IN,$file2) or die ("can not open $file2\n");
while (my $line=<IN>) {
	my @mid=split /\s+/,$line;
	if ($mid[10]<0.01 && $mid[2]>=70 ) {
		my @gene=split /:/,$mid[1];
		unless (exists $cds{$gene[0]}) {
			$cds{$gene[0]}=$mid[1];
			$seq{$gene[0]}=$mid[2];
		}else{
			if (exists $cds{$gene[0]} && $seq{$gene[0]} < $mid[2]) {
				$cds{$gene[0]}=$mid[1];
				$seq{$gene[0]}=$mid[2];
			}
		}
	}
}
close IN;


foreach my $name(keys %cds) {
	print "$cds{$name}\n";
}
