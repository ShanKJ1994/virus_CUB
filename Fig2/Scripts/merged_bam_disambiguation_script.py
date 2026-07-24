from Bio import SeqIO
import pysam
import datetime
import os
import sys

virus_genome_fasta_path = sys.argv[1]
input_dir = sys.argv[2]
output_dir = sys.argv[3]

print(f'{datetime.datetime.now()}: Task started.', flush=True)

input_bam_files = os.listdir(input_dir)
input_bam_files = [f for f in input_bam_files if f.endswith('.bam')]
os.makedirs(output_dir, exist_ok=True)
virus_chr = set([record.id for record in SeqIO.parse(virus_genome_fasta_path, 'fasta')])

for input_bam in input_bam_files:
    print(f'{datetime.datetime.now()}: Started processing {input_bam}', flush=True)
    output_bam = input_bam.replace('.bam', '.disambiguated.bam')
    input_bam_path = os.path.join(input_dir, input_bam)
    output_bam_path = os.path.join(output_dir, output_bam)
    samfile = pysam.AlignmentFile(input_bam_path, "rb")
    virus_read_set = set()
    host_read_set = set()

    samfile_records = []
    for read in samfile:
        samfile_records.append(read)
        if read.reference_name in virus_chr:
            virus_read_set.add(read.query_name)
        else:
            host_read_set.add(read.query_name)

    ambiguous_read_set = host_read_set.intersection(virus_read_set)
    print(f'{datetime.datetime.now()}: Found {len(ambiguous_read_set)} reads aligned to both host and virus genome.', flush=True)

    out_bam = pysam.AlignmentFile(output_bam_path, "wb", header=samfile.header)
    removed_alignment_num = 0
    for read in samfile_records:
        if read.query_name not in ambiguous_read_set:
            out_bam.write(read)
        else:
            removed_alignment_num = removed_alignment_num + 1

    samfile.close()
    out_bam.close()

    print(f'{datetime.datetime.now()}: Finished processing {input_bam}, a total of {removed_alignment_num} alignments are removed.\n', flush=True)

print(f'{datetime.datetime.now()}: Finished all bam files.', flush=True)