#mkdir -p ./mim/
rm -r /gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj20250709/
species=Hsap    #Hsap, Mmus,  Rnor, Scer, Spom, Dmel, Drer, Ecol, Atha
control_condition=Ctrl
outdir=/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/skj20250709
prefix=normal
input="/gpfs2/shankj/CodonUsage2/Ensembl/0Lib/01tRNA_seq/sample.txt"

mimseq --species $species \
        --cluster-id 0.97 \
        --threads 10 \
        --min-cov 0.0  \ #--min-cov 0.005  \ for tRNA modification
        --max-mismatches 0.1 \
        --control-condition $control_condition \
        --out-dir $outdir \
        --name $prefix \
        --max-multi 4 \
        --remap \
        --remap-mismatches 0.075 \
        $input

