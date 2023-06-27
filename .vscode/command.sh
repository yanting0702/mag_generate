#!/bin/bash

main_dir=/home/yanting/DATA/Mangrove_samples_YX/mag_generate

# run snakemake
# create a folder for fastp_report
mkdir $main_dir/fastp_report


# conda install
conda env export --name fastp \
    --file $main_dir/envs/fastp.yml

# snakmake
conda activate snakemake

# save snakemake worklow (--rulegraph show only the rule; --dag show indiv jobs)
snakemake -s $main_dir/../scripts/mag_generate/.vscode/reads2mags.smk --rulegraph | dot -Tpdf > $main_dir/reads2mags.pdf

# dry-run
snakemake -n -s $main_dir/../scripts/mag_generate/.vscode/reads2mags.smk --use-conda

# actual-run
snakemake --cores 90 -s $main_dir/../scripts/mag_generate/.vscode/reads2mags.smk --use-conda

# # using cluster (faster)
# cp $code_dir/$main_dir/../scripts/mag_generate/.vscode/reads2mags.smk $main_dir/snake_files
# cd $main_dir/snake_files

# snakemake -s $main_dir/../scripts/mag_generate/.vscode/reads2mags.smk --use-conda \
# --cluster 'qsub -q low -l ncpus={threads},mem={params.mem},walltime=48:10:00' \
# -j 1000 --latency-wait 120