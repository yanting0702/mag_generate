#!/bin/bash

main_dir=/home/yanting/DATA/Mangrove_samples_YX/mag_generate
snakemake_script_dir=$main_dir/../scripts/mag_generate/.vscode

# run snakemake
# create a folder for fastp_report
mkdir $main_dir/fastp_report

# produce conda enviroment information for snakemake
conda env export --name fastp \
    --file $main_dir/envs/fastp.yml

conda env export --name pigz \
    --file $main_dir/envs/pigz.yml

conda env export --name bwamen2 \
    --file $main_dir/envs/bwamen2.yml

conda env export --name spades \
    --file $main_dir/envs/spades.yml

conda env export --name megahit \
    --file $main_dir/envs/megahit.yml

# activate r srcipt
chmod +x $snakemake_script_dir/bam_ani_filter.r

# snakmake
#! the relative path in smk file is relative to the current working directory (i.e., pwd),
#! not the path where the smk file is located

conda activate snakemake

#! here we entered into the smk file directory, so the pwd = smk file path
cd $snakemake_script_dir

# save snakemake worklow (--rulegraph show only the rule; --dag show indiv jobs)
snakemake -s reads2mags.smk --rulegraph | dot -Tpdf > $main_dir/reads2mags.pdf

# dry-run
snakemake -n -s reads2mags.smk --use-conda

# actual-run
snakemake --cores 120 -s reads2mags.smk --use-conda

# # using cluster (faster)
# cp $code_dir/reads2mags.smk $main_dir/snake_files
# cd $main_dir/snake_files

# snakemake -s reads2mags.smk --use-conda \
# --cluster 'qsub -q low -l ncpus={threads},mem={params.mem},walltime=48:10:00' \
# -j 1000 --latency-wait 120


