# import pandas as pd

# run_raw = pd.read_csv('../../frag_recruit/metadata/selected_epi_metag.csv')
# run_ids = list(run_raw.run_accessions)

# dataset_ids = ["sunlit"]

#define_wildcards

reads_id, = glob_wildcards("../../../rawdata/{reads}_1.fq.gz")

# expand function: wait for all wildcards runs finished

rule all:
    input:
        expand("../../../mag_generate/fastp/{reads}_1_paired.fq.gz", reads=reads_id),
        expand("../../../mag_generate/fastp/{reads}_2_paired.fq.gz", reads=reads_id),
        rna="../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz",
        cds="../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz"

rule fastp_qc:
    input: 
        reads_f="../../../rawdata/{reads}_1.fq.gz",
        reads_r="../../../rawdata/{reads}_2.fq.gz"
    output:
        paired_1="../../../mag_generate/fastp/{reads}_1_paired.fq.gz",
        paired_2="../../../mag_generate/fastp/{reads}_2_paired.fq.gz",
        unpaired_1=temp("../../../mag_generate/fastp/{reads}_1_unpaired.fq.gz"),
        unpaired_2=temp("../../../mag_generate/fastp/{reads}_2_unpaired.fq.gz")
    params:
        mem="20G"
    threads: 8
    conda:
        "../../../mag_generate/envs/fastp.yml"
    shell:
        """
        fastp -i {input.reads_f} \
            -I {input.reads_r} \
            --out1 {output.paired_1} \
            --out2 {output.paired_2} \
            --unpaired1 {output.unpaired_1} \
            --unpaired2 {output.unpaired_2} \
            --compression 6  \
            --detect_adapter_for_pe \
            --thread {threads} \
            -q 20 -u 20 -g -c -W 5 -3 -l 50 \
            -j ../../../mag_generate/fastp_report/{wildcards.reads} \
            -h ../../../mag_generate/fastp_report/{wildcards.reads}
        """

rule download_hum_cds:
    output:
        cds="../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz"
    params:
        mem="5G",
        outdir="../../../mag_generate/download_hum_ref"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir} &&
        wget -P {params.outdir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz
        """

rule download_hum_rna:
    output:
        rna="../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz"
    params:
        mem="5G",
        outdir="../../../mag_generate/download_hum_ref"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir} &&
        wget -P {params.outdir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz
        """