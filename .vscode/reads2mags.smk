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
        "../../../mag_generate/hum_ref_index/ref_done"

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
        "../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz"
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
        "../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz"
    params:
        mem="5G",
        outdir="../../../mag_generate/download_hum_ref"
    threads: 1
    shell:
        """
        mkdir -p {params.outdir} &&
        wget -P {params.outdir} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz
        """

rule join_rna_and_cds:
    input: 
        cds="../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz",
        rna="../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz"
    output:
        cds_rna_genomics="../../../mag_generate/download_hum_ref/human_cds_rna_genomics.fna",
        rna_decompress=temp("../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna"),
        cds_decompress=temp("../../../mag_generate/download_hum_ref/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna")
    params:
        mem="5G"
    threads: 8
    conda:
        "../../../mag_generate/envs/pigz.yml"
    shell:
        """
        (pigz -d {input.cds} -p {threads} -k &&
        pigz -d {input.rna} -p {threads} -k) &&
        cat {output.rna_decompress} {output.cds_decompress} > {output.cds_rna_genomics}
        """

rule bwa2rm_hum_ref:
    input: 
        "../../../mag_generate/download_hum_ref/human_cds_rna_genomics.fna"
    output:
        touch("../../../mag_generate/hum_ref_index/ref_done")
    params:
        mem="5G"
    threads: 8
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    shell:
        """
        bwa-mem2 index -p cds_rna_genomics {input} &&
        mv cds_rna_genomics* ../../../mag_generate/hum_ref_index/
        """