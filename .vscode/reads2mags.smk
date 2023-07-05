# import pandas as pd

# run_raw = pd.read_csv('../../frag_recruit/metadata/selected_epi_metag.csv')
# run_ids = list(run_raw.run_accessions)

# dataset_ids = ["sunlit"]

#define_wildcards

reads_id, = glob_wildcards("../../../rawdata/{reads}_1.fq.gz")

# expand function: wait for all wildcards runs finished

rule all:
    input:
        expand("../../../mag_generate/spades_assembly_output/{reads}/spades_done", reads=reads_id)
        
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
        touch("../../../mag_generate/download_hum_ref/ref_done")
    params:
        mem="5G"
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    shell:
        """
        bwa-mem2 index {input}
        """
    
rule bwa2rm_hum_map:
    input: 
        paired_1="../../../mag_generate/fastp/{reads}_1_paired.fq.gz",
        paired_2="../../../mag_generate/fastp/{reads}_2_paired.fq.gz",
        bwa_index_done="../../../mag_generate/download_hum_ref/ref_done",
        cds_rna="../../../mag_generate/download_hum_ref/human_cds_rna_genomics.fna"
    output:
        temp("../../../mag_generate/hum_cds_rna_mapped/{reads}.bam")
    params:
        mem="20G"
    threads: 12
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    shell:
        """
        bwa-mem2 mem -t {threads} {input.cds_rna} {input.paired_1} {input.paired_2} | \
        samtools view -@ {threads} -Sb > {output}
        """

rule samtools2rm_hum_sort:
    input: 
        "../../../mag_generate/hum_cds_rna_mapped/{reads}.bam"
    output:
        "../../../mag_generate/hum_cds_rna_mapped_sorted/{reads}.bam"
    params:
        mem="15G"
    threads: 8
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    shell:
        """
        samtools sort -@ {threads} {input} -o {output}
        """

rule samtools2rm_hum_index:
    input: 
        "../../../mag_generate/hum_cds_rna_mapped_sorted/{reads}.bam"
    output:
        "../../../mag_generate/hum_cds_rna_mapped_sorted/{reads}.bam.bai"
    params:
        mem="15G"
    threads: 1
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    shell:
        """
        samtools index {input}
        """
        
rule rscript2rm_hum_bam_aln_filter:
    input:
        bam="../../../mag_generate/hum_cds_rna_mapped_sorted/{reads}.bam",
        bai="../../../mag_generate/hum_cds_rna_mapped_sorted/{reads}.bam.bai",
        rscript="bam_ani_filter.r"
    output:
        mapped="../../../mag_generate/hum_cds_rna_r_filter_mapped/{reads}.bam",
        unmapped="../../../mag_generate/hum_cds_rna_r_filter_unmapped/{reads}.bam"
    threads: 12
    params:
        mem="10G",
        aln_ani="0.99",
        aln_cov="0.9",
        aln_len="50"
    shell:
        """
        {input.rscript} \
        {input.bam} {params.aln_ani} {params.aln_cov} {params.aln_len} \
        {output.mapped} {output.unmapped}
        """
    
rule samtools_unmapped_bam_sort_name:
    input:
        "../../../mag_generate/hum_cds_rna_r_filter_unmapped/{reads}.bam"
    output:
        temp("../../../mag_generate/hum_cds_rna_unmapped_sort/{reads}.bam")
    threads: 12
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    params:
        mem="5g"
    shell:
        """
        samtools sort -n -@ {threads} {input} -o {output}
        """
    
rule samtools_mapped_bam_sort_name:
    input:
        "../../../mag_generate/hum_cds_rna_r_filter_mapped/{reads}.bam"
    output:
       temp("../../../mag_generate/hum_cds_rna_mapped_sort/{reads}.bam")
    threads: 12
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    params:
        mem="5g"
    shell:
        """
        samtools sort -n -@ {threads} {input} -o {output}
        """

#todo split bam into forward and reverse fq files
rule samtools_unmapped_bam2fq:
    input:
        "../../../mag_generate/hum_cds_rna_unmapped_sort/{reads}.bam"
    output:
        unmapped_forwards="../../../mag_generate/hum_cds_rna_unmapped_seq/{reads}_f.fq.gz",
        unmapped_reverse="../../../mag_generate/hum_cds_rna_unmapped_seq/{reads}_r.fq.gz"
    threads: 5
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    params: 
        mem="5g",
        compress="6"
    shell:
        """
        samtools fastq -1 {output.unmapped_forwards} -2 {output.unmapped_reverse} \
        -@ {threads} -c {params.compress} -f 4 {input}
        """

rule samtools_mapped_bam2fq:
    input:
        "../../../mag_generate/hum_cds_rna_mapped_sort/{reads}.bam"
    output:
        mapped_forwards="../../../mag_generate/hum_cds_rna_mapped_seq/{reads}_f.fq.gz",
        mapped_reverse="../../../mag_generate/hum_cds_rna_mapped_seq/{reads}_r.fq.gz"
    threads: 5
    conda:
        "../../../mag_generate/envs/bwamen2.yml"
    params: 
        mem="5g",
        compress="6"
    shell:
        """
        samtools fastq -1 {output.mapped_forwards} -2 {output.mapped_reverse} \
        -@ {threads} -c {params.compress} -F 4 {input} 
        """

rule metaspades2assembly:
    input:
        unmapped_f="../../../mag_generate/hum_cds_rna_unmapped_seq/{reads}_f.fq.gz",
        unmapped_r="../../../mag_generate/hum_cds_rna_unmapped_seq/{reads}_r.fq.gz"
    output:
        touch("../../../mag_generate/spades_assembly_output/{reads}/spades_done")
    threads: 32
    conda:
        "../../../mag_generate/envs/spades.yml"
    params:
        mem="250g"
    shell:
        """
        spades.py --meta --threads {threads} \
        --pe1-1 {input.unmapped_f} --pe1-2 {input.unmapped_r} \
        -o {output} \
        -k 21,33,55,77,99,127
        """
