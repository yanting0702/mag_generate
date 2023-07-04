#!/home/yanting/anaconda3/envs/mamba/envs/bam_R_filter/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(GenomicAlignments)

param <- ScanBamParam(
    what=c("qname","flag","rname","strand","pos","qwidth","mapq","cigar","mrnm","mpos","isize","seq","qual"),
    tag=c("NM","MD","MC","AS","XS")
    )

aln <- readGAlignments(args[1], param=param, use.names=TRUE)

all_qname <- unlist(
    scanBam(args[1], param=ScanBamParam(what="qname")),
        use.names=FALSE)

df <- aln %>%
    as_tibble() %>%
    select(qname, seqnames, qwidth, start, end, NM) %>%
    mutate(
        aln_len = end - start,
        aln_ani = (aln_len - NM) / aln_len,
        aln_cov = aln_len / qwidth
    )

satisfied_qname <-
    filter(df, aln_ani >= as.numeric(args[2]),
        aln_cov >= as.numeric(args[3]),
        aln_len >= as.numeric(args[4]))$qname

unsatisfied_qname <- setdiff(all_qname, satisfied_qname)

# filter bam file by read_ids

filter_factory <- function(want) {
    list(KeepQname = function(x) x$qname %in% want)
}

filter_satisfied <- FilterRules(filter_factory(satisfied_qname))
filter_unsatisfied <- FilterRules(filter_factory(unsatisfied_qname))

filterBam(args[1], args[5], filter=filter_satisfied,
    param=ScanBamParam(what="qname"))

filterBam(args[1], args[6], filter=filter_unsatisfied,
    param=ScanBamParam(what="qname"))
