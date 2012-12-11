#!/usr/bin/env Rscript
#####################################################################
#
# T. cruzi GOSeq Analysis
# Keith Hughitt <khughitt@umd.edu>
# 2012/12/10
#
# The aim of this script is to begin to explore the functionality
# provided by GOSeq for T. Cruzi RNA-Seq data.
#
# To generate a mapping of T. Cruzi gene ids and GO terms, the file
# tryp2go.py can be used.
#
#####################################################################
library(goseq)

# Use DESeq to determine which genes are differentially expressed
library(DESeq)

# Infection conditions
condition <- c("4 hours after infection", 
               "6 hours after infection",
               "12 hours after infection",
               "20 hours after infection",
               "24 hours after infection",
               "48 hours after infection",
               "72 hours after infection")

# Normalize counts
counts <- read.table('../../data/human_counts_replicates.csv', header=TRUE, row.names=1, sep=",")
cds <- newCountDataSet(counts, condition)
cds <- estimateSizeFactors(cds)

#Variance estimation
cds <- estimateDispersions(cds)

# DE Analysis
result <- nbinomTest(cds, "before_tc_infection", "48hrs_after_tc_infection")

#Filter out significantly DE genes with fdr < 0.1
significant_genes <- result[result$padj < 0.1,]

# Read in TriTryp gene ID/GO term mapping
go.terms <- read.delim('output/TcruziEsmeraldo_GOTerms.tsv', row.names=NULL)

