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
library(DESeq)
library(goseq)

# To begin, use DESeq to determine which genes are differentially expressed

# input file
input.dir <- '../data/DEseq_comparisons/hpgl0062vs64vs66vs68vs54vs59vs70vs56vs60vs72vs58vs61_tophatv2.0.3_deseq_esmer_oneloci'
input.file <- 'hpgl0062vs64vs66vs68vs54vs59vs70vs56vs60vs72vs58vs61_tophatv2.0.3_deseq_esmer_oneloci.counttable.sorted'

# Read in count count table
file.path <- paste(input.dir, input.file, sep='/')
count.table <- read.table(file.path, header=TRUE, row.names=1)

# Exclude ambiguous and no_feature rows
count.table <- count.table[!rownames(count.table) %in% c('ambiguous', 'no_feature'),]

# Infection conditions
condition <- c("Prior to infection",
               "4 hours after infection",
               "6 hours after infection",
               "12 hours after infection",
               "20 hours after infection",
               "20 hours after infection",
               "24 hours after infection",
               "48 hours after infection",
               "48 hours after infection",
               "48 hours after infection",
               "72 hours after infection",
               "72 hours after infection")

cds <- newCountDataSet(count.table, condition)

# Estimate gene abundance normalization factor
cds <- estimateSizeFactors(cds)

# Variance estimation
cds <- estimateDispersions(cds)

# DE Analysis
result <- nbinomTest(cds, "Prior to infection", "48 hours after infection")

# Filter out significantly DE genes with fdr < 0.1
significant_genes <- result[result$padj < 0.1,]

# Read in TriTryp gene ID/GO term mapping
go.terms <- read.delim('output/TcruziEsmeraldo_GOTerms.tsv', row.names=NULL)
