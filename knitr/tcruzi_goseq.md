Table of Contents
=================

1. Overview
2. Data Preparation
3. Differential Expression Analysis (DESeq)
4  GO Enrichment Analysis (GOSeq)
5. Meta

Overview
========

**STATUS**: In Progress (2012/12/11)

The goal of this document is to explore the use of [DESeq](http://www-huber.embl.de/users/anders/DESeq/)
and [GOSeq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html)
to look for functional enrichment in [T. Cruzi](http://en.wikipedia.org/wiki/Trypanosoma_cruzi)
RNA-Seq data.

Data Preparation
================

The first step necessary is to read in both the RNA-Seq data, and also a table
relating T. Cruzi gene id's and GO terms. GOSeq supports automatically mapping
gene id's to the necessary GO terms and gene lengths for certain organisms and
gene identifiers (e.g. human/Ensembl), however, it does not currently support
[TriTrypDB](http://tritrypdb.org/tritrypdb/) id's.

GOSeq can still be used for the enrichment analsys, however, it must be manually
provided with the necessary mapping and length info.

Python was used (see tryp2go.py) to process a data file from TriTrypDB which
contains all of the neccessary information: `TriTrypDB-4.2_TcruziEsmeraldo-LikeGene.txt`.
The output of this is a tab-delimited file (`TcruziEsmeraldo_GOTerms.tsv`) where
each row includes a gene id, gene length, GO id, GO ontology, GO term name, source
and evidence code.

**Future Work**: In addition to using GO terms listed in TriTrypDB, another
possibility might be to use [pfam2go](http://www.geneontology.org/external2go/pfam2go)
to scan the T. Cruzi transcriptome at the protein family level to look for new
GO term associations.

Differential Expression Analysis (DESeq)
========================================

To begin, let's read in the RNA-Seq count data generated from HTSeq.


```r

# Input files
input.dir <- "../../data/DEseq_comparisons/hpgl0062vs64vs66vs68vs54vs59vs70vs56vs60vs72vs58vs61_tophatv2.0.3_deseq_esmer_oneloci"
input.file <- "hpgl0062vs64vs66vs68vs54vs59vs70vs56vs60vs72vs58vs61_tophatv2.0.3_deseq_esmer_oneloci.counttable.sorted"

# Read in count count table
file.path <- paste(input.dir, input.file, sep = "/")
count.table <- read.table(file.path, header = TRUE, row.names = 1)

# Exclude ambiguous and no_feature rows
count.table <- count.table[!rownames(count.table) %in% c("ambiguous", "no_feature"), 
    ]
```



```r
library(DESeq)
```



```r
# Infection conditions
condition <- c("Prior to infection", "4 hours after infection", "6 hours after infection", 
    "12 hours after infection", "20 hours after infection", "20 hours after infection", 
    "24 hours after infection", "48 hours after infection", "48 hours after infection", 
    "48 hours after infection", "72 hours after infection", "72 hours after infection")

cds <- newCountDataSet(count.table, condition)
```


Next we will normalize for transcript abundance and estimate sample variance.


```r
# Estimate gene abundance normalization factor
cds <- estimateSizeFactors(cds)

# Variance estimation
cds <- estimateDispersions(cds)
```

```
## Warning: glm.fit: algorithm did not converge
```

```
## Warning: glm.fit: algorithm did not converge
```

```
## Warning: glm.fit: algorithm did not converge
```

```
## Warning: Dispersion fit did not converge.
```


Now we are ready to perform look for DE.


```r
# DE Analysis
result <- nbinomTest(cds, "Prior to infection", "48 hours after infection")

# Filter out significantly DE genes with fdr < 0.1
significant_genes <- result[result$padj < 0.1, ]
```


GO Enrichment Analysis (GOSeq)
==============================

```r
library(goseq)
```


Meta
====

System Information
------------------

```r
sessionInfo()
```

```
## R version 2.15.2 (2012-10-26)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C             
##  [3] LC_TIME=en_US.utf8        LC_COLLATE=en_US.utf8    
##  [5] LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8   
##  [7] LC_PAPER=C                LC_NAME=C                
##  [9] LC_ADDRESS=C              LC_TELEPHONE=C           
## [11] LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C      
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] goseq_1.10.0            geneLenDataBase_0.99.10
##  [3] BiasedUrn_1.04          DESeq_1.10.1           
##  [5] lattice_0.20-10         locfit_1.5-8           
##  [7] Biobase_2.18.0          BiocGenerics_0.4.0     
##  [9] knitr_0.9               colorout_1.0-0         
## 
## loaded via a namespace (and not attached):
##  [1] annotate_1.36.0        AnnotationDbi_1.20.3   biomaRt_2.14.0        
##  [4] Biostrings_2.26.2      bitops_1.0-5           BSgenome_1.26.1       
##  [7] DBI_0.2-5              digest_0.6.0           evaluate_0.4.3        
## [10] formatR_0.7            genefilter_1.40.0      geneplotter_1.36.0    
## [13] GenomicFeatures_1.10.1 GenomicRanges_1.10.5   grid_2.15.2           
## [16] IRanges_1.16.4         Matrix_1.0-10          mgcv_1.7-22           
## [19] nlme_3.1-106           parallel_2.15.2        RColorBrewer_1.0-5    
## [22] RCurl_1.95-3           Rsamtools_1.10.2       RSQLite_0.11.2        
## [25] rtracklayer_1.18.1     splines_2.15.2         stats4_2.15.2         
## [28] stringr_0.6.2          survival_2.36-14       tools_2.15.2          
## [31] XML_3.95-0.1           xtable_1.7-0           zlibbioc_1.4.0
```


