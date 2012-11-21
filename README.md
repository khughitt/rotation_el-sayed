Rotation Project: T. cruzi RNA-Seq Analysis
===========================================

Overview
--------
This repository includes code related to work I did as part of a first-year lab
rotation in the [El-Sayed laboratory](http://www.najibelsayed.org/).

The goal of the project was to get exposed to some of the kinds of problems
that are being worked on in the lab, in particular, with respect to 
[Trypanosoma Cruzi](http://en.wikipedia.org/wiki/Trypanosoma_cruzi)
transcriptomics.

...(find at end of lab rotation)

RNA-Seq Time Series Analysis
----------------------------

The first problem I worked on was the clustering and visualization of T. cruzi
infection RNA-Seq time series data. The started with DESeq counttables for 
genes mapped to either a Human or T. cruzi reference genome.

![screenshot](https://raw.github.com/khughitt/rotation_el-sayed/master/extra/screenshot.png)

Files:
* **rna-seq_clustering.py**
* **notebooks/T_Cruzi_Gene_Expression_Time_Series_Analysis.ipynb**

Although a significant amount of reserach over the past few years has focused
on the normalizing RNA-seq data, e.g. [Sun & Zhu (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22914217),
for the purposes of my rotation project I decided not to spend too much time
up front choosing the best approach and instead simply converted the DESeq
counts to relative fractions (% total for a given sample).

To group genes with similar patterns of expression across time, [k-means clustering](http://en.wikipedia.org/wiki/K-means_clustering)
was used. A few different arbitrary values of k were tested.

Finally, the results of the clustering were visualized using a [parallel coordinates plot](http://en.wikipedia.org/wiki/Parallel_coordinates),
and genes which appeared to show significant patterns of differential expression
were pulled out and mapped to functional annotations from [Ensembl](http://useast.ensembl.org/Homo_sapiens/Info/Index)
(human) and [TriTrypDB](http://tritrypdb.org/tritrypdb/) (T. cruzi).


