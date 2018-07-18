# Benchmarking pipeline for ES cell data

Prerequisites:

-virtualenv

-github account

-reference genome gtf and fasta files. The Mus musculus Ensembl release 89 genome was used in this study with ERCC sequences appended (see https://tools.thermofisher.com/content/sfs/manuals/cms_095048.txt).

-java version 1.8

-R version 3.4.4


To run the pipeline:

Execute ./wrapper.sh path/to/java path/to/ref/fasta path/to/ref/gtf

In practice it is unlikely that your machine will have the resources to run the entire pipeline in one go, so you will probably need to split up the wrapper script and run it in bits.

As part of the pipeline, quality control steps are automatically carried out. For reference, these are the statistics used to filter the raw data:

| Statistic | Name of statistic in table | Threshold |
-------------|--------|---------
|No. uniquely mapping reads|Unique    | >8000000 |
|No. of non-uniquely mapping reads|NonUnique|>350000|
|No. alignments|NumAlign|>8200000|
|No. of reads|NumReads|>4000000|

And the simulated data:

| Statistic | Name of statistic in table | Threshold |
-------------|--------|---------
|No. of non-uniquely mapping reads|NonUnique|>2500000|
|No. alignments|NumAlign|>32000000 or <4000000|
|No. of reads|NumReads|>12000000 or <3500000|

In addition, the scater package was used to filter cells in which more than 10% of reads mapped to mitochondrial genes in both the raw and simulated data.
