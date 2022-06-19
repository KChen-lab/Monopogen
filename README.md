# scPopGene
Single Cell Population Genetics and Association Analysis Toolkit

**scPopGene** is one python package for population analysis in single-cell studies, developed and maintained by [Ken chen's lab](https://sites.google.com/view/kchenlab/Home) in MDACC. `scPopMap` is developed to benefit the population-level association study for single cell studies. It can work on datasets generated from single cell RNA 10x 5', 10x 3', smartseq, single ATAC-seq technoloiges. 
It is composed of four modules: 
* Germinle variant identification from shallow 10x scRNA-seq or scATAC-seq profiles. Given the sparsity of single cell sequencing data, we leverage linkage disequilibrium (LD) from external reference panel(such as 1KG3, TopMed) to refine genotypes. 
* Poteintal Somatic variant identification. The candidated variants with high alternative allele supporeted in study sample are further classifed based on their allele frequency patten among cell clusters (cell type/ cell states). 
* Population ancestry analysis (TBD)
* Association study/ GWAS variant annotation in single-cell level (TBD). 

