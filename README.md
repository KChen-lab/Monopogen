# scPopGene
Single Cell Population Genetics and Association Analysis Toolkit

**scPopGene** is one python package for population analysis in single-cell studies, developed and maintained by [Ken chen's lab](https://sites.google.com/view/kchenlab/Home) in MDACC. `scPopMap` is developed to benefit the population-level association study for single cell studies. It can work on datasets generated from single cell RNA 10x 5', 10x 3', smartseq, single ATAC-seq technoloiges. 
It is composed of four modules: 
* Germinle variant identification from shallow 10x scRNA-seq or scATAC-seq profiles. Given the sparsity of single cell sequencing data, we leverage linkage disequilibrium (LD) from external reference panel(such as 1KG3, TopMed) to refine genotypes. 
* Poteintal Somatic variant identification. The candidated variants with high alternative allele supporeted in study sample are further classifed based on their allele frequency patten among cell clusters (cell type/ cell states). 
* Population ancestry analysis (TBD)
* Association study/ GWAS variant annotation in single-cell level (TBD). 


## 1. Dependencies
* python  (version >= 3.73)
* java (open JDK>=1.8.0)
* pandas>=1.2.3
* pysam>=0.16.0.1
* NumPy>=1.19.5
* sciPy>=1.6.3
* pillow>=8.2.0
## 2. Installation 
Right now scPopGene is avaiable on github, you can install it through github 

`git clone https://github.com/KChen-lab/scPopGene.git`  
`cd scPopGene`  
`pip install -e .`  

## 3. Usage 
You can type the following command to get the help information.

`python ./src/scPopGene.py  SCvarCall --help`

```
usage: scPopGene.py SCvarCall [-h] -b BAMFILE -c CHR [-o OUT] -r REFERENCE -p
                              IMPUTATION_PANEL [-d DEPTH_FILTER]
                              [-t ALT_RATIO] [-m MAX_MISMATCH]
                              [-s MAX_SOFTCLIPPED] -a APP_PATH -i CELL_CLUSTER

optional arguments:
  -h, --help            show this help message and exit
  -b BAMFILE, --bamFile BAMFILE
                        The bam file for the study sample, the bam file should be sorted (default: None)
  -c CHR, --chr CHR     The chromosome used for variant calling (default: None)
  -o OUT, --out OUT     The output director (default: None)
  -r REFERENCE, --reference REFERENCE
                        The human genome reference used for alignment(default: None)
  -p IMPUTATION_PANEL, --imputation-panel IMPUTATION_PANEL
                        The population-level variant panel for variant refinement such as 1000 Genome 3 (default: None)
  -d DEPTH_FILTER, --depth_filter DEPTH_FILTER
                        The sequencing depth filter for variants not overlapped with public database (default: 50)
  -t ALT_RATIO, --alt_ratio ALT_RATIO
                        The minina allele frequency for variants as potential somatic mutation (default: 0.1)
  -m MAX_MISMATCH, --max-mismatch MAX_MISMATCH
                        The maximal mismatch allowed in one reads for variant calling (default: 3)
  -s MAX_SOFTCLIPPED, --max-softClipped MAX_SOFTCLIPPED
                        The maximal soft-clipped allowed in one reads for variant calling (default: 1)
  -a APP_PATH, --app-path APP_PATH
                        The app library paths used in the tool (default: None)
  -i CELL_CLUSTER, --cell_cluster CELL_CLUSTER
                        The cell cluster csv file used for somatic variant calling (default: None)
  ```
  
  
## 4. Examples

Here we provide demo of variant calling  based on data provided in the `example/` folder, which include:
* chr20_2Mb.rh.filter.sort.bam (.bai)
  The bam file storing read alignment for one study sample. Current scPopGene does not support mulitple sample mode. 
* CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz
  The reference panel with over 2,500 samples in 1000 Genome database. Only variants located in chr20: 0-2Mb are extracted. 
* chr20_2Mb.hg38.fa (.fai)
  The genome reference used for aligment. Only seuqence in chr20:0-20Mb are extracted.
* cell_cluster.csv 
  The cell cluster information. In the csv file, the first column is the cellID and second column is the cluster (can be derived from Seurat)





