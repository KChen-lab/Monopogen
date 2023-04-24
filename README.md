# Monopogen: SNV calling from single cell sequencing data

  
  
## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)
* [Introduction](#introduction)
* [Installation](#installation)
* [Data preprocess](#preProcess)
* [Germline SNV calling](#germline)
* [Putative Somatic SNV calling](#somatic)
* [FAQ](#FAQ)
* [Citation](#Citation)

* [Input requirements](#input-requirements)
  * [Sequencing Data](#sequencing-data)
  * [Alignment Files](#alignment-files)
  * [VCF Files](#vcf-files)
* [Outputs](#outputs)
  * [Variant prediction](#variant-prediction)
    * [Germline](#germline)
    * [Somatic](#somatic)
  * [Statistics](#statistics)
* [Run configuration and Execution](#run-configuration-and-execution)
  * [Configuration](#configuration)
    * [Somatic configuration example](#somatic-configuration-example)
    * [Germline configuration example](#germline-configuration-example)
    * [General configuration options](#general-configuration-options)
    * [Advanced configuration options](#advanced-configuration-options)
  * [Execution](#execution)
    * [Advanced execution options](#advanced-execution-options)
  * [Extended use cases](#extended-use-cases)
    * [Improving runtime for references with many short contigs, such as GRCh38](#improving-runtime-for-references-with-many-short-contigs-such-as-grch38)
    * [Exome/Targeted](#exometargeted)
    * [De novo variant calling](#de-novo-variant-calling)
    * [RNA-Seq](#rna-seq)
    * [Heteroplasmic/pooled calling](#heteroplasmicpooled-calling)
    * [Somatic callability](#somatic-callability)
* [Special Topics](#special-topics)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)

## Introduction
**Monopogen** is an analysis package for SNV calling from single-cell sequencing, developed and maintained by [Ken chen's lab](https://sites.google.com/view/kchenlab/Home) in MDACC. `Monopogen` works on sequencing datasets generated from single cell RNA 10x 5', 10x 3', smartseq, single ATAC-seq technoloiges, scDNA-seq etc. 

<image src="./example/Fig1.png" width="600"> 
  
It is composed of three modules: 
* **Data preprocess**. This module removes reads with high alignment mismatches from single cell sequencing and also makes data formats compatiable with Monopongen.
* **Germline SNV calling**. Given the sparsity of single cell sequencing data, we leverage linkage disequilibrium (LD) from external reference panel(such as 1KG3, TopMed) to improve both SNV calling accuracy and detection sensitivity. 
* **Putative somatic SNV calling**. We extended the machinery of LD refinement from human population level to cell population level. We statistically phased the observed alleles with adjacent germline alleles to estimate the degree of LD, taking into consideration widespread sparseness and allelic dropout in single-cell sequencing data, and calculated a probabilistic score as an indicator of somatic SNVs.  The putative somatic SNVs were further genotyped at cell type/cluster level from `Monovar` developed in [Ken chen's lab](https://github.com/KChen-lab/MonoVar).

The output of `Monopogen` will enable 1) ancestry identificaiton on single cell samples; 2) genome-wide association study on the celluar level if sample size is sufficient, and 3) putative somatic SNV investigation.

## Installation
** Dependencies **
* python  (version >= 3.73)
* java (open JDK>=1.8.0)
* pandas>=1.2.3
* pysam>=0.16.0.1
* NumPy>=1.19.5
* sciPy>=1.6.3
* pillow>=8.2.0
** Installation **
Right now Monopogen is avaiable on github, you can install it through github 

`git clone https://github.com/KChen-lab/Monopogen.git`  
`cd Monopogen`  
`pip install -e .`  

## Data preprocess

You can type the following command to get the help information.

`python ./src/Monopogen.py  preProcess --help`

```
usage: Monopogen.py preProcess [-h] -b BAMFILE [-o OUT] -a APP_PATH
                               [-m MAX_MISMATCH] [-t NTHREADS]

optional arguments:
  -h, --help            show this help message and exit
  -b BAMFILE, --bamFile BAMFILE
                        The bam file for the study sample, the bam file should
                        be sorted. If there are multiple samples, each row
                        with each sample (default: None)
  -o OUT, --out OUT     The output director (default: None)
  -a APP_PATH, --app-path APP_PATH
                        The app library paths used in the tool (default: None)
  -m MAX_MISMATCH, --max-mismatch MAX_MISMATCH
                        The maximal alignment mismatch allowed in one reads
                        for variant calling (default: 3)
  -t NTHREADS, --nthreads NTHREADS
                        Number of threads used for SNVs calling (default: 1)
 ```

We provide one example dataset provided the `example/` folder, which includes:
* `A.bam (.bai)`  
  The bam file storing read alignment for sample A.
* `B.bam (.bai)`  
  The bam file storing read alignment for sample B. 
* `CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz`  
  The reference panel with over 3,000 samples in 1000 Genome database. Only SNVs located in chr20: 0-2Mb were extracted in this vcf file. 
* `chr20_2Mb.hg38.fa (.fai)`  
  The genome reference used for read aligments. Only seuqences in chr20:0-20Mb were extracted in this fasta file.

There is a bash script `./test/runPreprocess.sh` to run above example in the folder `test`. You need to prepare the bam file list for option `-b`. If you have multiple sample in this file, you can use more CPUs by setting `-t` to make `Monopogen` faster.  Run the test script as following:
  
```
path="XXy/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

python  ${path}/src/Monopogen.py  preProcess -b bam.lst -o out  -a ${path}/apps -t 8

```
After running the `preProcess` module, there will be bam files after quality controls in the folder `out/Bam/` which will be used for downstream SNV calling.
  
## Germline SNV calling  
 
You can type the following command to get the help information.

`python ./src/Monopogen.py  germline --help`

```
usage: Monopogen.py germline [-h] -r REGION -s
                             {varScan,varImpute,varPhasing,all} [-o OUT] -g
                             REFERENCE -p IMPUTATION_PANEL
                             [-m MAX_SOFTCLIPPED] -a APP_PATH [-t NTHREADS]

optional arguments:
  -h, --help            show this help message and exit
  -r REGION, --region REGION
                        The genome regions for variant calling (default: None)
  -s {varScan,varImpute,varPhasing,all}, --step {varScan,varImpute,varPhasing,all}
                        Run germline variant calling step by step (default:
                        all)
  -o OUT, --out OUT     The output director (default: None)
  -g REFERENCE, --reference REFERENCE
                        The human genome reference used for alignment
                        (default: None)
  -p IMPUTATION_PANEL, --imputation-panel IMPUTATION_PANEL
                        The population-level variant panel for variant
                        imputation refinement, such as 1000 Genome 3 (default:
                        None)
  -a APP_PATH, --app-path APP_PATH
                        The app library paths used in the tool (default: None)
  -t NTHREADS, --nthreads NTHREADS
                        Number of threads used for SNVs calling (default: 1)
 ```

There is a bash script `./test/runGermline.sh` to run above example. You need to prepare the genome region file list for option `-r` with an example shown in `test/region.lst`. Each region is in one row. If you want to call the whole chromosome, you can only specficy the chromosome ID in each row.  Also you can use more CPUs by setting `-t` to make `Monopogen` faster when there are many genome regions.  Run the test script as following:
  
```
python  ${path}/src/Monopogen.py  germline  \
    -a   ${path}/apps -t 8   -r  region.lst \
    -p  ../example/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz  \
    -g  ../example/chr20_2Mb.hg38.fa   -s all  -o out

```
The `germline` module will generate the phased VCF files with name `*.phased.vcf.gz` in the folder `out/germline`. If there are multiple samples in the bam file list from `-b` option in `preProcess` module, the phased VCF files will contain genotypes from multiple samples. The output of phased genotypes are as following:
  
```
##fileformat=VCFv4.2
##filedate=20230422
##source="beagle.27Jul16.86a.jar (version 4.1)"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=AR2,Number=1,Type=Float,Description="Allelic R-Squared: estimated squared correlation betwe
##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage R-Squared: estimated squared correlation betwee
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + P(AA)]">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  19D013_European_F_78    19D014_European_M_84
chr20   68303   .       T       C       .       PASS    .       GT      1|0     1|1
chr20   88108   .       T       C       .       PASS    .       GT      1|1     0|1
chr20   127687  .       A       C       .       PASS    .       GT      1|1     1|1
chr20   153835  .       T       C       .       PASS    .       GT      1|0     1|1
chr20   154002  .       C       T       .       PASS    .       GT      1|1     1|1
chr20   159104  .       T       C       .       PASS    .       GT      1|1     1|1
chr20   167839  .       T       C       .       PASS    .       GT      1|1     1|1
chr20   198814  .       A       T       .       PASS    .       GT      1|0     1|1
chr20   231710  .       T       G       .       PASS    .       GT      1|1     1|1
chr20   237210  .       T       C       .       PASS    .       GT      1|1     1|1
chr20   247326  .       G       A       .       PASS    .       GT      1|1     1|0
chr20   248854  .       T       C       .       PASS    .       GT      1|1     1|0
chr20   255081  .       G       A       .       PASS    .       GT      1|1     1|0
chr20   274893  .       G       C       .       PASS    .       GT      0|1     1|1
chr20   275122  .       G       T       .       PASS    .       GT      0|1     1|1
chr20   275241  .       G       A       .       PASS    .       GT      0|1     1|0
chr20   275361  .       C       T       .       PASS    .       GT      0|1     1|0
chr20   275932  .       A       G       .       PASS    .       GT      0|1     1|0
chr20   276086  .       T       A       .       PASS    .       GT      0|1     1|0

```
  
 
## Run on the HPC or supercomputing environments

If there are multiple single cell RNA samples and you want to use Monopogen on germline SNV calling, you can enable the `-norun` option.

```
python  ${path}/src/Monopogen.py  germline  \
    -a   ${path}/apps -t 8   -r  region.lst \
    -p  ../example/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz  \
    -g  ../example/chr20_2Mb.hg38.fa   -s all  -o out
    --norun

```
The `-norun` module will generate jobs from different regions and you can submit them to HPC based on your own preference. The generated job files will be in `out/Script/`

## Putative somatic SNV calling 



## FAQs 
* ***where to download 1KG3 reference panel (hg38)***
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
  
* ***how to perform downstream PCA-based projection or admixture analysis***  
  PCA-based projection analysis can be peformed using [LASER 2.0](http://csg.sph.umich.edu/chaolong/LASER/)
   
* ***bcftools: error while loading shared libraries: libbz2.so.1.0: not able to open shared object file: No such file or directly***  
  Adding the `apps` folder of `Monopogen` in your library environment 
  
  `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/xx/apps`  
  
* ***AssertionError: Program vcftools cannot be found!***  
  You may set the read/write permission on the folder `xx/apps` as  
  
  `chmod 770 -R  /xx/apps` 
 
## Citation
[Dou J, Tan Y, Wang J, Cheng X, Han KY, Hon CC, Park WY, Shin JW, Chen H, Prabhakar S, Navin N, Chen K. Monopogen : single nucleotide variant calling from single cell sequencing. bioRxiv. 2022 Jan 1](https://www.biorxiv.org/content/10.1101/2022.12.04.519058v1.abstract)




