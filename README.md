# Monopogen: SNV calling from single cell sequencing data
## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Introduction](#introduction)
* [Installation](#installation)
* [Quick Start](#quick-start)
  * [Data preprocess](#data-preprocess)
  * [Germline SNV calling](#germline-snv-calling)
  * [Putative Somatic SNV calling](#putative-somatic-snv-calling)
* [Germline SNV calling from snRNA-seq](#germline-snv-calling-from-snRNA-seq)
  * [variant calling](#variant-calling)
  * [genotyping accuracy evaluation](#genotyping-accuracy-evaluation)
  * [ancestry identification](#ancestry-identification)
* [Somatic SNV calling from scRNA-seq](#somatic-snv-calling-from-scrna-seq)
* [FAQs](#faqs)
* [Citation](#citation)

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
**Dependencies**
* python  (version >= 3.73)
* java (open JDK>=1.8.0)
* pandas>=1.2.3
* pysam>=0.16.0.1
* NumPy>=1.19.5
* sciPy>=1.6.3
* pillow>=8.2.0
  
**Installation**
  
Right now Monopogen is avaiable on github, you can install it through github 

`git clone https://github.com/KChen-lab/Monopogen.git`  
`cd Monopogen`  
`pip install -e .`  

## Quick Start 

We provide an example dataset provided the `example/` folder, which includes:
* `A.bam (.bai)`  
  The bam file storing read alignment for sample A.
* `B.bam (.bai)`  
  The bam file storing read alignment for sample B. 
* `CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz`  
  The reference panel with over 3,000 samples in 1000 Genome database. Only SNVs located in chr20: 0-2Mb were extracted in this vcf file. 
* `chr20_2Mb.hg38.fa (.fai)`  
  The genome reference used for read aligments. Only seuqences in chr20:0-2Mb were extracted in this fasta file.
There are three test scripts in the `test/` folder `test/runPreprocess.sh`, `test/runGermline.sh`, `test/runSomatic.sh` for quick start of `Monopogen`

### Data preprocess ###

You can type the following command to get the help information.

```
path="XXX/Monopogen"  # where to download Monopogen
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps
python ${path}/src/Monopogen.py  preProcess --help`
```
Output is 
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
 You need to prepare the bam file list for option `-b`. If you have multiple sample in this file, you can use more CPUs by setting `-t` to make `Monopogen` faster.  Run the test script `test/runPreprocess.sh`
  
```
path="XXy/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

python  ${path}/src/Monopogen.py  preProcess -b bam.lst -o out  -a ${path}/apps -t 8

```
After running the `preProcess` module, there will be bam files after quality controls in the folder `out/Bam/` used for downstream SNV calling.
  
### Germline SNV calling ###
 
You can type the following command to get the help information.

```
python ${path}/src/Monopogen.py  germline --help`
```
The output is 
 
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

You need to prepare the genome region file list for option `-r` with an example shown in `test/region.lst`. Each region is in one row. If you want to call the whole chromosome, you can only specficy the chromosome ID in each row.  Also you can use more CPUs by setting `-t` to make `Monopogen` faster when there are many genome regions.  Run the test script test/runGermline.sh as following:
  
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
### Run on the HPC ###
  
If there are multiple single cell RNA samples and you want to use Monopogen on germline SNV calling, you can enable the `-norun` option.

```
python  ${path}/src/Monopogen.py  germline  \
    -a   ${path}/apps -t 8   -r  region.lst \
    -p  ../example/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz  \
    -g  ../example/chr20_2Mb.hg38.fa   -s all  -o out
    --norun TRUE
```
The `-norun` module will generate jobs from different regions and you can submit them to HPC based on your own preference. The generated job files will be in `out/Script/`

### Putative somatic SNV calling ###

## Germline SNV calling from snRNA-seq
We demonstrate the utilization of Monopogen on germline SNV calling, ancestry identification on snRNA samples from human retina atlas. The 4 retina samples shown in Monopogen methodological paper are `19D013`, `19D014`, `19D015`, `19D016`. Thhe fastq files of these samples can be downloaded with from SRA database [SRR23617370](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/sra?term=SRX19501863), [SRR23617337](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/sra?term=SRX19501879), [SRR23617320](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/sra?LinkName=biosample_sra&from_uid=33441051) and [SRR23617310](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/sra?LinkName=biosample_sra&from_uid=33441045). Here we used `19D013` as an example (analysis on other samples is the same). 
### variant calling
For convenience, we skipped the read alignment step and shared the alignmed bam file (Only reads from chr20 were extracted) in [19D013.snRNA.chr20.bam](https://drive.google.com/file/d/18-vdGY9hxbGP-Mm06IBpCCF-NZI9ltAU/view?usp=share_link) and  [19D013.snRNA.chr20.bam.bai](https://drive.google.com/file/d/1HEozx6gmX2Z05R8nElEFOBUhnqayMgAi/view?usp=share_link). Users also need to prepare for the [genome reference]() used for read alignment (for GRCh38) and [imputation panel from 1KG3](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/). We can prepare for the `bam.lst` and `region.lst` as following

```
less bam.lst 
```
The output is 
```
19D013.snRNA.chr20.bam
```
```
less region.lst
```
The output is 
```
chr20
```
Please make sure all required files available
```
ls
19D013.snRNA.chr20.bam
19D013.snRNA.chr20.bam.bai
bam.lst
CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz
GRCh38.chr20.fa
region.lst
```
The data preprocess step can be run as (~3 mins)
```
path="/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen"
${path}/src/Monopogen.py  preProcess -b bam.lst -o retina  -a ${path}/apps  -t 1
```
The output is 
```
[2023-04-25 16:23:05,747] INFO     Monopogen.py Performing data preprocess before variant calling...
[2023-04-25 16:23:05,747] INFO     Monopogen.py Parameters in effect:
[2023-04-25 16:23:05,748] INFO     Monopogen.py --subcommand = [preProcess]
[2023-04-25 16:23:05,748] INFO     Monopogen.py --bamFile = [bam.lst]
[2023-04-25 16:23:05,748] INFO     Monopogen.py --out = [retina]
[2023-04-25 16:23:05,748] INFO     Monopogen.py --app_path = [/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/apps]
[2023-04-25 16:23:05,748] INFO     Monopogen.py --max_mismatch = [3]
[2023-04-25 16:23:05,748] INFO     Monopogen.py --nthreads = [1]
[2023-04-25 16:23:05,765] DEBUG    Monopogen.py PreProcessing sample 19D013
[2023-04-25 16:25:56,543] INFO     Monopogen.py Success! See instructions above.
```
The germline SNV calling can be run as (~80 mins).
 
```
${path}/src/Monopogen.py  germline  -a ${path}/apps  -r region.lst \
 -p CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz \
 -g  GRCh38.chr20.fa  -m 3 -s all  -o retina
```
The output is 
```
[2023-04-25 16:30:39,749] INFO     Monopogen.py Performing germline variant calling...
[2023-04-25 16:30:39,749] INFO     Monopogen.py Parameters in effect:
[2023-04-25 16:30:39,749] INFO     Monopogen.py --subcommand = [germline]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --region = [region.lst]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --step = [all]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --out = [retina]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --reference = [GRCh38.chr20.fa]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --imputation_panel = [CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --max_softClipped = [3]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --app_path = [/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/apps]
[2023-04-25 16:30:39,749] INFO     Monopogen.py --nthreads = [1]
[2023-04-25 16:30:39,750] INFO     Monopogen.py Checking existence of essenstial resource files...
[2023-04-25 16:30:39,754] INFO     Monopogen.py Checking dependencies...
['bash retina/Script/runGermline_chr20.sh']
[fai_load] build FASTA index.
[mpileup] 1 samples in 1 input files
(mpileup) Max depth is above 1M. Potential memory hog!
Lines   total/split/realigned/skipped:  56054517/437864/36916/0
beagle.27Jul16.86a.jar (version 4.1)
Copyright (C) 2014-2015 Brian L. Browning
Enter "java -jar beagle.27Jul16.86a.jar" for a summary of command line arguments.
Start time: 05:10 PM CDT on 25 Apr 2023

Command line: java -Xmx18204m -jar beagle.jar
  gl=retina/germline/chr20.gl.vcf.gz
  ref=CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz
  chrom=chr20
  out=retina/germline/chr20.gp
  impute=false
  modelscale=2
  nthreads=1
  gprobs=true
  niterations=0

No genetic map is specified: using 1 cM = 1 Mb

reference samples:    3202
target samples:          1

Window 1 [ chr20:60291-64332055 ]
reference markers:   31534
target markers:      31531

Starting burn-in iterations

Window=1 Iteration=1
Time for building model:         2 minutes 8 seconds
Time for sampling (singles):     24 seconds
DAG statistics
mean edges/level: 52     max edges/level: 123
mean edges/node:  1.169  mean count/edge: 123

Window=1 Iteration=2
Time for building model:         2 minutes 18 seconds
Time for sampling (singles):     26 seconds
DAG statistics
mean edges/level: 52     max edges/level: 124
mean edges/node:  1.157  mean count/edge: 123

Window=1 Iteration=3
Time for building model:         2 minutes 23 seconds
Time for sampling (singles):     25 seconds
DAG statistics
mean edges/level: 52     max edges/level: 125
mean edges/node:  1.158  mean count/edge: 123

Window=1 Iteration=4
Time for building model:         2 minutes 17 seconds
Time for sampling (singles):     24 seconds
DAG statistics
mean edges/level: 52     max edges/level: 129
mean edges/node:  1.158  mean count/edge: 123

Window=1 Iteration=5
Time for building model:         2 minutes 11 seconds
Time for sampling (singles):     25 seconds
DAG statistics
mean edges/level: 52     max edges/level: 126
mean edges/node:  1.158  mean count/edge: 123

Window=1 Iteration=6
Time for building model:         2 minutes 8 seconds
Time for sampling (singles):     37 seconds
DAG statistics
mean edges/level: 52     max edges/level: 124
mean edges/node:  1.157  mean count/edge: 123

Window=1 Iteration=7
Time for building model:         2 minutes 10 seconds
Time for sampling (singles):     35 seconds
DAG statistics
mean edges/level: 52     max edges/level: 122
mean edges/node:  1.158  mean count/edge: 123

Window=1 Iteration=8
Time for building model:         2 minutes 13 seconds
Time for sampling (singles):     37 seconds
DAG statistics
mean edges/level: 52     max edges/level: 124
mean edges/node:  1.156  mean count/edge: 123

Window=1 Iteration=9
Time for building model:         2 minutes 13 seconds
Time for sampling (singles):     35 seconds
DAG statistics
mean edges/level: 52     max edges/level: 128
mean edges/node:  1.158  mean count/edge: 123

Window=1 Iteration=10
Time for building model:         2 minutes 15 seconds
Time for sampling (singles):     36 seconds
DAG statistics
mean edges/level: 52     max edges/level: 120
mean edges/node:  1.158  mean count/edge: 123

Number of reference markers:     31534
Number of target markers:        31531
Total time for building model: 22 minutes 17 seconds
Total time for sampling:       5 minutes 4 seconds
Total run time:                29 minutes 44 seconds

End time: 05:40 PM CDT on 25 Apr 2023
beagle.27Jul16.86a.jar (version 4.1) finished
beagle.27Jul16.86a.jar (version 4.1)
Copyright (C) 2014-2015 Brian L. Browning
Enter "java -jar beagle.27Jul16.86a.jar" for a summary of command line arguments.
Start time: 05:40 PM CDT on 25 Apr 2023

Command line: java -Xmx18204m -jar beagle.jar
  gt=retina/germline/chr20.germline.vcf
  ref=CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz
  chrom=chr20
  out=retina/germline/chr20.phased
  impute=false
  modelscale=2
  nthreads=48
  gprobs=true
  niterations=0

No genetic map is specified: using 1 cM = 1 Mb

reference samples:    3202
target samples:          1

Window 1 [ chr20:60291-64331516 ]
reference markers:   23755
target markers:      23755

Starting burn-in iterations

Window=1 Iteration=1
Time for building model:         1 minute 29 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 51     max edges/level: 122
mean edges/node:  1.206  mean count/edge: 126

Window=1 Iteration=2
Time for building model:         1 minute 31 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 115
mean edges/node:  1.194  mean count/edge: 128

Window=1 Iteration=3
Time for building model:         1 minute 29 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 124
mean edges/node:  1.194  mean count/edge: 128

Window=1 Iteration=4
Time for building model:         1 minute 30 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 117
mean edges/node:  1.194  mean count/edge: 128

Window=1 Iteration=5
Time for building model:         1 minute 18 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 119
mean edges/node:  1.193  mean count/edge: 128

Window=1 Iteration=6
Time for building model:         1 minute 20 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 116
mean edges/node:  1.194  mean count/edge: 128

Window=1 Iteration=7
Time for building model:         1 minute 19 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 124
mean edges/node:  1.192  mean count/edge: 128

Window=1 Iteration=8
Time for building model:         1 minute 20 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 119
mean edges/node:  1.194  mean count/edge: 128

Window=1 Iteration=9
Time for building model:         1 minute 23 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 121
mean edges/node:  1.192  mean count/edge: 128

Window=1 Iteration=10
Time for building model:         1 minute 21 seconds
Time for sampling (singles):     0 seconds
DAG statistics
mean edges/level: 50     max edges/level: 117
mean edges/node:  1.193  mean count/edge: 128

Number of markers:               23755
Total time for building model: 14 minutes 0 seconds
Total time for sampling:       2 seconds
Total run time:                15 minutes 7 seconds

End time: 05:55 PM CDT on 25 Apr 2023
beagle.27Jul16.86a.jar (version 4.1) finished
[2023-04-25 17:55:37,771] INFO     Monopogen.py Success! See instructions above.

```
The final output of germline SNVs from `Monopogen` are in the folder `retina/germline/chr20.phased.vcf.gz`. These phased genotypes could be used for downstream ancestry identification, association study, and somatic SNV calling.

```
##fileformat=VCFv4.2
##filedate=20230425
##source="beagle.27Jul16.86a.jar (version 4.1)"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=AR2,Number=1,Type=Float,Description="Allelic R-Squared: estimated squared correlation b
##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage R-Squared: estimated squared correlation be
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + P(AA)]">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  19D013_European_F_78
chr20   60291   .       G       T       .       PASS    .       GT      1|0
chr20   68303   .       T       C       .       PASS    .       GT      1|0
chr20   75250   .       C       T       .       PASS    .       GT      1|0
chr20   88108   .       T       C       .       PASS    .       GT      1|1
chr20   101574  .       G       A       .       PASS    .       GT      1|0
chr20   101576  .       G       A       .       PASS    .       GT      1|1
chr20   159104  .       T       C       .       PASS    .       GT      1|1
chr20   175269  .       T       C       .       PASS    .       GT      1|1
chr20   186086  .       G       A       .       PASS    .       GT      1|1
chr20   186183  .       G       A       .       PASS    .       GT      1|1
chr20   198814  .       A       T       .       PASS    .       GT      1|0
chr20   203580  .       G       A       .       PASS    .       GT      1|1
chr20   213223  .       G       C       .       PASS    .       GT      0|1
chr20   213244  .       A       G       .       PASS    .       GT      0|1
chr20   231710  .       T       G       .       PASS    .       GT      1|1

```
### genotyping accuracy evaluation
We can validate the genotyping accuracy and sensitvity (recall) by comparing Monopogen outputs with matched WGS-based genotypes. Users can download the WGS-based genotypes from chr22 only [19D013.wgs.chr20.vcf](https://drive.google.com/file/d/1u55oZgNiwzj5PXeIHCn4NF9dAQlb9uwk/view?usp=share_link). We use [vcftools](https://vcftools.sourceforge.net/) to compare genotypes of Monopogen to the gold standard. 

```
vcftools --gzvcf  ./retina/germline/chr20.phased.vcf.gz    --diff  19D013.wgs.chr20.vcf   --diff-discordance-matrix --out  19D013  --chr chr20
```
 
The output is 
```
VCFtools - 0.1.15
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --gzvcf ./retina/germline/chr20.phased.vcf.gz
        --chr chr20
        --out 19D013
        --diff 19D013.wgs.chr20.vcf
        --diff-discordance-matrix

Using zlib version: 1.2.3
Versions of zlib >= 1.2.4 will be *much* faster when reading zipped VCF files.
After filtering, kept 1 out of 1 Individuals
Outputting Discordance Matrix
        For bi-allelic loci, called in both files, with matching alleles only...
Non-matching ALT. Skipping all such sites.
Non-matching REF. Skipping all such sites.
Found 23290 sites common to both files.
Found 464 sites only in main file.
Found 85853 sites only in second file.
After filtering, kept 23755 out of a possible 23755 Sites
Run Time = 0.00 seconds
```
`Monopogen` can detect `21.3% (23290/(23290+85853))` germline SNVs although the singel cell data is quite sparisty. Remarkably, the false positive rate is lower than `2% (464/(464+23290))`. The genotype concordance could be further examined based on the overlapped 23290 SNVs by looking at the output of `19D013.diff.discordance_matrix`. 

```
less 19D013.diff.discordance_matrix
-       N_0/0_file1     N_0/1_file1     N_1/1_file1     N_./._file1
N_0/0_file2     0       0       0       0
N_0/1_file2     0       13628   723     0
N_1/1_file2     0       60      8869    0
N_./._file2     0       0       0       0
```
The genotyping concordance is calculated as `97% ((60+723)/(60+723+13628+8869))`. The overall genotyping accuracy could be `95% (0.97*(1-0.02))`

### ancestry identification
Here we demonstrate how we can identify ancestry background on snRNA sample `19D013` based on the output of `Monopogen`. Users can use [LASER/TRACE](http://csg.sph.umich.edu/chaolong/LASER/) software to project `19D013` on HGDP reference panel. The HGDP genotyping panel was already included in the [LASER/TRACE](http://csg.sph.umich.edu/chaolong/LASER/) software. Before that, we need to liftover Monopogen output from GRCh38 to GRCh37 to match the HGDP genotyping coordinates. The fasta file of GRCh37 on chr20 could be downloaded as [GRCh37.chr20.fa](https://drive.google.com/file/d/194eSsL4xRLQwwL_3VcWsNxMziSM8SyMB/view?usp=share_link).

```
chain="${path}/resource/hg38ToHg19.over.chain.gz"
GRCh37_chr20="GRCh37.chr20.fa"
picard="${path}/apps/picard.jar"

java -jar ${picard} CreateSequenceDictionary R=${GRCh37_chr20} O="GRCh37.chr20.dict"
java -Xmx10g  -jar ${picard} LiftoverVcf I=./retina/germline/chr20.phased.vcf.gz    O=./retina/germline/chr20.phased.GRCh37.vcf.gz   R=${GRCh37_chr20}  CHAIN=${chain} REJECT="temp.vcf"   WARN_ON_MISSING_CONTIG=true

```
The output will be as following
```
INFO    2023-04-26 01:45:14     CreateSequenceDictionary

********** NOTE: Picard's command line syntax is changing.
**********
********** For more information, please see:
********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
**********
********** The command line looks like this in the new syntax:
**********
**********    CreateSequenceDictionary -R GRCh37.chr20.fa -O GRCh37.chr20.dict
**********

01:45:14.905 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/apps/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Wed Apr 26 01:45:14 CDT 2023] CreateSequenceDictionary OUTPUT=GRCh37.chr20.dict REFERENCE=GRCh37.chr20.fa    TRUNCATE_NAMES_AT_WHITESPACE=true NUM_SEQUENCES=2147483647 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
[Wed Apr 26 01:45:14 CDT 2023] Executing as jdou1@ldragon2 on Linux 3.10.0-1160.15.2.el7.x86_64 amd64; OpenJDK 64-Bit Server VM 1.8.0_312-b07; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: 2.26.10
[Wed Apr 26 01:45:14 CDT 2023] picard.sam.CreateSequenceDictionary done. Elapsed time: 0.00 minutes.
Runtime.totalMemory()=2058354688
To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp
Exception in thread "main" picard.PicardException: /rsrch3/scratch/bcb/jdou1/scAncestry/retina/bam_backup/monopogen_demo1/GRCh37.chr20.dict already exists.  Delete this file and try again, or specify a different output file.
        at picard.sam.CreateSequenceDictionary.doWork(CreateSequenceDictionary.java:220)
        at picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:308)
        at picard.cmdline.PicardCommandLine.instanceMain(PicardCommandLine.java:103)
        at picard.cmdline.PicardCommandLine.main(PicardCommandLine.java:113)
INFO    2023-04-26 01:45:15     LiftoverVcf

********** NOTE: Picard's command line syntax is changing.
**********
********** For more information, please see:
********** https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
**********
********** The command line looks like this in the new syntax:
**********
**********    LiftoverVcf -I ./retina/germline/chr20.phased.vcf.gz -O ./retina/germline/chr20.phased.GRCh37.vcf.gz -R GRCh37.chr20.fa -CHAIN /rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/resource/hg38ToHg19.over.chain.gz -REJECT temp.vcf -WARN_ON_MISSING_CONTIG true
**********


01:45:15.530 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/apps/picard.jar!/com/intel/gkl/native/libgkl_compression.so
[Wed Apr 26 01:45:15 CDT 2023] LiftoverVcf INPUT=./retina/germline/chr20.phased.vcf.gz OUTPUT=./retina/germline/chr20.phased.GRCh37.vcf.gz CHAIN=/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/resource/hg38ToHg19.over.chain.gz REJECT=temp.vcf WARN_ON_MISSING_CONTIG=true REFERENCE_SEQUENCE=GRCh37.chr20.fa    LOG_FAILED_INTERVALS=true WRITE_ORIGINAL_POSITION=false WRITE_ORIGINAL_ALLELES=false LIFTOVER_MIN_MATCH=1.0 ALLOW_MISSING_FIELDS_IN_HEADER=false RECOVER_SWAPPED_REF_ALT=false TAGS_TO_REVERSE=[AF] TAGS_TO_DROP=[MAX_AF] DISABLE_SORT=false VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
[Wed Apr 26 01:45:15 CDT 2023] Executing as jdou1@ldragon2 on Linux 3.10.0-1160.15.2.el7.x86_64 amd64; OpenJDK 64-Bit Server VM 1.8.0_312-b07; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: 2.26.10
INFO    2023-04-26 01:45:15     LiftoverVcf     Loading up the target reference genome.
INFO    2023-04-26 01:45:16     LiftoverVcf     Lifting variants over and sorting (not yet writing the output file.)
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chr3 that is not part of the target reference.
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chr22 that is not part of the target reference.
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chr9_gl000198_random that is not part of the target reference.
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chrUn_gl000241 that is not part of the target reference.
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chr4 that is not part of the target reference.
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chr14 that is not part of the target reference.
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chr8 that is not part of the target reference.
WARNING 2023-04-26 01:45:16     LiftoverVcf     Encountered a contig, chr12 that is not part of the target reference.
INFO    2023-04-26 01:45:16     LiftoverVcf     Processed 23755 variants.
INFO    2023-04-26 01:45:16     LiftoverVcf     184 variants failed to liftover.
INFO    2023-04-26 01:45:16     LiftoverVcf     99 variants lifted over but had mismatching reference alleles after lift over.
INFO    2023-04-26 01:45:16     LiftoverVcf     1.1913% of variants were not successfully lifted over and written to the output.
INFO    2023-04-26 01:45:16     LiftoverVcf     liftover success by source contig:
INFO    2023-04-26 01:45:16     LiftoverVcf     chr20: 23472 / 23755 (98.8087%)
INFO    2023-04-26 01:45:16     LiftoverVcf     lifted variants by target contig:
INFO    2023-04-26 01:45:16     LiftoverVcf     chr20: 23472
WARNING 2023-04-26 01:45:16     LiftoverVcf     99 variants with a swapped REF/ALT were identified, but were not recovered.  See RECOVER_SWAPPED_REF_ALT and associated caveats.
INFO    2023-04-26 01:45:16     LiftoverVcf     Writing out sorted records to final VCF.
[Wed Apr 26 01:45:16 CDT 2023] picard.vcf.LiftoverVcf done. Elapsed time: 0.02 minutes.
Runtime.totalMemory()=2058354688
```
Given SNVs from chr20 only are not enough to identify individual ancestry, we provided the VCF files [19D013.phased.GRCh37.vcf.gz]() by mering all 22 chromosomes. Users can run other chromosome using the same way as we did . Then we can run `TRACE` to project `19D013` on HGDP panel 
```
zless -S ./retina/germline/19D013.phased.GRCh37.vcf.gz  | awk '{gsub(/\chr/, "")}1'  > 19D013.trace.vcf
/rsrch1/bcb/kchen_group/ytan1/LASER-2.04/vcf2geno/vcf2geno --inVcf 19D013.trace.vcf  --out 19D013.trace
/rsrch1/bcb/kchen_group/ytan1/LASER-2.04/trace  -s 19D013.trace.geno  \
  -g /rsrch1/bcb/kchen_group/ytan1/LASER-2.04/HGDP/HGDP_938.geno \
  -c /rsrch1/bcb/kchen_group/ytan1/LASER-2.04/HGDP/HGDP_938.RefPC.coord  \
```
The output is
```
Analysis started at: Wed Apr 26 02:33:45 2023
The following parameters are available.  Ones with "[]" are in effect:

Available Options
           Input/Output : --inVcf [19D013.trace.vcf], --out [19D013.trace]
          People Filter : --peopleIncludeID [], --peopleIncludeFile []
                          --peopleExcludeID [], --peopleExcludeFile []
            Site Filter : --rangeList [], --rangeFile []
      Auxilary Function : --keepDuplication, --updateID []
Skip duplicated variant site:  [ 1      2257207 .       C       T ]
Skip duplicated variant site:  [ 1      42276367        .       A       T ]
Skip duplicated variant site:  [ 1      144952346       .       G       A ]
Skip duplicated variant site:  [ 1      156706559       .       A       G ]
Skip duplicated variant site:  [ 1      220986516       .       C       T ]
Skip duplicated variant site:  [ 1      237281500       .       G       T ]
Skip duplicated variant site:  [ 2      42518616        .       C       G ]
Skip duplicated variant site:  [ 2      191267889       .       G       C ]
Skip duplicated variant site:  [ 3      65996570        .       T       C ]
Skip duplicated variant site:  [ 4      42593642        .       G       T ]
Skip duplicated variant site:  [ 4      89850542        .       T       C ]
Skip duplicated variant site:  [ 4      162791383       .       T       C ]
Skip duplicated variant site:  [ 5      146778988       .       C       G ]
Skip duplicated variant site:  [ 6      77761603        .       G       A ]
Skip duplicated variant site:  [ 6      148867757       .       T       C ]
Skip duplicated variant site:  [ 7      746451  .       A       G ]
Skip duplicated variant site:  [ 7      2605905 .       C       G ]
Skip duplicated variant site:  [ 7      7845898 .       T       C ]
Skip duplicated variant site:  [ 7      16667686        .       A       T ]
Skip duplicated variant site:  [ 7      26868032        .       A       G ]
Skip duplicated variant site:  [ 7      39209642        .       C       G ]
Skip duplicated variant site:  [ 7      97890297        .       C       A ]
Skip duplicated variant site:  [ 7      129431052       .       C       A ]
Skip duplicated variant site:  [ 7      133500313       .       A       T ]
Skip duplicated variant site:  [ 8      1210501 .       T       G ]
Skip duplicated variant site:  [ 8      2703064 .       A       C ]
Skip duplicated variant site:  [ 8      19239747        .       A       G ]
Skip duplicated variant site:  [ 8      32449785        .       C       T ]
Skip duplicated variant site:  [ 8      107142923       .       G       C ]
Skip duplicated variant site:  [ 8      119585345       .       G       C ]
Skip duplicated variant site:  [ 9      95864623        .       C       A ]
Skip duplicated variant site:  [ 10     89084826        .       C       G ]
Skip duplicated variant site:  [ 11     746992  .       T       G ]
Skip duplicated variant site:  [ 11     7441194 .       T       A ]
Skip duplicated variant site:  [ 11     9904151 .       C       T ]
Skip duplicated variant site:  [ 11     12934969        .       A       C ]
Skip duplicated variant site:  [ 11     82934026        .       C       T ]
Skip duplicated variant site:  [ 12     123882977       .       G       C ]
Skip duplicated variant site:  [ 12     133493988       .       C       G ]
Skip duplicated variant site:  [ 13     107948877       .       C       G ]
Skip duplicated variant site:  [ 15     35071084        .       A       T ]
Skip duplicated variant site:  [ 16     3508704 .       C       G ]
Skip duplicated variant site:  [ 17     15420960        .       T       A ]
Skip duplicated variant site:  [ 18     56371446        .       C       G ]
Skip duplicated variant site:  [ 18     77922913        .       A       C ]
Skip duplicated variant site:  [ 19     1489460 .       C       T ]
Skip duplicated variant site:  [ 22     50273174        .       A       C ]
Total 830699 VCF records have converted successfully
Total 1 people and 830652 markers are outputted

=====================================================================
====    TRACE: fasT and Robust Ancestry Coordinate Estimation    ====
====          Version 1.03, Last updated on Dec/30/2016          ====
====          (C) 2013-2016 Chaolong Wang, GNU GPL v3.0          ====
=====================================================================
Started at: Wed Apr 26 02:33:48 2023

1 individuals are detected in the STUDY_FILE.
830652 loci are detected in the STUDY_FILE.
938 individuals are detected in the GENO_FILE.
Warning: Two datasets have different alleles at locus [8:2929436]: [A,G] vs [A,T].
Warning: Two datasets have different alleles at locus [12:5734319]: [A,G] vs [A,C].
Warning: Two datasets have different alleles at locus [13:109351901]: [T,G] vs [T,A].
632958 loci are detected in the GENO_FILE.
938 individuals are detected in the COORD_FILE.
100 PCs are detected in the COORD_FILE.

Parameter values used in execution:
-------------------------------------------------
STUDY_FILE (-s) 19D013.trace.geno
GENO_FILE (-g)  /rsrch1/bcb/kchen_group/ytan1/LASER-2.04/HGDP/HGDP_938.geno
COORD_FILE (-c) /rsrch1/bcb/kchen_group/ytan1/LASER-2.04/HGDP/HGDP_938.RefPC.coord
OUT_PREFIX (-o) trace
DIM (-k)        2
DIM_HIGH (-K)   20
THRESHOLD (-t)  1e-06
MIN_LOCI (-l)   100
FIRST_IND (-x)  1
LAST_IND (-y)   1
KNN_ZSCORE (-knn)       10
RANDOM_SEED (-seed)     0
NUM_THREADS (-nt)       8
-------------------------------------------------

Wed Apr 26 02:33:54 2023
Identify 85497 loci shared by STUDY_FILE and GENO_FILE.
Exclude 3 loci that have different alleles in two datasets.
The analysis will base on the remaining 85494 shared loci.

Wed Apr 26 02:33:54 2023
Reading reference genotype data ...

Wed Apr 26 02:35:05 2023
Calculating reference covariance matrix ...

Wed Apr 26 02:35:06 2023
Reading reference PCA coordinates ...

Wed Apr 26 02:35:06 2023
Analyzing study individuals ...
Procrustean PCA coordinates are output to 'trace.ProPC.coord'.

Finished at: Wed Apr 26 02:35:06 2023
=====================================================================
```
The PCA coordinates of `19D013` is in the file `trace.ProPC.coord`. 

```
less trace.ProPC.coord
```
 
```
popID   indivID L       K       t       Z       PC1     PC2
19D013_European_F_78    19D013_European_F_78    2546    20      0.98445 7.45623 93.869  164.372
```
We can show it on the HGDP PCA plot as 

```
Rscript ${path}/resource/plotTrace.R  ${path}/resource/HGDP.PC.csv  trace.ProPC.coord 19D013_onHGDP
```
The PCA projection plot will be shown generated in `19D013_onHGDP.pdf` 

 
<image src="./example/19D013.onHGDP.png" width="600"> 


 


## Somatic SNV calling from scRNA-seq ##
We demonstrate how the LD refinement model implemented in Monopogen can improve somatic SNV detection from scRNA-seq profiles without matched bulk WGS data available. Please see the [Vignette] []

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




