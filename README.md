# Monopogen
SNV calling from single cell sequencing data

<image src="./example/Fig1.png" width="400"> 

**Monopogen** is an analysis package for SNV calling from single-cell sequencing, developed and maintained by [Ken chen's lab](https://sites.google.com/view/kchenlab/Home) in MDACC. `Monopogen` works on sequencing datasets generated from single cell RNA 10x 5', 10x 3', smartseq, single ATAC-seq technoloiges, scDNA-seq etc. 
It is composed of three modules: 
* **Data preprocess**. This module removes reads with high alignment mismatches from single cell sequencing and also makes data formats compatiable with Monopongen.
* **Germline SNV calling**. Given the sparsity of single cell sequencing data, we leverage linkage disequilibrium (LD) from external reference panel(such as 1KG3, TopMed) to improve both SNV calling accuracy and detection sensitivity. 
* **Putative somatic SNV calling**. We extended the machinery of LD refinement from human population level to cell population level. We statistically phased the observed alleles with adjacent germline alleles to estimate the degree of LD, taking into consideration widespread sparseness and allelic dropout in single-cell sequencing data, and calculated a probabilistic score as an indicator of somatic SNVs.  The putative somatic SNVs were further genotyped at cell type/cluster level from `Monovar` developed in [Ken chen's lab](https://github.com/KChen-lab/MonoVar).

The output of `Monopogen` will enable 1) ancestry identificaiton on single cell samples; 2) genome-wide association study on the celluar level if sample size is sufficient, and 3) putative somatic SNV investigation.


## 1. Dependencies
* python  (version >= 3.73)
* java (open JDK>=1.8.0)
* pandas>=1.2.3
* pysam>=0.16.0.1
* NumPy>=1.19.5
* sciPy>=1.6.3
* pillow>=8.2.0
## 2. Installation 
Right now Monopogen is avaiable on github, you can install it through github 

`git clone https://github.com/KChen-lab/Monopogen.git`  
`cd Monopogen`  
`pip install -e .`  

## 3. Usage of Monopogen
  
## 3.1 Data preprocess

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

There is a bash script `./test/runPreprocess.sh` to run above example in the folder `test`. You need to prepare the bam file list for option `-b`. If you have multiple sample in the list file, `Monopogen` will run the joint calling which can increase the SNV calling accuracy and sensitivity. Run the test script as following:
  
```
path="XXy/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

python  ${path}/src/Monopogen.py  preProcess -b bam.lst -o out  -a ${path}/apps -t 8

```

  

```
path="XX/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

python  ../src/Monopogen.py    germline  \
      -b  ../example/chr20_2Mb.rh.filter.sort.bam  \
      -y  single  \
      -t  all  \
      -a  ../apps  \
      -c  chr20  \
      -o  out \
      -d  10 \
      -p  ../example/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz  \
      -r  ../example/chr20_2Mb.hg38.fa   -m 3 -s 5


```

## Output

The most important results are in the folder `out/germline`: 

* `chr20.germline.vcf`  
In this vcf file, the germline SNVs were genotyped for the study sample 
```
##fileformat=VCFv4.2
##filedate=20220808
##source="beagle.27Jul16.86a.jar (version 4.1)"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=AR2,Number=1,Type=Float,Description="Allelic R-Squared: estimated squared correlation between most probable REF dose and true REF dose">
##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + P(AA)]">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  test1   test2
chr20   60291   .       G       T       .       PASS    AR2=0.00;DR2=0.00;AF=0.50       GT:DS:GP        0/1:1:0,1,0     0/1:1:0,1,0
chr20   64506   .       C       T       .       PASS    AR2=0.00;DR2=0.00;AF=0.0051     GT:DS:GP        0/0:0.01:0.99,0.01,0    0/0:0.01:0.99,0.01,0
chr20   68303   .       T       C       .       PASS    AR2=0.00;DR2=0.00;AF=0.50       GT:DS:GP        0/1:1:0,1,0     0/1:1:0,1,0
chr20   75250   .       C       T       .       PASS    AR2=0.00;DR2=0.00;AF=0.50       GT:DS:GP        0/1:1:0,1,0     0/1:1:0,1,0
chr20   88108   .       T       C       .       PASS    AR2=0.00;DR2=0.00;AF=0.87       GT:DS:GP        1/1:1.73:0,0.26,0.73    1/1:1.73:0,0.26,0.73
chr20   159104  .       T       C       .       PASS    AR2=0.00;DR2=0.00;AF=0.9960     GT:DS:GP        1/1:1.99:0,0.01,0.99    1/1:1.99:0,0.01,0.99
chr20   175269  .       T       C       .       PASS    AR2=0.00;DR2=0.00;AF=0.9903     GT:DS:GP        1/1:1.98:0,0.02,0.98    1/1:1.98:0,0.02,0.98
chr20   186086  .       G       A       .       PASS    AR2=0.00;DR2=0.00;AF=0.983      GT:DS:GP        1/1:1.97:0,0.03,0.97    1/1:1.97:0,0.03,0.97
chr20   186183  .       G       A       .       PASS    AR2=0.00;DR2=0.00;AF=0.983      GT:DS:GP        1/1:1.97:0,0.03,0.97    1/1:1.97:0,0.03,0.97
chr20   198814  .       A       T       .       PASS    AR2=0.00;DR2=0.00;AF=0.50       GT:DS:GP        0/1:1:0.05,0.89,0.06    0/1:1:0.05,0.89,0.06
chr20   231710  .       T       G       .       PASS    AR2=0.00;DR2=0.00;AF=0.70       GT:DS:GP        0/1:1.4:0,0.6,0.4       0/1:1.4:0,0.6,0.4
chr20   240166  .       T       G       .       PASS    AR2=0.00;DR2=0.00;AF=0.49       GT:DS:GP        0/1:0.97:0.03,0.97,0    0/1:0.97:0.03,0.97,0
chr20   240377  .       T       G       .       PASS    AR2=0.00;DR2=0.00;AF=0.49       GT:DS:GP        0/1:0.97:0.03,0.97,0    0/1:0.97:0.03,0.97,0
chr20   247326  .       G       A       .       PASS    AR2=0.00;DR2=0.00;AF=0.9964     GT:DS:GP        1/1:1.99:0,0.01,0.99    1/1:1.99:0,0.01,0.99
chr20   248854  .       T       C       .       PASS    AR2=0.00;DR2=0.00;AF=0.73       GT:DS:GP        0/1:1.47:0,0.53,0.47    0/1:1.47:0,0.53,0.47
chr20   254891  .       C       T       .       PASS    AR2=0.00;DR2=0.00;AF=0.51       GT:DS:GP        0/1:1.02:0,0.97,0.03    0/1:1.02:0,0.97,0.03
chr20   255081  .       G       A       .       PASS    AR2=0.00;DR2=0.00;AF=0.85       GT:DS:GP        1/1:1.7:0,0.3,0.7       1/1:1.7:0,0.3,0.7
chr20   260997  .       G       A       .       PASS    AR2=0.00;DR2=0.00;AF=0.40       GT:DS:GP        0/1:0.8:0.2,0.8,0       0/1:0.8:0.2,0.8,0
chr20   265294  .       G       T       .       PASS    AR2=0.00;DR2=0.00;AF=0.918      GT:DS:GP        1/1:1.84:0,0.16,0.84    1/1:1.84:0,0.16,0.84

```
## Run on multiple chromosomes and multiple samples 

Users can submit jobs with multiple chromosomes in the parallele fashion as following:

```
path="XX/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

for chr in {1..22}
do 
  python  ../src/Monopogen.py    germline  \
        -b  ../example/chr${chr}_2Mb.rh.filter.sort.bam  \
        -y  single  \
        -t  all  \
        -a  ../apps  \
        -c  chr${chr}  \
        -o  out \
        -d  10 \
        -p  ../example/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz  \
        -r  ../example/chr${chr}_2Mb.hg38.fa   -m 3 -s 5
done

```



## 7. FAQs 
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
 
  
## 8. Citation
[Dou J, Tan Y, Wang J, Cheng X, Han KY, Hon CC, Park WY, Shin JW, Chen H, Prabhakar S, Navin N, Chen K. Monopogen : single nucleotide variant calling from single cell sequencing. bioRxiv. 2022 Jan 1](https://www.biorxiv.org/content/10.1101/2022.12.04.519058v1.abstract)




