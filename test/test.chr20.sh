
python  ../src/scPopGene.py    SCvarCall \
      --bamFile  bam.lst \
      --mode multiple \
      --step germline \
      --appPath  ../apps  \
      --chr 20  \
      --out out \
      --cellCluster  ../example/cell_cluster.csv   \
      --imputation-panel  ../example/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz  \
      --reference  ../example/chr20_2Mb.hg38.fa 
