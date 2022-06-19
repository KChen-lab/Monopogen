python  ../src/scPopGene.py    SCvarCall \
      -b  ../example/chr20_2Mb.rh.filter.sort.bam  \
      -a  ../apps  \
      -c chr20  \
      -o out \
      -i  ../example/cell_cluster.csv   \
      -p  ../example/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz  \
      -r  ../example/chr20_2Mb.hg38.fa  -d 200  -t 0.1  -m 3 -s 5
     
