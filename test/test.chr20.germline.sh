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
