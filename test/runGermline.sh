path="/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps


python  ${path}/src/Monopogen.py  preProcess -b bam.lst -o out  -a ${path}/apps -t 8


python  ${path}/src/Monopogen.py  germline  \
    -a   ${path}/apps -t 8   -r  region.lst \
    -p  ../example/ \
    -g  ../example/chr20_2Mb.hg38.fa   -m 3 -s all  -o out


     
