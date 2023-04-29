path="/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps


python  ${path}/src/Monopogen.py  preProcess -b bam.lst -o out  -a ${path}/apps -t 8



     
