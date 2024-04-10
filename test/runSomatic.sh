#BSUB -W 10:00
#BSUB -L /bin/bash
#BSUB -o    /rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen_v1.5/test
#BSUB -cwd  /rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen_v1.5/test
#BSUB -q medium
#BSUB -n 16
#BSUB -M 10
#BSUB -R rusage[mem=40]
#BSUB -J job
#BSUB -P project_name
#PBS -S /bin/bash



path="/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen_v1.5"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${path}/apps

 python  ${path}/src/Monopogen.py  preProcess -b bam.somatic.lst -o somatic  -a ${path}/apps -t 8

python  ${path}/src/Monopogen.py  germline  \
    -a   ${path}/apps -t 8   -r  region.lst \
    -p  ../example/ \
    -g  ../example/chr20_2Mb.hg38.fa   -m 3 -s all  -o somatic



python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  region.lst  -t 22 \
    -i somatic  -l test.csv  -s featureInfo    \
    -g  ../example/chr20_2Mb.hg38.fa   




python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  region.lst  -t 22 \
    -i   somatic  -l test.csv  -s cellScan    \
    -g  ../example/chr20_2Mb.hg38.fa


####### Not work due to fewer marker size 

module load R 
python  ${path}/src/Monopogen.py  somatic  \
    -a   ${path}/apps  -r  region.lst  -t 22 \
    -i   somatic  -l test.csv  -s LDrefinement   \
    -g  ../example/chr20_2Mb.hg38.fa

