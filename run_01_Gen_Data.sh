#!/bin/tcsh
#BSUB -J gd_pd  #job name 
#BSUB -n 5          #number of threads
#BSUB -W 03:00      #walltime limit: hh:mm
#BSUB -R "rusage[mem=30GB]"
#BSUB -o /share/hmmrs/ngierty/pd/cluster_out/Output_%J_%I.out 
#BSUB -e /share/hmmrs/ngierty/pd/cluster_out/Error_%J_%I.err  #error - %J is the job-id %I is the job-array index 

setenv XDG_RUNTIME_DIR $TMP

conda activate /usr/local/usrapps/hmmrs/ngierty/pd

python code/01_Gen_Data.py


conda deactivate

