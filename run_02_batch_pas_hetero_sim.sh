#!/bin/tcsh
#BSUB -x
#BSUB -W 10:00
#BSUB -R span[hosts=1] 
#BSUB -R rusage[mem=42GB]
#BSUB -J hetero
#BSUB -o /share/hmmrs/ngierty/pd/cluster_out/Output_%J_%I.out #error - %J is the job-id %I is the job-array index
#BSUB -e /share/hmmrs/ngierty/pd/cluster_out/Error_%J_%I.err 


conda activate /usr/local/usrapps/hmmrs/ngierty/pd

python code/02_PDM_AOV_SL_multi.py --dname "null_hetero_sim"

conda deactivate
