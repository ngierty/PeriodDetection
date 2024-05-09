#!/bin/tcsh
#BSUB -W 3:00
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=10GB]
#BSUB -J ls_sims
#BSUB -o /share/hmmrs/ngierty/pd/cluster_out/Output_%J_%I.out #error - %J is the job-id %I is the job-array index
#BSUB -e /share/hmmrs/ngierty/pd/cluster_out/Error_%J_%I.err 

setenv XDG_RUNTIME_DIR $TMP

conda activate /usr/local/usrapps/hmmrs/ngierty/pd

python code/03_LSPeriodogramSims.py

conda deactivate
