#!/bin/sh

#SBATCH --time=250

#SBATCH --partition=gpu

export KMP_AFFINITY="granularity=fine,compact,1,0"
##export MKL_NUM_THREADS=16

##echo $1/sp_matrix $2 wr $3

srun perl $1/multi_run.pl $1/sp_matrix $2 $3 $4
##echo perl $1/multi_run.pl $1/sp_matrix $2 $3 $4
