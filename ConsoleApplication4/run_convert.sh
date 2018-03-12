#!/bin/sh

#SBATCH --time=1440

#SBATCH --partition=gpu

export KMP_AFFINITY="granularity=fine,compact,1,0"

echo $1/sp_matrix $2 $4 $3

 $1/sp_matrix $2 $4 $3 
