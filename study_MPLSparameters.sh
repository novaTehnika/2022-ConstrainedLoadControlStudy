#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --mem=16gb
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p small
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

cd ~/2022Q3-Load-scheduling
module load matlab
matlab -nodisplay -r \
"parpool('local',$SLURM_JOB_CPUS_PER_NODE); \
study_MPLSparameters(${SLURM_ARRAY_TASK_ID},$SS)"


# sbatch --array=1-100 --export=SS=7 study_MPLSparameters.sh
# dos2unix study_MPLSparameters.sh