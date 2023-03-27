#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p small
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

cd ~/2022Q3-Load-scheduling_yearlyAve
module load matlab
matlab -nodisplay -r \
"startParPool(${SLURM_JOB_CPUS_PER_NODE}); \
study_yearlyAve(${SLURM_ARRAY_TASK_ID})"


# sbatch --array=1-4560 study_yearlyAve.sh
# dos2unix  study_yearlyAve.sh
