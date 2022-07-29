#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p small
#SBATCH -o %j.out
#SBATCH -e %j.err

cd ~/MPCloadControl
module load matlab
matlab -nodisplay -r "study_MPCparameters(${SLURM_ARRAY_TASK_ID})"

# sbatch --array=1-100 study_MPCparameters.sh
# dos2unix  study_MPCparameters.sh
