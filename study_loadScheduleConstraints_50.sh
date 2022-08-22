#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8gb
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=simmo536@umn.edu
#SBATCH -p small
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

cd ~/MPCloadControl_50
module load matlab
matlab -nodisplay -r "study_loadScheduleConstraints(${SLURM_ARRAY_TASK_ID},7,0.5)"

# sbatch --array=1-100 study_loadScheduleConstraints_50.sh
# dos2unix  study_loadScheduleConstraints_50.sh
