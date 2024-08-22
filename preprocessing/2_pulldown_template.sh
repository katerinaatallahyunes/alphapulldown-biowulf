#!/bin/bash


#SBATCH --job-name=               # put your job name
#SBATCH --time=10-00:00:00

#Adjust this depending on the node
#SBATCH --partition=gpu
#SBATCH --gres=lscratch:100,gpu:a100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128g
#SBATCH --mail-user=atallahyuneska@nih.gov    # change to your email address
#SBATCH --mail-type=ALL
#SBATCH --array=1-$count%10
#SBATCH -o folder_name/folder_name_%j_%A_%a.out
#SBATCH -e folder_name/folder_name_%j_%A_%a.err

module load alphapulldown

run_multimer_jobs.py \
  --mode=pulldown \
  --num_cycle=3 \
  --num_predictions_per_model=1 \
  --output_path= \
  --protein_lists= \
  --monomer_objects_dir= \
  --job_index=$SLURM_ARRAY_TASK_ID
