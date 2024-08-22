#!/bin/bash
#SBATCH --job-name=  # Put your job name
#SBATCH --time=10-00:00:00

#log files:
#SBATCH -e alphafold_logs/create_individual_features_%A_%a_err.txt  
#SBATCH -o alphafold_logs/create_individual_features_%A_%a_out.txt   

#Adjust this depending on the node
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --mail-user= # change to your email address
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:100

module load alphapulldown/0.30.7

create_individual_features.py \
  --fasta_paths= \
  --save_msa_files=False \
  --output_dir= \
  --use_precomputed_msas=False \
  --max_template_date= \
  --skip_existing=True \
  --seq_index=$SLURM_ARRAY_TASK_ID
