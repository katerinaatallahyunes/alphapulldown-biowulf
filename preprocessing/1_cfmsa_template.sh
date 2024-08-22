#!/bin/bash

#SBATCH --job-name=           # put your job name
#SBATCH --time=10-00:00:00

#Adjust this depending on the node
#SBATCH --gres=lscratch:100
#SBATCH --cpus-per-task=32
#SBATCH --mem=128g
#SBATCH --mail-user=  # change to your email address
#SBATCH --mail-type=ALL

module load colabfold alphapulldown

colabfold_search \
    --threads $SLURM_CPUS_PER_TASK \
    fasta.fasta $COLABFOLD_DB pulldown
pushd pulldown
rename_colab_search_a3m.py
popd

create_individual_features.py \
  --fasta_paths= \
  --output_dir= \
  --use_precomputed_msas=True \
  --max_template_date= \
  --use_mmseqs2=True \
  --skip_existing=True \
