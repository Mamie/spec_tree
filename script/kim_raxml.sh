#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=kim
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --error=log/kim.%A.err
#SBATCH --output=log/kim.%A.out
#SBATCH --mail-user=mamie.wang@yale.edu
# reference script from https://rcc.uchicago.edu/docs/running-jobs/array/index.html
module load miniconda/4.7.10
source activate r_env

python kim_raxml.py
