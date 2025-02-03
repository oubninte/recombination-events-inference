#!/bin/bash

#SBATCH --account=def-bureau
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-20:30

#SBATCH --output=/home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/TunningRDL/TuningDL_%A.out
#SBATCH --mail-user=samir.oubninte.1@ulaval.ca
#SBATCH --mail-type=FAIL

# to run cmd example :
  #for chrom in 22; do for rep in 1; do sbatch --export=ALL,chrom=${chrom},rep=${rep} **.sh ; done; done

#don't forget to not exceed 1000 jobs

   

echo "**************************** Import module ****************************"
module load StdEnv/2020
module load gcc/9.3.0
module load r-bundle-bioconductor/3.16
module load python/3.10


echo "****************************  Run python script **************************** " 

python -u /home/oubninte/projects/rrg-girardsi/genealogy_sims/results/Samir/P1/RN/DeepLearning_conception_to_infer_Recombination.py
