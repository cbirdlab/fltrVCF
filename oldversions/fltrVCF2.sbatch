#!/bin/bash

#SBATCH --job-name=fltrVCF
#SBATCH -p cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=tamucc.edu
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent
module load perl/5.22.0
module load vcftools/0.1.15
module load rad_haplotyper/1.1.5
module load parallel
module load vcflib/1.0
module load R

bash fltrVCF2.bash -s config.fltr.clean.ind