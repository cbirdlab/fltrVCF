#!/bin/bash

#SBATCH --job-name=fltrRadHap
#SBATCH --output=fltrRadHap%j.out
#SBATCH -p cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=llopezdemesa@tamucc.edu
#SBATCH --mail-type=end    # email me when the job finish

module load vcftools/0.1.15
module load ddocent

#You need to make a list using the results from the RadHaplotyper filter (filter 19). Decide cut offs for possible paralogs, max of haplotypes, missing 
#genotype and low coverage. The cut offs will change according to your specific data set. List: CgsToFilter.txt

#Write the cut offs in the echo bellow
echo 'date' "Filtering contigs using resuts from RadHaplotyper: possible paralogs >10, missGen+LowCov >100, haplotypes >100"

#Replace the input file, which should be the input file that was used to run the radhaplotyper filter.
zgrep '^dDocent' opihiSK2014A.10.25.Fltr07.vcf.gz | cut -f 1,2 | grep -w -f CgsToFilter.txt > contigsRadHap.txt

vcftools --gzvcf  opihiSK2014A.10.25.Fltr07.vcf.gz --exclude-positions contigsRadHap.txt --recode --recode-INFO-all --out opihiSK2014A.10.25.Fltr19 

#Adjust input files

bgzip -c opihiSK2014A.10.25.Fltr19.recode.vcf > opihiSK2014A.10.25.Fltr19.recode.vcf.gz
tabix -p vcf opihiSK2014A.10.25.Fltr19.recode.vcf.gz
