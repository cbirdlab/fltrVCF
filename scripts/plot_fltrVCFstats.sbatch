#!/bin/bash

#SBATCH --job-name=pltFltrStats
#SBATCH -p normal
#SBATCH --time=96:00:00

module load R/gcc7/3.5.3

#to run:
#sbatch plot_fltrVCFstats.sbatch filePrefix

PREFIX=$1

#check if file has FilterID and FilterOrder cols
Col2=$(head -1 $PREFIX.fltrStats.dat | cut -f2)
Col3=$(head -1 $PREFIX.fltrStats.dat | cut -f3)

if [ "$Col2" == "FilterID" ] && [ "$Col3" == "FilterOrder" ]; then
echo "$PREFIX.fltrStats.dat has FilterID and FilterOrder columns, creating plots..."
else
echo "$PREFIX.fltrStats.dat is missing FilterID and FilterOrder columns, creating columns..."
paste <(cut -f1 $PREFIX.fltrStats.dat) \
<(cat <(echo -e "FilterID\tFilterOrder") <(cut -f1 $PREFIX.fltrStats.dat | tail -n+2 | tr "." "\t" | cut -f5-6)) \
<(cut -f2- $PREFIX.fltrStats.dat) > $PREFIX.fltrStats.2.dat
mv $PREFIX.fltrStats.2.dat $PREFIX.fltrStats.dat
echo "Columns created, creating plots..."

fi

Rscript plot_fltrVCFstats.R $PREFIX.fltrStats.dat $PREFIX.fltrStats.pdf
