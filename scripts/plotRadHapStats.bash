#!/bin/bash

#this script will plot output from radhaplotyper run from within fltrVCF.bash and is included in newer versions of fltrVCF.bash

#to run:
#bash plotRadHapStats.bash <FileNamePrefix>

#example when filter19 (radhaplotyper) produced the following vcf: Hlobatus.D.5.5.Fltr19.Haplotyped.vcf
#bash plotRadHapStats.bash Hlobatus.D.5.5

$VCF_OUT=$1
$FILTER_ID=19

#histogram of # SNPs
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f2 | sort -rg > sites.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of number of sites per contig."
set ylabel "Number of Contigs"
set xlabel "Number of Sites"
set yrange [0:*]
xmax="`head -1 sites.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'sites.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of # Sites per Contig."
set ylabel "# Sites"
set xlabel "Contigs"
plot 'sites.dat' pt "*" 
pause -1
EOF

#histogram of # haplotypes
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f3 | sort -rg > haplotypes.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of number of haplotypes per contig."
set ylabel "Number of Contigs"
set xlabel "Number of Haplotypes"
set yrange [0:*]
xmax="`head -1 haplotypes.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'haplotypes.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of # Haplotypes per Contig."
set ylabel "# Haplotypes"
set xlabel "Contigs"
plot 'haplotypes.dat' pt "*" 
pause -1
EOF

#scatter plot of sites vs haplotypes
paste <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f3) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f2) > haplo_sites.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Sites vs Haplotypes per Contig"
set ylabel "# Sites"
set xlabel "# Haplotypes"
plot 'haplo_sites.dat' pt "*" 
pause -1
EOF

#histogram of proportion of individuals haplotyped
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f6 | sort -rg > propindhaplotyped.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of the proportion of individuals haplotyped per contig"
set ylabel "Number of Contigs"
set xlabel "Proportion Haplotyped"
set yrange [0:*]
set xrange [0:1]
binwidth= 0.05
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'propindhaplotyped.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of the proportion of individuals haplotyped per contig"
set ylabel "Proportion Haplotyped"
set xlabel "Contigs"
plot 'propindhaplotyped.dat' pt "*" 
pause -1
EOF

#histogram of number of paralogous individuals
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f8 | sort -rg > paralogs.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of Number of Paralogous Individuals per Contig"
set ylabel "Number of Contigs"
set xlabel "Number of Paralogous Individuals per Contig (0-20)"
set yrange [0:*]
set xrange [0:20]
binwidth= 1
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'paralogs.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of # Paralogous Individuals per Contig"
set ylabel "# Paralogs"
set xlabel "Contigs"
plot 'paralogs.dat' pt "*"
pause -1
EOF

#histogram of number of individuals with low cov/geno errs
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f9 | sort -rg > lowcovgenoerr.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of Number of Individuals With Low Cov/ Geno Errs per Contig"
set ylabel "Number of Contigs"
set xlabel "Number of Individuals per Contig"
set yrange [0:*]
xmax="`head -1 lowcovgenoerr.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'lowcovgenoerr.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Number of Individuals With Low Cov/ Geno Errs per Contig"
set ylabel "Number of Individuals per Contig"
set xlabel "Contigs"
plot 'lowcovgenoerr.dat' pt "*"
pause -1
EOF

#histogram of number of individuals with missing genotypes
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f10 | sort -rg > missinggeno.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of Number of Individuals With Missing Genotypes per Contig"
set ylabel "Number of Contigs"
set xlabel "Number of Individuals per Contig"
set yrange [0:*]
xmax="`head -1 missinggeno.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'missinggeno.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Number of Individuals With Missing Genotypes per Contig"
set ylabel "Number of Individuals per Contig"
set xlabel "Contigs"
plot 'missinggeno.dat' pt "*"
pause -1
EOF

#scatter plots of paralogs vs X
paste haplo_sites.dat <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f6) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f8) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f9) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f10) > haplo_sites_prophap_paralogs.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Paralogous Individuals vs Haplotypes per Contig"
set ylabel "# Paralogs"
set xlabel "# Haplotypes"
plot 'haplo_sites_prophap_paralogs.dat' using 1:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs Sites per Contig"
set ylabel "# Paralogs"
set xlabel "# Sites"
plot 'haplo_sites_prophap_paralogs.dat' using 2:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs Proportion of Individuals Haplotyped per Contig"
set ylabel "# Paralogs"
set xlabel "Proportion Haplotyped"
plot 'haplo_sites_prophap_paralogs.dat' using 3:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs # Individuals with Low Cov/ Geno Err per Contig"
set ylabel "# Paralogs"
set xlabel "# Individuals"
plot 'haplo_sites_prophap_paralogs.dat' using 5:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs # Individuals with Missing Genotypes per Contig"
set ylabel "# Paralogs"
set xlabel "# Individuals"
plot 'haplo_sites_prophap_paralogs.dat' using 6:4 with points pt "*"

set title "Scatter plot of # Individuals with Low Cov/ Geno Err vs # Individuals with Missing Genotypes per Contig"
set ylabel "# Individuals w Low Cov / Geno Err"
set xlabel "# Individuals w Missing Geno"
plot 'haplo_sites_prophap_paralogs.dat' using 6:5 with points pt "*"

set title "Scatter plot of # Individuals with Low Cov/ Geno Err vs Proportion of Individuals Haplotyped per Contig"
set ylabel "# Individuals with Low Cov/ Geno Err"
set xlabel "Proportion Individuals Haplotyped"
plot 'haplo_sites_prophap_paralogs.dat' using 3:5 with points pt "*"

set title "Scatter plot of # Individuals with Missing Genotypes vs Proportion of Individuals Haplotyped per Contig"
set ylabel "# Individuals with Missing Genotypes"
set xlabel "Proportion Individuals Haplotyped"
plot 'haplo_sites_prophap_paralogs.dat' using 3:6 with points pt "*"

pause -1
EOF