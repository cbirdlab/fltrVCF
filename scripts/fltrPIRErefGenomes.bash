PREFIX=PIRE_SiganusSpinus.L.5.5
REF=../mkVCF2/reference.5.5.fasta
THRESHOLD=1
THRESHOLDb=2

# list contigs with >= X paralogs and Y rescued paralogs, then save them into a variable, note you may have to change the column identifiers ($10, $7) because of changes I've made to the `prlgStats.sbatch` script since creating my `*paralog_test_results.tsv`
PrlgCntgFltr1=$(tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 > X {print $1;}' )
PrlgCntgFltr2=$(tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v Y=$THRESHOLDb '$7 > Y {print $1;}' )

# Make a third variable with the list of haps with extreme amounts of missing data
PrlgCntgFltr3=$(tail -n+3 PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | grep 'Missing data$' | cut -f1 )

# Concatenate the contigs in the PrlgCntgFltr vars and remove duplicates to make a list of contigs to remove and convert to regex patterns for vcf file
cat <(echo $PrlgCntgFltr1 | tr " " "\n") <(echo $PrlgCntgFltr2 | tr " " "\n") <(echo $PrlgCntgFltr3 | tr " " "\n") | sort | uniq | sed -e 's/^/\^/g' -e 's/$/\$/g' > PIRE_SiganusSpinus.L.5.5.Fltr19.remove.contigs

# Get a list of the contigs that passed Filter 32 and remove those that failed Filter 19 and convert the contig list to a regex patterns for fasta file
grep '^dDocent_Contig' PIRE_SiganusSpinus.L.5.5.Fltr32.6.vcf | cut -f1 | uniq | grep -vf PIRE_SiganusSpinus.L.5.5.Fltr19.remove.contigs | sed -e 's/^/\^>/g' -e 's/$/\$/g' > PIRE_SiganusSpinus.L.5.5.Fltr19.keep.contigs

# Filter the reference genome for the contigs that passed filters
grep -A1 -f PIRE_SiganusSpinus.L.5.5.Fltr19.keep.contigs ../mkVCF2/reference.5.5.fasta | grep -v '^--$' > PIRE_SiganusSpinus.L.5.5.probes4development.fasta

# Filter reference for microsatellites

# Filter reference for NNNNNNNNNNNNNNNNNNNN

# Filter reference by contig length
THRESHOLD=0.01
THRESHOLDb=0.99

paste <(grep '^>'  PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.fasta) \
	<(grep -v '^>' PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.fasta | awk '{ print length }') | \
	sort -nk2 > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv
THRESHOLD=$(cut -f2 PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | awk -v PCT=$THRESHOLD '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
THRESHOLDb=$(cut -f2 PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | awk -v PCT=$THRESHOLDb '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
awk -v LENb=$THRESHOLDb -v LEN=$THRESHOLD '$2 > LENb || $2 < LEN {print $1;}' PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | sed -e 's/^/\^/' -e 's/$/\t/' > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.remove.contigs
cat PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.fasta | paste - - > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.tsv
grep -vf PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.remove.contigs PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.tsv | tr "\t" "\n" > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.${THRESHOLD}-${THRESHOLDb}bp.fasta

# Remove bp prior to restriction site motif in contigs from ref genome
BAR=8

grep  -B1 "^.\{0,$BAR\}TGCAGG" PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.371-618bp.fasta | grep -v '^--' | sed "s/^.\{0,$BAR\}\(TGCAGG\)/N\1/" > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.371-618bp.0-${BAR}TGCAGG.fasta

