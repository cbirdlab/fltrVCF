PREFIX=PIRE_SiganusSpinus.L.5.5
REF=../mkVCF2/reference.5.5.fasta
THRESHOLD=1
THRESHOLDb=2

# list contigs with >= X paralogs and Y rescued paralogs, then save them into a variable, note you may have to change the column identifiers ($10, $7) because of changes I've made to the `prlgStats.sbatch` script since creating my `*paralog_test_results.tsv`
PrlgCntgFltr1=$(tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 > X {print $1;}' )
PrlgCntgFltr2=$(tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v Y=$THRESHOLDb '$7> Y {print $1;}' )

# Make a third variable with the list of haps with extreme amounts of missing data
PrlgCntgFltr3=$(tail -n+3 PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | grep 'Missing data$' | cut -f1 )

# Concatenate the contigs in the PrlgCntgFltr vars and remove duplicates to make a list of contigs to remove and convert to regex patterns for vcf file
cat <(echo $PrlgCntgFltr1 | tr " " "\n") <(echo $PrlgCntgFltr2 | tr " " "\n") <(echo $PrlgCntgFltr3 | tr " " "\n") | sort | uniq | sed -e 's/^/\^/g' -e 's/$/\$/g' > PIRE_SiganusSpinus.L.5.5.Fltr19.remove.contigs

# Get a list of the contigs that passed Filter 32 and remove those that failed Filter 19 and convert the contig list to a regex patterns for fasta file
grep '^dDocent_Contig' PIRE_SiganusSpinus.L.5.5.Fltr32.6.vcf | cut -f1 | uniq | grep -vf PIRE_SiganusSpinus.L.5.5.Fltr19.remove.contigs | sed -e 's/^/\^>/g' -e 's/$/\$/g' > PIRE_SiganusSpinus.L.5.5.Fltr19.keep.contigs

# Filter the reference genome for the contigs that passed filters
grep -A1 -f PIRE_SiganusSpinus.L.5.5.Fltr19.keep.contigs ../mkVCF2/reference.5.5.fasta | grep -v '^--$' > PIRE_SiganusSpinus.L.5.5.probes4development.fasta
