#!/bin/bash 

#SBATCH --job-name=prlgStats
#SBATCH -p normal
#SBATCH --time=96:00:00

#This is a script to summarize the paralog output from rad_haplotyper by contig with slightly more detail than `*stats.out`

#to run:                 sbatch prlgStats.sbatch <threads> <prefix of input/output>
#to run on vcf files:    sbatch prlgStats.sbatch 40 SpeciesX.A

#enable_lmod
module load parallel
#module load R/gcc7/3.5.3

THREADS=$1
PREFIX=$2

if [[ -z $THREADS ]]; then
	THREADS=8
fi

# Make list of Contigs that are not Complex rad_haplotyper filter, store in array variable
CONTIGS=($(tail -n+3 $PREFIX.Fltr19.stats.out | grep -v 'Complex' | cut -f1))
#CONTIGS=($(tail -n+3 $PREFIX.Fltr19.stats.out | grep -v 'Complex' | cut -f1 | head -n 100))

# Count how many individuals passed paralog test outright for each contig
PassedOutright=($(parallel --no-notice -k -j $THREADS -q grep -c "{}: .*Looks good$" $PREFIX.Fltr19.hap_log.out ::: ${CONTIGS[@]} ))

# Count how many individuals passed paralog test outright with haps for each contig
PassedHAP=($(parallel --no-notice -k -j $THREADS -q grep -c "{}: Looks good$" $PREFIX.Fltr19.hap_log.out ::: ${CONTIGS[@]} ))

# Count how many ind passed paralog test outright with a 1 SNP haps
PassedSNP=($(parallel --no-notice -k -j $THREADS -q grep -c "{}: Single SNP, Looks good$" $PREFIX.Fltr19.hap_log.out ::: ${CONTIGS[@]} ))

# Count how many ind initially failed paralog test, but ultimately passed
PassedHAPinitFail=($(parallel --no-notice -k -j $THREADS -q grep -c "^{}: Problem fixed$" $PREFIX.Fltr19.hap_log.out ::: ${CONTIGS[@]} ))

# Count how many ind initially failed paralog test, but passed after all haplotype reads considered due to minimum count threshold
PassedHAP2Fail=($(parallel --no-notice -k -j $THREADS -q pcregrep -Mc "^{}: Problem fixed\nFixed... huzzah" $PREFIX.Fltr19.hap_log.out ::: ${CONTIGS[@]} ))
PassedHAP2Fail=($(parallel --no-notice -k -j $THREADS "grep -B1 '^Fixed... huzzah' $PREFIX.Fltr19.hap_log.out | grep '^{}: Problem fixed$' " ::: ${CONTIGS[@]} ))

# Calculate how many ind failed paralog test, but passed due to 1st minimum count threshold to consider a paralog
PassedHAP1Fail=($(paste  <(echo ${PassedHAPinitFail[@]} | tr " " "\n") <(echo ${PassedHAP2Fail[@]} | tr " " "\n" ) | awk '{print($1 - $2)}'))

# Calculate how many ind passed paralog test
PASSED=($(paste  <(echo ${PassedOutright[@]} | tr " " "\n") <(echo ${PassedHAPinitFail[@]} | tr " " "\n" ) | awk '{print($1 + $2)}'))

# Count how many ind failed paralog test
FAILED=($(parallel --no-notice -k -j $THREADS "grep -B22 '^Failed again\.\.\.' $PREFIX.Fltr19.hap_log.out | grep -c '^{}$' " ::: ${CONTIGS[@]} ))

# Calculate how many ind took paralog test
NumInd=($(paste  <(echo ${PASSED[@]} | tr " " "\n") <(echo ${FAILED[@]} | tr " " "\n" ) | awk '{print($1 + $2)}'))


# Save to summary file, tab delimited
cat <(echo -e "CONTIG\tNumInd\tPASSED\tPassedOutrightHAP\tPassedOutrightSNP\tPassedHapAfterFail\tPassedHapAfter1Fail\tPassedHapAfter2Fail\tFAILED") \
<(paste <(echo ${CONTIGS[@]} | tr " " "\n") \
	<(echo ${NumInd[@]} | tr " " "\n") \
	<(echo ${PASSED[@]} | tr " " "\n") \
	<(echo ${PassedHAP[@]} | tr " " "\n") \
	<(echo ${PassedSNP[@]} | tr " " "\n") \
	<(echo ${PassedHAPinitFail[@]} | tr " " "\n") \
	<(echo ${PassedHAP1Fail[@]} | tr " " "\n") \
	<(echo ${PassedHAP2Fail[@]} | tr " " "\n") \
	<(echo ${FAILED[@]} | tr " " "\n") \
	) > $PREFIX.paralog_test_results.tsv


###################################################################################################################################################

# Compressed output from `*haps_log.out` for each possible outcome. Used to make patterns above

# # Passed Paralog Test
# dDocent_Contig_81893: Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Unique Observed Haps:
# $VAR1 = [
# ### dDocent_Contig_81893: Looks good

# # Passed Due to 1 SNP in Haplotype
# ### dDocent_Contig_74899: Single SNP, Looks good

# # Failed Paralog Test
# dDocent_Contig_81893: Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Problem- trying to fix...
# dDocent_Contig_81893: haplotype count threshold:2.6
# dDocent_Contig_81893: Corrected Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Unable to rescue
# dDocent_Contig_81893
# Expected haplotypes: 2
# Observed haplotypes: 3
# Failed, trying to recover...
# dDocent_Contig_81893: Observed Haps:
# dDocent_Contig_81893: Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Problem- trying to fix...
# dDocent_Contig_81893: haplotype count threshold:2.6
# dDocent_Contig_81893: Corrected Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Unable to rescue
# dDocent_Contig_81893
# $VAR1 = {
# Expected haplotypes: 2
# Observed haplotypes: 3
# ### Failed again...

# # Failed Initial Paralog Test, but Hap Count Threshold Allows Pass
# dDocent_Contig_81893: Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_81893: Problem- trying to fix...
# dDocent_Contig_81893: haplotype count threshold:4.7
# dDocent_Contig_81893: Corrected Unique Observed Haps:
# $VAR1 = [
# ## dDocent_Contig_81893: Problem fixed


# # Failed Initial Paralog Test, Failed 1st Hap Count Threshold Test, Fails 2nd Paralog Test, but Passes 2nd Hap Count Threshold Because All Haps Considered
# dDocent_Contig_61565: Observed Haps:
# $VAR1 = [
# dDocent_Contig_61565: Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_61565: Problem- trying to fix...
# dDocent_Contig_61565: haplotype count threshold:4.7
# dDocent_Contig_61565: Corrected Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_61565: Unable to rescue
# dDocent_Contig_61565
# $VAR1 = {
# Expected haplotypes: 2
# Observed haplotypes: 3
# Failed, trying to recover...
# dDocent_Contig_61565: Observed Haps:
# $VAR1 = [
# dDocent_Contig_61565: Unique Observed Haps:
# $VAR1 = [
# dDocent_Contig_61565: Problem- trying to fix...
# dDocent_Contig_61565: haplotype count threshold:9.4
# dDocent_Contig_61565: Corrected Unique Observed Haps:
# $VAR1 = [
# ## dDocent_Contig_61565: Problem fixed
# ### Fixed... huzzah!

