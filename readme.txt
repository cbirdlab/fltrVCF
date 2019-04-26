fltrVCF.bash v4.2 -- a program to filter vcf files with RAD data

Dependencies:
        vcftools
        vcflib
        samtools
        perl
        mawk
        parallel
        rad_haplotyper.pl https://github.com/cbirdlab/rad_haplotyper.git 
        filter_hwe_by_pop_HPC

Reading options from command line:

NAME
        fltrVCF.bash v4.2  -- a program to filter vcf files with RAD data

SYNOPSIS
        fltrVCF.bash [filter settings] [input files] [output file prefix] [parallelization]

DESCRIPTION
        fltrVCF is a tool to filter VCF files created by dDocentHPC. The filters can be run in any order.

        Arguments can be controlled from either the command line or a configuration file.  Filter
        thresholds can only be altered in the config.fltr file. Filters are described and defined
        in the config.fltr file.

        fltrVCF is parallelized where possible, but only runs on one node or computer. MPI is not
        supported.

        fltrVCF requires minor modification to work with dDocent output (nonHPC).  To do so, remove
        ".$CutoffCode" "$CutoffCode." and "$CutoffCode" in order from the script). Forks of both
        filter_hwe_by_pop_HPC.pl and rad_haplotyper.pl can be obtained from cbirdlab on github and
        are tested with fltrVCF and work.

OPTIONS
        [filter settings]
                -f <arg>        if set, controls filters to be run, in order. Argument should be 2
                                 digit numbers separated by spaces. 
				 -f "01 04 02 rmContig"  or  -f 01\ 04\ 02\ 86 
                                 will specify that filters 01, 04, and 02 will be run in succession.
			         Then, filter 86 will remove the contigs that had SNPs filtered by 02.
                                 Filter 86 should only be called after filters that remove SNPs. If 
				 Filter 86 is the only filter called, it will compare the newest vcf to
			         the penultimate vcf, and remove contigs based upon the differences.
				 Filters are described in the config files. If -f is not set, the
                                 config file is used to determine the filters and order. If -f is
                                 set, it will override the config file. The results of each filter will
				 be saved in a separate vcf file.[]
                -s <arg>        file with filter settings [config.fltr.ind]. Should be used.

        [input files]
                -c <arg>        cutoff values used for reference genome [3.3]
                -b <arg>        path to mapping directory with *.bam and *.bed files [../mapping]
                -d <arg>        bed file describing complete data set. Required only if -P is set.
                                 [${b}/mapped.${c}.bed]
                -v <arg>        vcf file to be filtered [${b}/TotalRawSNPs.${c}.vcf]
                -g <arg>        reference genome fasta file [${b}/reference.${c}.fasta]
                -p <arg>        popmap file to use for defining population affiliation
                                 [${b}/popmap.${c}]
                -w <arg>        filter_hwe perl script [filter_hwe_by_pop_HPC.pl]
                -r <arg>        rad_haplotyper perl script [rad_haplotyper.pl v1.19]

        [output file prefix]
                -o <arg>        optional, all output files will be prefixed with this argument []

        [parallelization]
                -P              run every filter in parallel using GNU parallel. Requires *.bed files.
                                 If not set, then only natively-parallel filters will use multiple
                                 threads if -t > 1. Requires -d. [not set]
                -t <arg>        number of threads available for parallel processing [1]

EXAMPLES
        The following command is recommended for most users
                
		fltrVCF.bash -s config.fltr.ind

        The following two commands are the same, the first takes advantage of the defaults,
        the second does not.

                fltrVCF.bash -f "01 02 03" -c 25.10 -o ProjectX.A -t 40

                fltrVCF.bash -f "01 02 03" -c 25.10 -m ../mapping -v ../mapping/TotalRawSNPs.3.6.vcf
                        -p ../mapping/popmap.25.10 -s config.fltr.clean -w filter_hwe_by_pop.pl
                        -r rad_haplotyperHPC116.pl -o ProjectX.A -t 40
SCRIPTS
	Additional scripts for filtering are provided in the scripts subdirectory
	
	fltrVCFstats
		Collects summary statisics from all vcf files in a directory with a matching prefix and
		returns a tidy table where each row is a vcf file and columns are: vcf file name, 
		number of individuals, number of contigs, number of snps, and number of genotype calls 
		based on 0, 1-9, 10-19, etc...  reads.
