fltrVCF.bash v3.2 -- a program to filter vcf files with RAD data

Dependencies:
        vcftools
        vcflib
        samtools
        perl
        mawk
        parallel
        rad_haplotyper116HPC rad
        filter_hwe_by_pop_HPC

Reading options from command line:

NAME
        fltrVCF.bash v3.2  -- a program to filter vcf files with RAD data

SYNOPSIS
        fltrVCF.bash [filter settings] [input files] [output file prefix] [parallelization]

DESCRIPTION
        fltrVCF is a tool to filter VCF files created by dDocentHPC. The filters can be run in any order.

        Arguments can be controlled from either the command line or a configuration file.  Filter
        thresholds can only be altered in the config.fltr file. Filters are described  and defined
                in the config.fltr file.

        fltrVCF is parallelized where possible, but only runs on one node or computer. MPI is not
        supported.

        fltrVCF requires minor modification to work with dDocent output.  To do so, remove
        ".$CutoffCode" "$CutoffCode." and "$CutoffCode" in order from the script). Both
        filter_hwe_by_pop_HPC.pl and rad_haplotyperHPC116.pl are tested with fltrVCF and work. Use of
        other versions is possible, and will be neccessary if filtering data created by dDocent rather
        than dDocentHPC, but is not supported.

OPTIONS
        [filter settings]
                -f <arg>        if set, controls filters to be run, in order. Argument should be 2
                                 digit numbers separated by spaces. -f "01 04 02"  or  -f 01 04 02
                                 will specify that filters 01, 04, and 02 will be run in succession.
                                 Filters are described in the config files. If -f is not set, the
                                 config file is used to determine the filters and order. If -f is
                                 set, it will override the config file. []
                -s <arg>        file with filter settings [config.fltr.clean.ind]

                [input files]
                -c <arg>        cutoff values used for reference genome [3.3]
                -b <arg>        path to mapping directory with *.bam and *.bed files [../mapping] -d <arg>
		                bed file describing complete data set. Required only if -P is set.
                                 [${b}/mapped.${c}.bed]
                -v <arg>        vcf file to be filtered [${b}/TotalRawSNPs.${c}.vcf]
                -g <arg>        reference genome fasta file [${b}/reference.${c}.fasta]
                -p <arg>        popmap file to use for defining population affiliation
                                 [${b}/popmap.${c}]
                -w <arg>        filter_hwe perl script [filter_hwe_by_pop_HPC.pl]
                -r <arg>        rad_haplotyper perl script [rad_haplotyperHPC116.pl]

        [output file prefix]
                -o <arg>        optional, all output files will be prefixed with this argument []

        [parallelization]
                -P              run every filter in parallel using GNU parallel. Requires *.bed files.
                                If not set, then only natively-parallel filters will use multiple
                                threads if -t > 1. Requires -d. [not set]
                -t <arg>        number of threads available for parallel processing [1]

EXAMPLES
        The following command is recommended for most users
                fltrVCF.bash -P -s config.fltr.ind

        The following two commands are the same, the first takes advantage of the defaults,
        the second does not.

                fltrVCF.bash -f "01 02 03" -c 25.10 -o ProjectX.A -t 40

                fltrVCF.bash -f "01 02 03" -c 25.10 -m ../mapping -v ../mapping/TotalRawSNPs.3.6.vcf
                        -p ../mapping/popmap.25.10 -s config.fltr.clean -w filter_hwe_by_pop.pl
                        -r rad_haplotyperHPC116.pl -o ProjectX.A -t 40
