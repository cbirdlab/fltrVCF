## fltrVCF v4.4 -- a program to filter VCF files

*Dependencies*

[vcftools](https://vcftools.github.io/index.html)
	* [vcflib](https://github.com/vcflib/vcflib)
	* [samtools](http://www.htslib.org/)
	* perl
	* [mawk](https://invisible-island.net/mawk/)
	* [GNU parallel](https://www.gnu.org/software/parallel/)
	* R (tidyverse, gridExtra)
	* [rad_haplotyper 1.1.5 fork](https://github.com/cbirdlab/rad_haplotyper)

### NAME
```
fltrVCF.bash v4.3  -- a program to filter vcf files with RAD data
```

### SYNOPSIS
```
fltrVCF.bash [filter settings] [input files] [output file prefix] [parallelization]
```

### DESCRIPTION
```
fltrVCF is a tool to filter VCF files created by dDocentHPC. The filters can be run in any order.

Arguments can be controlled from either the command line or a configuration file.  Filter thresholds 
can only be altered in the config.fltr file. Filters are described and defined in the config.fltr file.

fltrVCF is parallelized where possible, but only runs on one node or computer. MPI is not supported.

fltrVCF requires minor modification to work with non dDocentHPC output.  To do so, remove ".$CutoffCode" 
"$CutoffCode." and "$CutoffCode" in order from the script). 
	
Forks of both filter_hwe_by_pop_HPC.pl and rad_haplotyper.pl can be obtained from cbirdlab on github, 
are tested with fltrVCF and work. 
	
Newer filters (identified by the "custom bash" label in the config file) employ R to output plots as 
*.pdf.  For this to work properly, the R scripts in the fltrVCF/scripts dir need to be copied to the 
working directory.
```

### OPTIONS
        [filter settings]
                -f <arg>        if set, controls filters to be run, in order. Argument should be 2-3
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
                -P              (depricated a new method of parallelizaiton is applied by default) 
				 run every filter in parallel using GNU parallel. Requires *.bed files.
                                 If not set, then only natively-parallel filters will use multiple
                                 threads if -t > 1. Requires -d. [not set]
                -t <arg>        number of threads available for parallel processing [1]


### DOWNLOADING & PREPARING TO RUN

I will assume that you can find and install the dependencies listed at the top of this doc. What 
follows is a description of how to get fltrVCF up and running, assuming that the other dependencies are 
in place.

Assumed directory structure:
```
home			either your home dir, or where you keep your ProjectDir
+--fltrVCF		cloned repo
+--rad_haplotyper	cloned repo
+--ProjectDir		directory for your RAD project
|  +--mkBAM		directory where BAM files were made
|  +--mkVCF		directory where VCF file was made (could be mkBAM)
|  +--filterVCF		directory where VCF file is filtered
```

Goto home (or where you keep your ProjectDir.  (The `$` indicates the `bash` cmd prompt, don't type it.)

```bash
$ cd ~
```

Clone this repo and my fork of rad_haplotyper to your computer
		
```bash
$ git clone https://github.com/cbirdlab/fltrVCF.git
$ git clone https://github.com/cbirdlab/rad_haplotyper.git
```

Copy the fltrVCF/fltrVCF.bash and fltrVCF/fltrVCF.sbatch files to your working directory

```bash
$ cp fltrVCF/fltrVCF* ProjectDir/filterVCF
```

For filters with R plotting, copy fltrVCF/scripts/*R to your working directory, 
	
```bash
$ cp fltrVCF/scripts/*R ProjectDir/filterVCF
```

Then run R and install required packages  (The `>` indicates the `R` cmd prompt, don't type it.)
	
```bash
$ R
```
		
```R
> install.packages("tidyverse")
> install.packages("gridExtra")
```

and exit R (hit ctrl+d)   

```R
>
Save workspace image? [y/n/c]: y
```
	
For the default config to work with minimal editing, copy the following to your working dir
	
```bash
$ cp rad_haplotyper/rad_haplotyper.pl ProjectDir/filterVCF
$ cp fltrVCF/filter_hwe_by_pop_HPC.pl ProjectDir/filterVCF
```
		
		
Copy the config file to your working dir and modify it match your directory structure
	
```bash
$ cp fltrVCF/config* ProjectDir/filterVCF
$ nano ProjectDir/filterVCF
```
	
Move to your filterVCF dir

```bash
$ cd ProjectDir/filterVCF
```

### EXAMPLES
The following command is recommended for most users and requires the config file to have the 
correct paths.
                
```bash
$ fltrVCF.bash -s config.fltr.ind
```

The following two commands are the same, the first takes advantage of the defaults,
the second does not.

```bash
# Note that the escape character `\` continues the present line on the next line
fltrVCF.bash \
  -f "01 02 03" \
  -c 25.10 -o ProjectX.A \
  -t 40
```

```bash
$ fltrVCF.bash \
  -f "01 02 03" \
  -c 25.10 \
  -m ../mapping \
  -v ../mapping/TotalRawSNPs.3.6.vcf \
  -p ../mapping/popmap.25.10 -s config.fltr.clean -w filter_hwe_by_pop.pl \
  -r rad_haplotyperHPC116.pl -o ProjectX.A -t 40
```
		
### SCRIPTS

Additional scripts for filtering and viewing output are provided in the scripts subdirectory

* `fltrVCFstats2.sbatch`
	Collects summary statisics from all vcf files in a directory with a matching prefix and
	returns a tidy table where each row is a vcf file and columns are: vcf file name, 
	number of individuals, number of contigs, number of snps, and number of genotype calls 
	based on 0, 1-9, 10-19, etc...  reads.
	
---

### WHAT IS NEW

#### 4.4

* handle vcf files not made by dDocent

---

## CITATION

If you use fltrVCF in a publication, please cite the following: 

Biesack, E. E., Dang, B. T., Ackiss, A. S., Bird, C. E., Chheng, P., Phounvisouk, L., Truong, O.T. & Carpenter, K. E. (2020). Evidence for population genetic structure in two exploited Mekong River fishes across a natural riverine barrier. Journal of Fish Biology, 1-12.

Bird, C.E. (InsertYearOfMostRecentUpdateWhenCloned) fltrVCF. https://github.com/cbirdlab/fltrVCF

