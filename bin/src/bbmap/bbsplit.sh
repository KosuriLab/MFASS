#!/bin/bash -l

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
calcXmx () {
	x=$(ulimit -v)
	#echo "x=$x"
	HOSTNAME=`hostname`
	y=1
	if [[ $x == unlimited ]] || [[ $HOSTNAME == gpint* ]]; then
		#echo "ram is unlimited"
		echo "This system does not have ulimit set, so max memory cannot be determined.  Attempting to use 4G." 1>&2
		echo "If this fails, please set ulimit or run this program qsubbed or from a qlogin session on Genepool." 1>&2
		y=4
	else
		mult=75;
		if [ $x -ge 1000000000 ]; then
			mult=85
			#echo "ram is 1000g+"
		elif [ $x -ge 500000000 ]; then
			mult=85
			#echo "ram is 500g+"
		elif [ $x -ge 250000000 ]; then
			mult=85
			#echo "ram is 250g+"
		elif [ $x -ge 144000000 ]; then
			mult=85
			#echo "ram is 144g+"
		elif [ $x -ge 120000000 ]; then
			mult=85
			#echo "ram is 120g+"
		elif [ $x -ge 40000000 ]; then
			mult=80
			#echo "ram is 40g+"
		else
			mult=85
			#echo "ram is under 40g"
		fi
		y=$(( ((x-500000)*mult/100)/1000000 ))
	fi
	#echo "y=$y"
	z="-Xmx${y}g"
	
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
		fi
	done
}
calcXmx "$@"

function bbsplit() {
	module unload oracle-jdk
	module unload samtools
	module load oracle-jdk/1.7_64bit
	module load pigz
	module load samtools
	local CMD="java -ea $z -cp $CP align2.BBSplitter build=1 overwrite=true match=long fastareadlen=500 minhits=2 minratio=0.9 maxindel=20 trim=both untrim=true $@"
	echo $CMD >&2
	$CMD
}

function usage(){
	echo "BBSplit / BBMap v30"
	echo "Written by Brian Bushnell, from Dec. 2010 - present"
	echo "Last modified December 17, 2013"
	echo ""
	echo "Description:  Maps reads to multiple references simultaneously."
	echo "Outputs reads to a file for the reference they best match, with multiple options for dealing with ambiguous mappings."
	echo ""
	echo "To index:  	bbsplit.sh build=<1> ref_x=<reference fasta> ref_y=<another reference fasta>"
	echo "To map:    	bbsplit.sh build=<1> in=<reads> out_x=<output file> out_y=<another output file>"
	echo ""
	echo "To be concise, and do everything in one command:"
	echo "bbsplit.sh ref=x.fa,y.fa in=reads.fq basename=o%.sam"
	echo "which is equivalent to"
	echo "bbsplit.sh build=1 in=reads.fq ref_x=x.fa ref_y=y.fa out_x=ox.sam out_y=oy.sam"
	echo ""
	echo "in=stdin will accept reads from standard in, and out=stdout will write to standard out,"
	echo "but file extensions are still needed to specify the format of the input and output files."
	echo "e.g. in=stdin.fa.gz will read gzipped fasta from standard in; out=stdout.sam.gz will write gzipped sam."
	echo ""	
	echo "Indexing Parameters (required when building the index):"
	echo "ref_<name>=<ref.fasta>  	Specify the reference sequence for the given name; e.g., ref_ecoli=ecoli.fasta"
	echo "                      	These can also be comma-delimited lists of files; e.g., ref_a=a1.fa,a2.fa,a3.fa"
	echo "build=<1>        		If multiple references are indexed in the same directory, each needs a unique build ID."
	echo "k=<13>           		Kmer length, range 8-15.  Longer is faster but uses more memory.  Shorter is more sensitive."
	#echo "                 		I suggest 13 for most cases; 14 for large genomes >3GB; and 12 for PacBio or cross-species mapping."
	#echo "                 		You can have multiple kmer lengths per build number in the same directory."
	#echo "startpad=<2000> 		Pad the beginning of the reference array this many Ns prior to the first scaffold."
	#echo "stoppad=<8000>  		Pad the end of the reference array this many Ns after the end of the last scaffold."
	#echo "midpad=<300>    		Pad this many Ns between adjacent scaffolds.  Higher is better, but wastes memory with tons of tiny scaffolds."
	#echo "colorspace=<false>		Set to true to build a SOLiD colorspace index.  Probably does not work any more."
	#echo "path=<.>         		Specify the location to write the index, if you don't want it in the current directory."
	#echo "minscaf=<1>      		Throw away scaffolds shorter than this when indexing."
	echo ""
	echo "Input Parameters:"
	echo "build=<1>        		Designate index to use.  Corresponds to the number specified when building the index."
	echo "in=<reads.fq>    		Primary reads input; required parameter."
	echo "in2=<reads2.fq>  		For paired reads in two files."
	echo "qin=<auto>       		Set to 33 or 64 to specify input quality value ASCII offset."
	echo "interleaved=<auto>  		True forces paired/interleaved input; false forces single-ended mapping."
	echo "                 		If not specified, interleaved status will be autodetected from read names."
	#echo "fastareadlen=<500>		Break up FASTA reads longer than this.  Max is 500.  Only works for FASTA input."
	#echo "parsecustom=<false>		Specially process read headers from my random read generator, to determine true and false positive rates."
	#echo ""
	#echo "Sampling Parameters:"
	#echo "reads=<-1>       		Set to a positive number N to only process the first N reads, then quit.  -1 means use all reads."
	#echo "idmodulo=<1>     		Set to a number N to only map every Nth read (for deterministic sampling)."
	#echo "samplerate=<1>   		Set to a number from 0-1 to randomly select that fraction of reads for mapping."
	#echo "sampleseed=<1>	  	Set to a positive number N set the RNG seed for sampling at the samplerate,"
	#echo "                 		or a negative number to select a random seed (for nondeterministic sampling)."
	echo ""
	echo "Mapping Parameters:"
	echo "maxindel=<16000> 		Don't look for indels longer than this.  Lower is faster.  Set to >=100k for RNA-seq."
	echo "minratio=<0.56>  		Fraction of max alignment score required to keep a site.  Higher is faster."
	echo "minhits=<1>      		Minimum number of seed hits required for candidate sites.  Higher is faster."
	echo "k=<13>           		Key length for index.  Higher is faster (for large genomes) but uses more RAM.  Max is 15."
	echo "local=<f>      		Set to true to use local, rather than global, alignments.  This will soft-clip ugly ends of poor alignments."
	#echo "perfectmode=<false>		Allow only perfect mappings when set to true (very fast)."
	#echo "semiperfectmode=<false>  	Allow only perfect and semiperfect (perfect except for N's and off-end-of-contig) mappings."
	#echo "threads=<?>      		Set to number of threads desired.  By default, uses all cores available."
	echo "ambiguous=<best> 		Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations)."
	echo "                 			best	(use the first best site)"
	echo "                 			toss	(consider unmapped)"
	echo "                 			random	(select one top-scoring site randomly)"
	echo "                 			all	(retain all top-scoring sites.  Does not work yet with SAM output)"
	echo "ambiguous2=<best> 		Set behavior only for reads that map ambiguously to multiple different references."
	echo "                 		Normal 'ambiguous=' controls behavior on all ambiguous reads;"
	echo "                 		Ambiguous2 excludes reads that map ambiguously within a single reference."
	echo "                 			best	(use the first best site)"
	echo "                 			toss	(consider unmapped)"
	echo "                 			all	(write a copy to the output for each reference to which it maps)"
	echo "                 			split	(write a copy to the AMBIGUOUS_ output for each reference to which it maps)"
	#echo "kfilter=<-1>    		Set to a positive number N to require minimum N contiguous matches for a mapped read."
	#echo "samestrandpairs=<false>	Specify whether paired reads should map to the same strand or opposite strands."
	#echo "requirecorrectstrand=<true>	Forbid pairing of reads without correct strand orientation."
	#echo "rcompmate=<false>		Reverse complement second read in each pair prior to mapping."
	#echo "pairlen=<24000>  		Set max allowed distance between paired reads.  (insert size)=(pairlen)+(read1 length)+(read2 length)"
	echo "trim=<false>		 	Quality-trim ends to Q5 before mapping.  Options are 'l' (left), 'r' (right), and 'lr' (both)."
	echo "untrim=<false> 		Undo trimming after mapping.  Untrimmed bases will be soft-clipped in cigar strings."
	echo ""
	echo "Output Parameters:"
	echo "out_<name>=<file>		Output reads that map to the reference <name> to <file>."
	echo "outputunmapped=<true>		Set to false if unmapped reads should not be printed (saves time and disk space)."
	echo "mdtag=<false>    		Set to true for applications that need an MD tag in SAM files.  Not recommended for RNAseq on euks."
	echo "basename=prefix%suffix	Equivalent to multiple out_%=prefix%suffix expressions, in which each % is replaced by the name of a reference file."
	#echo "ordered=<false>  		Set to true to output reads in same order as input.  Slower and uses more memory."
	#echo "                 		Only relevant with multiple mapping threads."
	#echo "showprogress=<0> 		Set to a positive number N to print a '.' once per N reads processed."
	#echo "match=<short>    		Set to 'none' to skip generation of cigar strings. " 
	#echo "                 		'none' is faster, prevents generation of match and error rate statistics."
	#echo "overwrite=<false>		Allow process to overwrite existing files."
	echo "secondary=<f>    		Print secondary alignments."
	echo "maxsites=<5>     		Maximum number of total alignments to print per read.  Only relevant when secondary=t."
	echo "quickmatch=<f>    		Generate cigar strings more quickly.  Must be true to generate secondary site cigar strings."
	echo "bs=<file>        		Write a shell script to 'file' that will turn the sam output into a sorted, indexed bam file."
	echo "scafstats=<file>     		Write statistics on how many reads mapped to which scaffold to this file."
	echo "refstats=<file>      		Write statistics on how many reads mapped to which reference to this file (for BBSplitter)."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "For more information, please consult /global/projectb/sandbox/gaag/bbtools/docs/readme.txt or the bbmap.sh usage information."
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbsplit "$@"
