#!/bin/bash -l
#bbmap in=<infile> out=<outfile>

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
			mult=84
			#echo "ram is 1000g+"
		elif [ $x -ge 500000000 ]; then
			mult=84
			#echo "ram is 500g+"
		elif [ $x -ge 250000000 ]; then
			mult=84
			#echo "ram is 250g+"
		elif [ $x -ge 144000000 ]; then
			mult=84
			#echo "ram is 144g+"
		elif [ $x -ge 120000000 ]; then
			mult=84
			#echo "ram is 120g+"
		elif [ $x -ge 40000000 ]; then
			mult=80
			#echo "ram is 40g+"
		else
			mult=84
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


bbmap() {
	module unload oracle-jdk
	module unload samtools
	module load oracle-jdk/1.7_64bit
	module load pigz
	module load samtools
	local CMD="java -ea $z -cp $CP align2.BBMap build=1 overwrite=true match=long  fastareadlen=500 $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "BBMap v31"
	echo "Written by Brian Bushnell, from Dec. 2010 - present"
	echo "Last modified January 6, 2014"
	echo ""
	echo "Description:  Fast and accurate short-read aligner for DNA and RNA."
	echo ""
	echo "To index:  	bbmap.sh ref=<reference fasta>"
	echo "To map:    	bbmap.sh in=<reads> out=<output sam>"
	echo "To map without an index:  	bbmap.sh ref=<reference fasta> in=<reads> out=<output sam> nodisk"
	echo ""
	echo "in=stdin will accept reads from standard in, and out=stdout will write to standard out,"
	echo "but file extensions are still needed to specify the format of the input and output files."
	echo "e.g. in=stdin.fa.gz will read gzipped fasta from standard in; out=stdout.sam.gz will write gzipped sam."
	echo ""	
	echo "Indexing Parameters (required when building the index):"
	echo "nodisk=f         		Set to true to build index in memory and write nothing to disk except output."
	echo "ref=<file>      		Specify the reference sequence.  Only do this ONCE, when building the index (unless using 'nodisk')."
	echo "build=1          		If multiple references are indexed in the same directory, each needs a unique numeric ID (unless using 'nodisk')."
	echo "k=13             		Kmer length, range 8-15.  Longer is faster but uses more memory.  Shorter is more sensitive."
	#echo "                 		I suggest 13 for most cases; 14 for large genomes >3GB; and 12 for PacBio or cross-species mapping."
	#echo "                 		You can have multiple kmer lengths per build number in the same directory."
	#echo "startpad=2000   		Pad the beginning of the reference array this many Ns prior to the first scaffold."
	#echo "stoppad=8000   		Pad the end of the reference array this many Ns after the end of the last scaffold."
	#echo "midpad=300      		Pad this many Ns between adjacent scaffolds.  Higher is better, but wastes memory with tons of tiny scaffolds."
	#echo "colorspace=f      		Set to true to build a SOLiD colorspace index.  Probably does not work any more."
	echo "path=<.>         		Specify the location to write the index, if you don't want it in the current working directory."
	#echo "minscaf=1        		Throw away scaffolds shorter than this when indexing."
	echo ""
	echo "Input Parameters:"
	echo "build=1          		Designate index to use.  Corresponds to the number specified when building the index."
	echo "in=<file>        		Primary reads input; required parameter."
	echo "in2=<file>      		For paired reads in two files."
	echo "qin=auto         		Set to 33 or 64 to specify input quality value ASCII offset."
	echo "interleaved=auto  		True forces paired/interleaved input; false forces single-ended mapping."
	echo "                 		If not specified, interleaved status will be autodetected from read names."
	echo "fastareadlen=500  		Break up FASTA reads longer than this.  Max is 500.  Only works for FASTA input."
	echo "fakequality=-1     		Set to a positive number 1-50 to generate fake quality strings for fasta input reads."
	#echo "parsecustom=f      		Specially process read headers from my random read generator, to determine true and false positive rates."
	echo ""
	echo "Sampling Parameters:"
	echo "reads=-1         		Set to a positive number N to only process the first N reads (or pairs), then quit.  -1 means use all reads."
	#echo "idmodulo=1       		Set to a number N to only map every Nth read (for deterministic sampling)."
	echo "samplerate=1     		Set to a number from 0 to 1 to randomly select that fraction of reads for mapping. 1 uses all reads."
	#echo "sampleseed=1  	  	Set to a positive number N set the RNG seed for sampling at the samplerate,"
	#echo "                 		or a negative number to select a random seed (for nondeterministic sampling)."
	echo "skipreads=0      		Set to a number N to skip the first N reads (or pairs), then map the rest."
	echo ""
	echo "Mapping Parameters:"
	echo "fast=f  			This flag is a macro which sets other paramters to run faster, at reduced sensitivity.  Bad for RNA-seq."
	echo "maxindel=16000   		Don't look for indels longer than this.  Lower is faster.  Set to >=100k for RNAseq with long introns like mammals."
	echo "strictmaxindel=f 		When enabled, do not allow indels longer than 'maxindel'.  By default these are not sought, but may be found anyway."
	#echo "minratio=0.56    		Fraction of max alignment score required to keep a site.  Higher is faster."
	echo "minid=0.76    			Approximate minimum alignment identity to look for.  Higher is faster and less sensitive."
	echo "minhits=1        		Minimum number of seed hits required for candidate sites.  Higher is faster."
	echo "k=13             		Kmer length of index.  Higher is faster (for large genomes), less sensitive, and uses more RAM.  Max is 15."
	echo "local=f         		Set to true to use local, rather than global, alignments.  This will soft-clip ugly ends of poor alignments."
	echo "perfectmode=f      		Allow only perfect mappings when set to true (very fast)."
	echo "semiperfectmode=f        	Allow only perfect and semiperfect (perfect except for N's in reference) mappings."
	echo "threads=auto      		(t) Set to number of threads desired.  By default, uses all cores available."
	echo "ambiguous=<best> 		(ambig) Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations)."
	echo "                 			best	(use the first best site)"
	echo "                 			toss	(consider unmapped)"
	echo "                 			random	(select one top-scoring site randomly)"
	echo "                 			all	(retain all top-scoring sites.  Does not work yet with SAM output)"
	#echo "kfilter=-1      		Set to a positive number N to require minimum N contiguous matches for a mapped read."
	echo "samestrandpairs=f  		(ssp) Specify whether paired reads should map to the same strand or opposite strands."
	echo "requirecorrectstrand=t     	(rcs) Forbid pairing of reads without correct strand orientation.  Set to false for long-mate-pair libraries."
	echo "killbadpairs=f  		(kbp) If a read pair is mapped with an inappropriate insert size or orientation, "
	echo "                		the read with the lower mapping quality is marked unmapped."
	echo "pairedonly=f  	 		(po) Treat unpaired reads as unmapped.  Thus they will be sent to 'outu' but not 'outm'."
	echo "rcompmate=f      		Reverse complement second read in each pair prior to mapping."
	echo "pairlen=32000    		Set max allowed distance between paired reads.  (insert size)=(pairlen)+(read1 length)+(read2 length)"
	echo "trim=f      		 	Quality-trim ends to Q5 before mapping.  Options are 'l' (left), 'r' (right), and 'lr' (both)."
	echo "untrim=f         		Undo trimming after mapping.  Untrimmed bases will be soft-clipped in cigar strings."
	echo "bandwidthratio=0  		(bwr) If above zero, restrict alignment band to this fraction of read length.  Faster but less accurate."
	echo ""
	echo "Output Parameters:"
	echo "outputunmapped=t     		Set to false if unmapped reads should not be printed to 'out=' target (saves time and disk space)."
	echo "out=<file>      		Write all reads to this file (unless outputunmapped=t)."
	echo "outu=<file>      		Write only unmapped reads to this file.  Does not include unmapped paired reads with a mapped mate."
	echo "outm=<file>      		Write only mapped reads to this file.  Includes unmapped paired reads with a mapped mate."
	echo "scafstats=<file>     		Write statistics on how many reads mapped to which scaffold to this file."
	echo "refstats=<file>      		Write statistics on how many reads mapped to which reference to this file (for BBSplitter)."
	echo "qualityhistogram=<file> 	(qhist) Write histogram of quality score by read location to this file."
	echo "matchhistogram=<file>	 	(mhist) Write histogram of base match, substitution, deletion, and insertion rates by read location."
	echo "inserthistogram=<file>  	(ihist) Write histogram of insert sizes (for paired reads)."
	echo "bamscript=<file>        	(bs) Write a shell script to <file> that will turn the sam output into a sorted, indexed bam file."
	echo "ordered=f        		Set to true to output reads in same order as input.  Slower and uses more memory."
	#echo "                 		Only relevant with multiple mapping threads."
	#echo "showprogress=0   		Set to a positive number N to print a '.' once per N reads processed."
	echo "cigar=t      			Set to 'f' to skip generation of cigar strings (faster)."
	#echo "                 		'none' is faster, but prevents generation of match and error rate statistics."
	echo "overwrite=f      		(ow) Allow process to overwrite existing files."
	echo "secondary=f      		Print secondary alignments."
	echo "maxsites=5       		Maximum number of total alignments to print per read.  Only relevant when secondary=t."
	echo "quickmatch=f      		Generate cigar strings more quickly.  Must be true to generate secondary site cigar strings."
	echo "keepnames=f  			Keep original names of paired reads, rather than ensuring both reads have the same name."
	echo "trimreaddescriptions=f  	(trd) Truncate read names at the first whitespace, assuming that the remaineder is a comment or description."
	echo "sam=1.3          		Set to 1.4 to write Sam version 1.4 cigar strings, with = and X instead of M."
	echo "md=f             		Write MD tags."
	echo "xs=f             		Set to 'xs=fs', 'xs=ss', or 'xs=us' to write XS tags for RNAseq using firststrand,"
	echo "                 		secondstrand, or unstranded libraries.  Needed by Cufflinks.  JGI mainly uses 'firststrand'."
	echo "stoptag=t        		Write a tag indicating read stop location, prefixed by YS:i:"
	echo "idtag=t          		Write a tag indicating percent identity, prefixed by YI:f:"
	echo "ziplevel=2          		(zl) Set to true to write a tag indicating percent identity, prefixed by YI:f:"
	echo "machineout=f         		Set to true to output statistics in machine-friendly 'key=value' format."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx                  	This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "					-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "This list is not complete.  For more information, please consult $DIR""docs/readme.txt"
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbmap "$@"
