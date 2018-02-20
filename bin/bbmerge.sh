#!/bin/bash -l
#merge in=<infile> out=<outfile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx200m"
calcXmx () {
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
		fi
	done
}
calcXmx "$@"

function merge() {
	module unload oracle-jdk
	module unload samtools
	module load oracle-jdk/1.7_64bit
	module load pigz
	module load samtools
	local CMD="java -ea $z -cp $CP jgi.MateReadsMT $@"
	echo $CMD >&2
	$CMD
}

function usage(){
	echo "BBMerge v1.4"
	echo "This script is designed for Genepool nodes."
	echo "Last modified February 12, 2014"
	echo ""
	echo "Description:  Merges paired reads into single reads by overlap detection."
	echo "With sufficient coverage, can also merge nonoverlapping reads using gapped kmers."
	echo ""
	echo "Usage:	bbmerge.sh in=<input> out=<merged reads> outbad=<unmerged reads>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in2=null  		Second input file for paired reads."
	echo "extra=null		Additional files to use for input (generating hash table) but not for output."
	echo "interleaved=auto  	May be set to true or false to force the input read file to override autodetection of the input file as paired interleaved."
	echo "reads=-1  		Only process this number of reads, then quit (-1 means all)."
	echo ""
	echo "Output parameters:"
	echo "out=<file>        	File for merged reads."
	echo "outbad=<file>       	File for unmerged reads."
	echo "outinsert=<file>    	File list of read names and their insert sizes."
	echo "hist=null    		Insert length histogram output file."
	echo "ziplevel=2   		Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo ""
	echo "Trimming parameters:"
	echo "qtrim=f          	Trim read ends to remove bases with quality below minq.  Performed BEFORE merging."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=10           	Trim quality threshold."
	echo "minlength=20     	(ml) Reads shorter than this after trimming (before merging) will be discarded.  Pairs will be discarded only if both are shorter."
	echo "trimonfailure=t     	(tof) If detecting insert size by overlap fails, the reads will be trimmed and this will be re-attempted."
	echo ""
	echo "Other parameters:"
	echo "join=t			Create merged reads.  If set to false, you can still generate an insert histogram."
	echo "useoverlap=t		Attempt merge based on paired ead overlap."
	echo "minoverlapbases=12	Minimum number of overlapping bases to merge reads."
	echo "mininsert=0		Reads with insert sizes less than this (after merging) will be discarded."
	echo "gap=null      	   	Sets gap size for merging via gapped kmers."
	echo "                  	'gap=50;50,100' would run one pass with a gap size of 50 and another with both 50 and 100."
	echo "                  	This script sets memory appropriately for ungapped merging only, though."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

if [ -z "$1" ]; then
	usage
	exit
fi

merge "$@"
