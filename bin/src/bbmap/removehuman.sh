#!/bin/bash -l
#removehuman in=<infile> out=<outfile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

function removehuman() {
	module unload oracle-jdk
	#module unload samtools
	module load oracle-jdk/1.7_64bit
	module load pigz
	#module load samtools
	local CMD="java -ea -Xmx23g -cp $CP align2.BBMap minratio=0.75 maxindel=20 bwr=0.18 bw=20 quickmatch minhits=2 outputunmapped=f path=/global/projectb/sandbox/gaag/bbtools/hg19 pigz unpigz $@"
	echo $CMD >&2
	$CMD
}

function usage(){
	echo "This script is designed for Genepool nodes.  Requires at least 24GB RAM."
	echo "Last modified January 29, 2014"
	echo ""
	echo "Description:  Removes all reads that map to the human genome with at least 88% identity."
	echo ""
	echo "Usage:	removehuman.sh in=<input file> outu=<clean output file>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a fasta, fastq, or sam file."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "overwrite=t      	(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "interleaved=auto 	(int) If true, forces fastq input to be paired and interleaved."
	#echo "kfilter=47           	Require at least this many contiguous exact matches to the reference."
	echo "trim=t           	Trim read ends to remove bases with quality below minq."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "untrim=t           	Undo the trimming after mapping."
	echo "minq=4           	Trim quality threshold."
	echo "ziplevel=2       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "outm=<file>       	File to output the reads that mapped to human."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

if [ -z "$1" ]; then
	usage
	exit
fi

removehuman "$@"
