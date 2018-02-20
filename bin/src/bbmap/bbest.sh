#!/bin/bash -l
#bbest in=<infile> out=<outfile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

function bbest() {
	module unload oracle-jdk
	module unload samtools
	module load oracle-jdk/1.7_64bit
	module load samtools
	module load pigz
	local CMD="java -ea -Xmx64m -cp $CP jgi.SamToEst $@"
#	echo $CMD >&2
	$CMD
}

function usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 5, 2013"
	echo ""
	echo "Description:  Calculates EST (expressed sequence tags) capture by an assembly from a sam file."
	echo "Designed to use BBMap output generated with these flags: k=13 maxindel=100000 customtag ordered"
	echo ""
	echo "Usage:	bbest.sh in=<sam file> out=<stats file>"
	echo ""
	echo "Parameters:"
	echo "in=<file>      	Specify a sam file (or stdin) containing mapped ests."
	echo "out=<file>     	Specify the output stats file (default is stdout)."
	echo "ref=<file>     	Specify the reference file (optional)."
	echo "est=<file>     	Specify the est fasta file (optional)."
	echo "fraction=<0.98>      	Min fraction of bases mapped to ref to be considered 'all mapped'."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbest "$@"
