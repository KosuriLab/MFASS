#!/bin/bash -l
#stats in=<infile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

stats() {
	module unload oracle-jdk
	module load oracle-jdk/1.7_64bit
	local CMD="java -ea -Xmx120m -cp $CP jgi.MakeLengthHistogram $@"
#	echo $CMD >&2
	$CMD
}

usage(){
	echo "Generates a length histogram of input reads."
	echo "Written by Brian Bushnell"
	echo "Last modified December 11, 2013"
	echo ""
	echo "Usage:	readlength.sh in=<input file>"
	echo ""
	echo "in=<file>        	The 'in=' flag is needed only if the input file is not the first parameter.  'in=stdin.fq' will pipe from standard in."
	echo "in2=<file>       	Use this if 2nd read of pairs are in a different file."
	echo "out=<file>       	Write the histogram to this file.  Default is stdout."
	echo "bin=10        	  	Set the histogram bin size.  Default is 10."
	echo "max=4000         	Set the histogram's max readlength bin.  Default is 4000."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

if [ -z "$1" ]; then
	usage
	exit
fi

stats "$@"
