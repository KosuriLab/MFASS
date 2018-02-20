#!/bin/bash -l
#printtime in=<infile> out=<outfile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

function printtime() {
	module load oracle-jdk/1.7_64bit
	local CMD="java -ea -Xmx8m -cp $CP align2.PrintTime $@"
	echo $CMD >&2
	$CMD
}

function usage(){
	echo "Prints time elapsed since last called on the same file."
	echo "Written by Brian Bushnell"
	echo "Last modified October 24, 2013"
	echo ""
	echo "Usage:	printtime.sh <filename>"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

if [ -z "$1" ]; then
	usage
	exit
fi

printtime "$@"
