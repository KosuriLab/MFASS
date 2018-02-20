#!/bin/bash -l
#countgc in=<infile> out=<outfile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx120m"
calcXmx () {
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
		fi
	done
}
calcXmx "$@"

countgc() {
	module unload oracle-jdk
	module load oracle-jdk/1.7_64bit
	local CMD="java -ea $z -cp $CP jgi.CountGC $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 11, 2013"
	echo ""
	echo "Description:  Counts GC content of reads or scaffolds."
	echo ""
	echo "Usage:	countgc in=<input> out=<output> format=<format>"
	echo ""
	echo "Input may be stdin or a fasta or fastq file, compressed or uncompressed."
	echo "Output (which is optional) may be stdout or a file."
	echo "format=1:	name	start	stop	A	C	G	T	N"
	echo "format=2:	name	GC"
	echo "format=4:	name	length	GC"
	echo "Note that in format 1, A+C+G+T=1 even when N is nonzero."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

if [ -z "$1" ]; then
	usage
	exit
fi

countgc "$@"
