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

randomreads() {
	module unload oracle-jdk
	module load oracle-jdk/1.7_64bit
	module load pigz
	local CMD="java -ea $z -cp $CP align2.RandomReads3 build=1 $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "This script is designed for Genepool nodes.  It will by default attempt to use all memory on the target node."
	echo "Last modified December 11, 2013."
	echo ""
	echo "Description:  Generates random synthetic reads from a reference genome.  A read's name indicates its genomic location."
	echo "Allows precise customization of things like insert size and synthetic mutation type, sizes, and rates."
	echo ""
	echo "Usage:	randomreads.sh ref=<reference fasta> out=<output file> minlen=<min length> maxlen=<max length> reads=<number of reads>"
	echo ""
	echo "Optional parameters:"		
	echo "paired=<false>        	Set to true for paired reads."
	echo "interleaved=<false>    	Set to true if paired output is interleaved (rather than in two files)."
	echo "build=<1>             	If multiple references will be used when running in the same working directory, each needs a unique build ID."
	echo "replacenoref=<false>		Set to true to replace N in the reference sequence with random letters."
	echo "seed=<number>			Use this to set the random number generator seed; use -1 for a random seed."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       			This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "					-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
}

if [ -z "$1" ]; then
	usage
	exit
fi

randomreads "$@"
