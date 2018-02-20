#!/bin/bash -l
#stats in=<infile>

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

stats() {
	module unload oracle-jdk
	module load oracle-jdk/1.7_64bit
	local CMD="java -ea $z -cp $CP jgi.AssemblyStatsWrapper $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "Last modified December 11, 2013."
	echo ""
	echo "Description:  Runs stats.sh on multiple assemblies to produce one ouput line per file."
	echo ""
	echo "Usage:	statswrapper.sh in=<input file>"
	echo ""
	echo "Parameters:"
	echo "in=<file>      	Specify the input fasta file, or stdin.  For multiple files a, b, and c: 'statswrapper.sh in=a in=b in=c'."
	echo "               	'in=' may be omitted if this is the first arg, and asterisks may be used; e.g. statswrapper.sh *.fa"
	echo "gc=<file>      	Writes ACGTN content per scaffold to a file."
	echo "gchist=<file>  	Filename to output scaffold gc content histogram."
	echo "gcbins=<200>   	Number of bins for gc histogram."
	echo "n=<10>         	Number of contiguous Ns to signify a break between contigs."
	echo "k=<13>         	Estimate memory usage of BBMap with this kmer length."
	echo "minscaf=<0>    	Ignore scaffolds shorter than this."
	echo "n_=<t>         	This flag will prefix the terms 'contigs' and 'scaffolds' with 'n_' in formats 3-6."
	echo ""
	echo "format=<1 through 6>	Format of the stats information."
	echo "	format=1 uses variable units like MB and KB, and is designed for compatibility with existing tools."
	echo "	format=2 uses only whole numbers of bases, with no commas in numbers, and is designed for machine parsing."
	echo "	format=3 outputs stats in 2 rows of tab-delimited columns: a header row and a data row."
	echo "	format=4 is like 3 but with scaffold data only."
	echo "	format=5 is like 3 but with contig data only."
	echo "	format=6 is like 3 but the header starts with a #."
	echo ""
	echo "gcformat=<1 or 2>	Select GC output format."
	echo "	gcformat=1:	name	start	stop	A	C	G	T	N	GC"
	echo "	gcformat=2:	name	GC"
	echo "	Note that in gcformat 1, A+C+G+T=1 even when N is nonzero."
	echo ""
}

if [ -z "$1" ]; then
	usage
	exit
fi

stats "$@"
