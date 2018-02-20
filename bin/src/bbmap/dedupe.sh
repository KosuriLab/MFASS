#!/bin/bash -l
#khist in=<infile> out=<outfile>

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
z2="-Xms1g"
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
	z2="-Xms${y}g"
	
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
		elif [[ "$arg" == -Xms* ]]; then
			z2="$arg"
		fi
	done
}
calcXmx "$@"

dedupe() {
	module unload oracle-jdk
	module load oracle-jdk/1.7_64bit
	module load pigz
	local CMD="java -ea $z $z2 -cp $CP jgi.Dedupe $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "This script is designed for Genepool nodes."
	echo "Last modified December 11, 2013"
	echo ""
	echo "Description:  Accepts one or more files containing sets of sequences (reads or scaffolds)."
	echo "Removes duplicate sequences, which may be specified to be exact matches, subsequences, or sequences within some percent identity."
	echo "Can also find overlapping sequences and group them into clusters."
	echo ""
	echo "Usage:	dedupe.sh in=<file or stdin> out=<file or stdout>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file.  With no output parameter, data will be written to stdout."
	echo "If 'out=null', there will be no output, but statistics will still be printed."
	echo "You can also use 'dedupe <infile> <outfile>' without the 'in=' and 'out='."
	echo ""
	echo "I/O Parameters"
	echo ""
	echo "in=<file,file>   	A single file or a comma-delimited list of files."
	echo "out=<file>       	Destination for all output contigs."
	echo "threads=auto          (t) Set number of threads to use; default is number of logical processors."
	echo "overwrite=t           (ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "showspeed=t           (ss) Set to 'f' to suppress display of processing speed."
	echo "minscaf=0             (ms) Ignore contigs/scaffolds shorter than this."
	echo "interleaved=auto      If true, forces fastq input to be paired and interleaved."
	echo "ziplevel=2            Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "storename=t           (sn) Store scaffold names (set false to save memory)."
	echo "storequality=t        (sq) Store quality values for fastq assemblies (set false to save memory)."
	echo "uniquenames=t         (un) Ensure all output scaffolds have unique names.  Uses more memory."
	echo "sort=f                Sort output by scaffold length (otherwise it will be random)."
	echo "                      'a' for ascending, 'd' for descending, 'f' for false (no sorting)."
	echo ""
	echo "Processing Parameters"
	echo ""
	echo "absorbrc=t            (arc) Absorb reverse-complements as well as normal orientation."
	echo "absorbmatch=t         (am) Absorb exact matches of contigs."
	echo "absorbcontainment=t   (ac) Absorb full containments of contigs."
#	echo "absorboverlap=f       (ao) Absorb (merge) non-contained overlaps of contigs (TODO)."
	echo "findoverlap=f         (fo) Find overlaps between contigs (containments and non-containments).  Necessary for clustering."
	echo "cluster=f             (c) Group overlapping contigs into clusters."
	echo ""
	echo "fixmultijoins=t       (fmj) Remove redundant overlaps between the same two contigs."
	echo "removecycles=t        (rc) Remove all cycles so clusters form trees."
	echo "renameclusters=f      (rnc) Rename contigs to indicate which cluster they are in."
	echo "cc=t                  (canonicizeclusters) Flip contigs so clusters have a single orientation."
	echo "fcc=f                 (fixcanoncontradictions) Truncate graph at nodes with canonization disputes."
	echo "foc=f                 (fixoffsetcontradictions) Truncate graph at nodes with offset disputes."
	echo ""
	echo "Overlap Detection Parameters"
	echo ""
	echo "exact=t               (ex) Only allow exact symbol matches.  When false, an 'N' will match any symbol."
	echo "touppercase=f         (tuc) Change all input bases to upper case."
	echo "maxsubs=0             (s) Allow up to this many mismatches (substitutions only, no indels).  May be set higher than maxedits."
	echo "maxedits=0            (e) Allow up to this many edits (subs or indels).  Higher is slower."
	echo "minidentity=100       (mid) Absorb contained sequences with percent identity of at least this (includes indels)."
	echo "minlengthpercent=0    (mlp) Smaller contig must be at least this percent of larger contig's length to be absorbed."
	echo "minoverlappercent=0   (mop) Overlap must be at least this percent of smaller contig's length to cluster and merge."
	echo "minoverlap=200        (mo) Overlap must be at least this long to cluster and merge."
	echo "depthratio=0          (dr) When non-zero, overlaps will only be formed between reads with a depth ratio of at most this."
	echo "                      Should be above 1.  Depth is determined by parsing the read names; this information can be added"
	echo "                      by running KmerNormalize (khist.sh, bbnorm.sh, or ecc.sh) with the flag 'rename'"
	echo "k=31                  Seed length used for finding containments and overlaps.  Anything shorter than k will not be found."
#	echo "numaffixmaps=1        (nam) Set to 2 to index two prefixes and suffixes per contig."
#	echo "ignoreaffix1=f        (ia1) Ignore first affix (for testing)."
#	echo "storesuffix=f         (ss) Store suffix as well as prefix.  Automatically set to true when doing inexact matches."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo "-Xms       		If you use the -Xmx flag, also set -Xms to the same value."
	echo ""
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

if [ -z "$1" ]; then
	usage
	exit
fi

dedupe "$@"
