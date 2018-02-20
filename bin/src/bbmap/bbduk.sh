#!/bin/bash -l
#bbduk in=<infile> out=<outfile>

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
		echo "This system does not have ulimit set, so max memory cannot be determined.  Attempting to use 2G." 1>&2
		echo "If this fails, please set ulimit or run this program qsubbed or from a qlogin session on Genepool." 1>&2
		y=2000
	fi
	
	mult=42;

	y=$(( ((x-20000)*mult/100)/1000 ))

	if [ $y -ge 15000 ]; then
		y=15000
	elif [ 100 -ge $y ]; then
		y=100
	fi
	
	#echo "y=$y"
	z="-Xmx${y}m"
	
	for arg in "$@"
	do
		if [[ "$arg" == -Xmx* ]]; then
			z="$arg"
		fi
	done
}
calcXmx "$@"

bbduk() {
	module unload oracle-jdk
	module unload samtools
	module load oracle-jdk/1.7_64bit
	module load pigz
	module load samtools
	local CMD="java -ea $z -cp $CP jgi.BBDukF $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "This script is designed for Genepool nodes."
	echo "Last modified December 11, 2013"
	echo ""
	echo "Description:  Compares reads to the kmers in a reference dataset, optionally allowing an edit distance."
	echo "Splits the reads into two outputs - those that match the reference, and those that don't."
	echo "Can also trim (remove) the matching parts of the reads rather than binning the reads."
	echo ""
	echo "Usage:	bbduk.sh in=<input file> out=<output file> ref=<contaminant files>"
	echo ""
	echo "Input may be stdin or a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file."
	echo "If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz"
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in=<file>        	The 'in=' flag is needed only if the input file is not the first parameter.  'in=stdin.fq' will pipe from standard in."
	echo "in2=<file>       	Use this if 2nd read of pairs are in a different file."
	echo "ref=<file,file>  	Comma-delimited list of reference files."
	echo "touppercase=f    	(tuc) Change all letters in reads and reference to upper-case."
	echo "interleaved=auto 	(int) If true, forces fastq input to be paired and interleaved."
	echo "qin=auto         	ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto."
	echo "reads=-1         	If set to a positive number, only process this many reads (or pairs), then quit."
#	echo "skipreads=0      	Ignore this many initial reads (or pairs) and process the rest."
	echo ""
	echo "Output parameters:"
	echo "out=<file>       	(outnonmatch) Write reads here that do not contain kmers matching the database.  'out=stdout.fq' will pipe to standard out."
	echo "out2=<file>      	(outnonmatch2) Use this to write 2nd read of pairs to a different file."
	echo "outmatch=<file>    	(outm or outb) Write reads here that contain kmers matching the database."
	echo "outmatch2=<file>   	(outm2 or outb2) Use this to write 2nd read of pairs to a different file."
	echo "stats=<file>     	Write statistics about which contamininants were detected."
	echo "duk=<file>       	Write statistics in duk's format."
	echo "overwrite=t      	(ow) Set to false to force the program to abort rather than overwrite an existing file."
	echo "showspeed=t      	(ss) Set to 'f' to suppress display of processing speed."
	echo "ziplevel=2       	(zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster."
	echo "fastawrap=80     	Length of lines in fasta output."
	echo "qout=auto        	ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input)."
	echo ""
	echo "Processing parameters:"
	echo "threads=auto     	(t) Set number of threads to use; default is number of logical processors."
	echo "k=28             	Kmer length used for finding contaminants.  Contaminants shorter than k will not be found.  k must be at least 1."
	echo "rcomp=t            	Look for reverse-complements of kmers in addition to forward kmers."
	echo "maskmiddle=t     	(mm) Treat the middle base of a kmer as a wildcard, to increase sensitivity in the presence of errors."
	echo "maxbadkmers=0    	(mbk) Reads with more than this many contaminant kmers will be discarded."
	echo "hammingdistance=0	(hdist) Maximum Hamming distance from ref kmers (subs only).  Memory use is proportional to (3*K)^hdist."
	echo "editdistance=0   	(edist) Maximum edit distance from ref kmers (subs and indels).  Memory use is proportional to (8*K)^edist."
	echo "forbidn=f     	  	(fn) Forbids matching of read kmers containing N.  By default, these will match a reference 'A' if hdist>0 or edist>0, to increase sensitivity."
	echo "minskip=1    	  	(mns) Force minimal skip interval when indexing reference kmers.  1 means use all, 2 means use every other kmer, etc."
	echo "maxskip=99    	  	(mxs) Restrict maximal skip interval when indexing reference kmers."
	echo "                 	Normally all are used for scaffolds<100kb, but with longer scaffolds, up to K-1 are skipped."
	echo "removeifeitherbad=t	(rieb) Paired reads get sent to 'outmatch' if either is match (or either is trimmed shorter than minlen).  Set to false to require both."
	echo ""
	echo "Trimming parameters:"
	echo "ktrim=f          	Trim reads to remove bases matching reference kmers."
	echo "                 	Values: f (don't trim), r (trim right end), l (trim left end), n (convert to N instead of trimming)."
	echo "                 	Any single non-whitespace character other than t, f, r, l, n: convert to that symbol rather than trimming."
	echo "useshortkmers=f  	(usk) Look for shorter kmers at read tips (only for k-trimming).  Enabling this will disable maskmiddle."
	echo "mink=4           	Minimum length of short kmers.  Setting this automatically sets useshortkmers=t."
	echo "qtrim=f          	Trim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=4           	Trim quality threshold."
	echo "minlength=20     	(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter."
	echo "minavgquality=0  	(maq) Reads with average quality (before trimming) below this will be discarded."
	echo "otm=f            	(outputtrimmedtomatch) Output reads trimmed to shorter than minlength to outm rather than discarding."
	echo ""
#	echo "Other parameters:"
#	echo "array=t    	  	Use HashArray data structure."
#	echo "forest=f    	  	Use HashForest data structure."
#	echo "table=f    	  	Use KmerTable data structure."
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "There is a changelog at /global/projectb/sandbox/gaag/bbtools/docs/changelog_bbduk.txt"
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

if [ -z "$1" ]; then
	usage
	exit
fi

bbduk "$@"
