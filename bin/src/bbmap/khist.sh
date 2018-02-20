#!/bin/bash -l
#khist in=<infile> out=<outfile>

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

khist() {
	module unload oracle-jdk
	module load oracle-jdk/1.7_64bit
	module load pigz
	local CMD="java -ea $z -cp $CP jgi.KmerNormalize bits=32 ecc=f passes=1 keepall dr=f prefilter hist=stdout minprob=0 minqual=0 mindepth=0 minkmers=1 hashes=3 $@"
	echo $CMD >&2
	$CMD
}

usage(){
	echo "This script is designed for Genepool nodes."
	echo "Last modified February 12, 2014"
	echo ""
	echo "Description:  Generates a histogram of kmer counts for the input reads or assemblies."
	echo "Can also normalize, error-correct, and/or bin reads by kmer depth."
	echo ""
	echo "Usage:	khist.sh in=<input> out=<reads to keep> outt=<reads to toss> hist=<histogram output>"
	echo ""
	echo "Input may be a fasta, fastq, or sam file, compressed or uncompressed."
	echo "Output may be stdout or a file.  'out' and 'hist' are both optional."
	echo ""
	echo "Optional parameters (and their defaults)"
	echo ""
	echo "Input parameters:"
	echo "in2=null		Second input file for paired reads"
	echo "extra=null		Additional files to use for input (generating hash table) but not for output"
	echo "fastareadlen=2^31	Break up FASTA reads longer than this.  Can be useful when processing scaffolded genomes"
	echo "tablereads=-1		Use at most this many reads when building the hashtable (-1 means all)"
	echo "kmersample=1		Process every nth kmer, and skip the rest"
	echo "readsample=1		Process every nth read, and skip the rest"
	echo "interleaved=auto	May be set to true or false to force the input read file to ovverride autodetection of the input file as paired interleaved."
	echo ""
	echo "Output parameters:"
	echo "out=<file>        	File for normalized reads"
	echo "outt=<file>      	(outtoss) File for reads that were excluded from primary output"
	echo "reads=-1		Only process this number of reads, then quit (-1 means all)"
	echo "sampleoutput=t		Use sampling on output as well as input (not used if sample rates are 1)"
	echo "keepall=f		Set to true to keep all reads (e.g. if you just want error correction)."
	echo "zerobin=f		Set to true if you want kmers with a count of 0 to go in the 0 bin instead of the 1 bin in histograms."
	echo "             		Default is false, to prevent confusion about how there can be 0-count kmers."
	echo "             		The reason is that based on the 'minq' and 'minprob' settings, some kmers may be excluded from the bloom filter."
	echo "tmpdir=$TMPDIR  	This will specify a directory for temp files (only needed for multipass runs).  If null, they will be written to the output directory."
	echo "usetempdir=t    	Allows enabling/disabling of temporary directory; if disabled, temp files will be written to the output directory."
	echo ""
	echo "Hashing parameters:"
	echo "k=31			Kmer length (values under 32 are most efficient, but arbitrarily high values are supported)"
	echo "bits=32			Bits per cell in bloom filter; must be 2, 4, 8, 16, or 32.  Maximum kmer depth recorded is 2^cbits.  Automatically reduced to 16 in 2-pass."
	echo " 			Large values decrease accuracy for a fixed amount of memory, so use the lowest number you can that will still capture highest-depth kmers."
	echo "hashes=3		Number of times each kmer is hashed and stored.  Higher is slower."
	echo "  			Higher is MORE accurate if there is enough memory, and LESS accurate if there is not enough memory."
	echo "prefiliter=f		True is slower, but generally more accurate; filters out low-depth kmers from the main hashtable.  The prefilter is more memory-efficient because it uses 2-bit cells."
	echo "prehashes=2		Number of hashes for prefilter."
	echo "buildpasses=1		More passes can sometimes increase accuracy by iteratively removing low-depth kmers"
	echo "minq=6			Ignore kmers containing bases with quality below this"
	echo "minprob=0.5		Ignore kmers with overall probability of correctness below this"
	echo "threads=X		Spawn exactly X hashing threads (default is number of logical processors).  Total active threads may exceed X due to I/O threads."
	echo "rdk=t			(removeduplicatekmers) When true, a kmer's count will only be incremented once per read pair, even if that kmer occurs more than once."
	echo ""
	echo "Normalization parameters:"
	echo "fixspikes=f	  	(fs) Do a slower, high-precision bloom filter lookup of kmers that appear to have an abnormally high depth due to collisions."
	echo "target=40     		(tgt) Target normalization depth.  NOTE:  All depth parameters control kmer depth, not read depth."
	echo "         		For kmer depth Dk, read depth Dr, read length R, and kmer size K:  Dr=Dk*(R/(R-K+1))"
	echo "maxdepth=-1		(max) Reads will not be downsampled when below this depth, even if they are above the target depth."        		
	echo "mindepth=6        	(min) Kmers with depth below this number will not be included when calculating the depth of a read."
	echo "minkmers=15		(mgkpr) Reads must have at least this many kmers over min depth to be retained.  Aka 'mingoodkmersperread'."
	echo "percentile=54.0	(dp) Read depth is by default inferred from the 54th percentile of kmer depth, but this may be changed to any number 1-100."
	echo "deterministic=t	(dr) Generate random numbers deterministically to ensure identical output between multiple runs.  May decrease speed with a huge number of threads."
	echo "passes=2		(p) 1 pass is the basic mode.  2 passes (default) allows greater accuracy, error detection, better contol of output depth."
	echo ""
	echo "Error detection parameters:"
	echo "hdp=90.0 		(highdepthpercentile) Position in sorted kmer depth array used as proxy of a read's high kmer depth."
	echo "ldp=25.0    		(lowdepthpercentile) Position in sorted kmer depth array used as proxy of a read's low kmer depth."
	echo "tossbadreads=f		(tbr) Throw away reads detected as containing errors.  Only controls behavior of final pass."
	echo "errordetectratio=125	(edr) Reads with a ratio of at least this much between their high and low depth kmers will be classified as error reads."
	echo "highthresh=12		(ht) Threshold for high kmer.  A high kmer at this or above are considered non-error."
	echo "lowthresh=3		(lt) Threshold for low kmer.  Kmers at this and below are always considered errors."
	echo ""
	echo "Error correction parameters:"
	echo "ecc=f     		Set to true to correct errors."
	echo "ecclimit=3		Correct up to this many errors per read.  If more are detected, the read will remain unchanged."
	echo "errorcorrectratio=140	(ecr) Adjacent kmers with a depth ratio of at least this much between will be classified as an error."
	echo "echighthresh=22   	(echt) Threshold for high kmer.  A kmer at this or above may be considered non-error."
	echo "eclowthresh=2		(eclt) Threshold for low kmer.  Kmers at this and below are considered errors."
	echo "eccmaxqual=127		Do not correct bases with quality above this value."
	echo "meo=f			(markErrorsOnly) Marks errors by reducing quality value of suspected errors; does not correct anything."
	echo "mue=t			(markUncorrectableErrors) Marks errors only on uncorrectable reads; requires 'ecc=t'."
	echo ""	
	echo "Depth binning parameters:"
	echo "lowbindepth=10		(lbd) Cutoff for low depth bin."
	echo "highbindepth=80	(hbd) Cutoff for high depth bin."
	echo "outlow=<file>		Pairs in which both reads have a median below lbd go into this file."
	echo "outhigh=<file>		Pairs in which both reads have a median above hbd go into this file."
	echo "outmid=<file>		All other pairs go into this file."
	echo ""
	echo "Histogram parameters:"
	echo "hist=<file>		Specify a file to write the input kmer depth histogram"
	echo "histout=<file>		Specify a file to write the output kmer depth histogram"
	echo "pzc=f			(printzerocoverage) Print lines in the histogram with zero coverage."
	echo "histlen=1048576   	Max kmer depth displayed in histogram.  Also affects statistics displayed, but does not affect normalization."
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

khist "$@"
