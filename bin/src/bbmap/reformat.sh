#!/bin/bash -l
#reformat in=<infile> out=<outfile>

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

function reformat() {
	module unload oracle-jdk
	module unload samtools
	module load oracle-jdk/1.7_64bit
	module load pigz
	module load samtools
	local CMD="java -ea $z -cp $CP jgi.ReformatReads $@"
	echo $CMD >&2
	$CMD
}

function usage(){
	echo "Written by Brian Bushnell"
	echo "Last modified December 11, 2013"
	echo ""
	echo "Description:  Reformats reads to change ASCII quality encoding, interleaving, file format, or compression format."
	echo "Optionally performs additional functions such as quality trimming, subsetting, and subsampling."
	echo "Supports sam, fastq, fasta, fasta+qual, gzip, zip, and bz2."
	echo ""
	echo "Usage:  reformat.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2>"
	echo ""
	echo "in2 and out2 are for paired reads and are optional."
	echo "If input is paired and there is only one output file, it will be written interleaved."
	echo "Other parameters and their defaults:"
	echo ""
	echo "ow=f         		(overwrite) Overwrites files that already exist."
	echo "zl=4            	(ziplevel) Set compression level, 1 (low) to 9 (max)."
	echo "int=f		 	(interleaved) Determines whether INPUT file is considered interleaved."
	echo "fastawrap=100    	Length of lines in fasta output."
	echo "fastareadlen=0   	Set to a non-zero number to break fasta files into reads of at most this length."
	echo "minscaf=1       	Ignore fasta reads shorter than this."
	echo "tuc=f    		(touppercase) Change lowercase letters in reads to uppercase."
	echo "rcm=f    		(rcompmate) Reverse-complement read 2 in pairs."
	echo "qin=auto         	ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto."
	echo "qout=auto        	ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input)."
	echo "qfin=<.qual file>	Read qualities from this qual file, for the reads coming from 'in=<fasta file>'"
	echo "qfin2=<.qual file>	Read qualities from this qual file, for the reads coming from 'in2=<fasta file>'"
	echo "qfout=<.qual file>	Write qualities from this qual file, for the reads going to 'out=<fasta file>'"
	echo "qfout2=<.qual file>	Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'"
	echo "verifyinterleaved=f	(vint) When true, checks a file to see if the names look paired.  Prints an error message if not."
	echo "tossbrokenreads=f	(tbr) Discard reads that have different numbers of bases and qualities.  By default this will be detected and cause a crash."
	echo ""
	echo "Sampling parameters:"
	echo "reads=-1 		Set to a positive number to only process this many INPUT reads (or pairs), then quit."
	echo "samplerate=1  		Randomly output only this fraction of reads; 1 means sampling is disabled."
	echo "sampleseed=-1 		Set to a positive number to use that prng seed for sampling (allowing deterministic sampling)."
	echo "samplereadstarget=0	(srt) Exact number of OUTPUT reads (or pairs) desired."
	echo "samplebasestarget=0	(sbt) Exact number of OUTPUT bases desired."
	echo "                      Important: srt/sbt flags should not be used with stdin, samplerate, qtrim, minlength, or minavgquality."
	echo ""
	echo "Trimming parameters:"
	echo "qtrim=f          	Trim read ends to remove bases with quality below minq."
	echo "                 	Values: t (trim both ends), f (neither end), r (right end only), l (left end only)."
	echo "trimq=4           	Trim quality threshold."
	echo "minlength=20     	(ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter."
	echo "minavgquality=0     	(maq) Reads with average quality (before trimming) below this will be discarded."
	echo "forcetrimleft=0     	(ftl) If nonzero, trim left bases of the read to this position (exclusive, 0-based)."
	echo "forcetrimright=0    	(ftr) If nonzero, trim right bases of the read after this position (exclusive, 0-based)."
	echo ""
	echo "Sam file reformatting options.  Note that most of these will require an indexed reference."
	echo "build=<integer> 	Assign a genome's build id.  You can index like this: bbmap.sh ref=<file> build=1"
	echo "sam=1.3         	Set to 1.4 to write Sam version 1.4 cigar strings, with = and X instead of M."
	echo "md=f            	Set to true to write MD tags."
	echo "xs=f            	Set to 'ss', 'fs', or 'us' to write XS tags for RNAseq using secondstrand, firststrand,"
	echo "                	or unstranded libraries.  Needed by Cufflinks.  JGI mainly uses 'firststrand'."
	echo "stoptag=t       	Set to true to write a tag indicating read stop location, prefixed by YS:i:"
	echo "idtag=t         	Set to true to write a tag indicating percent identity, prefixed by YI:f:"
	echo ""
	echo "Java Parameters:"
	echo "-Xmx       		This will be passed to Java to set memory usage, overriding the program's automatic memory detection."
	echo "				-Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory."
	echo ""
	echo "Supported input formats are fastq, fasta, fast+qual, scarf, and bread (BBMap's native format)"
	echo "Supported output formats are fastq, fasta, fast+qual, bread, sam, and bam (bam only if samtools is installed)"
	echo "Supported compression formats are gz, zip, and bz2"
	echo "To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'"
	echo "To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'"
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
}

if [ -z "$1" ]; then
	usage
	exit
fi

reformat "$@"
