#!/bin/bash -l
#mapPacBio in=<infile> out=<outfile> ref=<reference>

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
			mult=84
			#echo "ram is 1000g+"
		elif [ $x -ge 500000000 ]; then
			mult=84
			#echo "ram is 500g+"
		elif [ $x -ge 250000000 ]; then
			mult=84
			#echo "ram is 250g+"
		elif [ $x -ge 144000000 ]; then
			mult=84
			#echo "ram is 144g+"
		elif [ $x -ge 120000000 ]; then
			mult=84
			#echo "ram is 120g+"
		elif [ $x -ge 40000000 ]; then
			mult=80
			#echo "ram is 40g+"
		else
			mult=84
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

mapPacBio() {
	module unload oracle-jdk
	module unload samtools
	module load oracle-jdk/1.7_64bit
	module load pigz
	module load samtools
	local CMD="java -ea $z -cp $CP align2.BBMapPacBio build=1 overwrite=true minratio=0.40 match=long fastareadlen=6000 dprr=false ambiguous=best minscaf=100 startpad=10000 stoppad=10000 midpad=6000 $@"
	echo $CMD >&2
	$CMD
}

usage(){
	bash "$DIR"bbmap.sh
}

if [ -z "$1" ]; then
	usage
	exit
fi

mapPacBio "$@"
