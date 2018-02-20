#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

function tf() {
	module load oracle-jdk/1.7_64bit
	local CMD="java -ea -Xmx120m -cp $CP jgi.GetReads $@"
	echo $CMD >&2
	$CMD
}

function usage(){
	echo "Selects reads with designated numeric IDs."
	echo "Last modified December 18, 2013."
	echo ""
	echo "Usage:	getreads.sh in=<file> id=<number,number,number...> out=<file>"
	echo ""
	echo "The first read (or pair) has ID 0, the second read (or pair) has ID 1, etc."
	echo ""
	echo "Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems."
	echo ""
}

if [ -z "$1" ]; then
	usage
	exit
fi

tf "$@"