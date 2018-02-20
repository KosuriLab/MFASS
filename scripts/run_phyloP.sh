infile=$1
outfile=$2

# create new empty file
echo > $outfile
echo 
for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    grep -w chr$i $infile > chr_infile.bed
    bigWigAverageOverBed ../../ref/phyloP/chr$i.phyloP100way.bigWig chr_infile.bed "tmp.bed" 
    cat tmp.bed >> $outfile
done

rm chr_infile.bed
rm tmp.bed