# Generate bed file of all genomic intron and exon regions

echo "Extracting non-overlapping exon coordinates..."
gzcat ../../ref/gencode.v26.annotation.gtf.gz | 
# extract chromosome, start, end coordinates from a tab-separated file if the
# third field is an exon
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' |
# use mergeBed, make sure coordinates are sorted first
sortBed |
mergeBed -i - | 
awk 'BEGIN{OFS="\t";} {print $1,$2,$3,"exon"}' | # add dummy ID column
gzip > ../../ref/gencode_v26_exon_merged.bed.gz

# to define intronic regions, subtract exonic region from genic region
echo "Calculating intronic regions..."
gzcat ../../ref/gencode.v26.annotation.gtf.gz | 
awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' |
sortBed |
subtractBed -a stdin -b ../../ref/gencode_v26_exon_merged.bed.gz |
awk 'BEGIN{OFS="\t";} {print $1,$2,$3,"intron"}' | # dummy ID column
gzip > ../../ref/gencode_v26_intron.bed.gz

# combine intronic bed file with exon bed file
gzcat ../../ref/gencode_v26_exon_merged.bed.gz ../../ref/gencode_v26_intron.bed.gz | 
gzip > ../../ref/exon_intron_all.bed.gz