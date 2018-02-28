# get all SNVs from ExAC annotation, this takes awhile
zgrep 'SNV' ../../ref/snv/ExAC.r0.3.1.sites.vep.vcf.gz > ../../ref/snv/ExAC_SNVs.vcf.gz

# reformat to bed format, add ID and strand column. Easier to just do this in Python
python snv_vcf2bed.py ../../ref/snv/ExAC_SNVs.vcf.gz ../../ref/snv/ExAC_SNVs.bed.gz

# intersect with bed file containing exonic and intronic regions to get location
# of SNVs, print the original feature a and original feature b
intersectBed -a ../../ref/snv/ExAC_SNVS.bed.gz -b ../../ref/exon_intron_all.bed.gz -wa -wb | 
gzip > ../../ref/snv/ExAC_SNVs_genomic_overlap.bed.gz