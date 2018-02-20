# Runs UCSC liftOver command to convert coordinates from hg19 to hg38
# called in splicemod_format_ref.R and depends on file produced in first chunk

liftOver <(cat ../../ref/sre/sre_ref_formatted.txt | \
	sed 's/\"//g' | \
	awk '{OFS="\t"; print "chr"$4,$5,$6,$1}' | grep -v 'NA' | sort | uniq) \
  ../../ref/hg19ToHg38.over.chain \
  ../../processed_data/sre/sre_ref_liftover.bed \
  ../../processed_data/sre/sre_unlifted.bed
