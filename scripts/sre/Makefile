all: format_ref rescore_ref clean_data

format_ref: ../../ref/sre/sre_ref_formatted_converted.txt

rescore_ref: ../../ref/sre/sre_ref_rescored.txt

clean_data: ../../processed_data/sre/sre_data_clean.txt

Rscript = /Library/Frameworks/R.framework/Resources/Rscript


../../ref/sre/sre_ref_formatted_converted.txt: ../../ref/sre/sre_ref.txt sre_format_ref.R sre_ref_liftover.sh
	@ echo "Formatting reference..."
	@ $(Rscript) sre_format_ref.R


../../ref/sre/sre_ref_rescored.txt: ../../ref/sre/sre_ref_formatted_converted.txt ../score_lib.py
	@ echo "Re-scoring reference file..."
	@ cd ../splicemod_src && \
	python ../score_lib.py $< seq ./data/motifs/Ke2011/ \
	../../ref/Rosenberg_2015/exonic_mer6_scores.series 2 $@

# process sequencing data
../../processed_data/sre/sre_data_clean.txt: ../../ref/sre/sre_ref_formatted_converted.txt sre_format_ref.R
	@ echo "Cleaning up data..."
	@ $(Rscript) sre_data_clean.R
