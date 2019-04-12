all: format_ref clean_data add_annot run_HAL run_SPANR pr_curve

format_ref: ../../ref/snv/snv_ref_formatted_converted_original_seq.txt

rescore_ref: ../../ref/snv/snv_ref_rescored.txt

clean_data: ../../processed_data/snv/snv_data_clean.txt

add_annot: ../../processed_data/snv/snv_func_annot.txt

run_HAL: ../../processed_data/snv/snv_HAL_scores.txt

run_SPANR: ../../processed_data/snv/snv_SPANR_scores.txt

pr_curve: ../../processed_data/snv/snv_models_pr_curves_all.txt

fig2: ../../figs/snv/fig2C_snv_density.tiff


Rscript = /Library/Frameworks/R.framework/Resources/Rscript


../../ref/snv/snv_ref_formatted_converted_original_seq.txt: ../../ref/snv/snv_ref_all.txt
	@ echo "Formatting reference..."
	@ $(Rscript) snv_format_ref.R


../../ref/snv/snv_ref_rescored.txt: ../../ref/snv/snv_ref_formatted_converted_original_seq.txt
	@ echo "Re-scoring reference file..."
	@ cd ../splicemod_src && \
	python ../score_lib.py $< original_seq ./data/motifs/Ke2011/ \
	../../ref/Rosenberg_2015/exonic_mer6_scores.series 30 $@


../../processed_data/snv/snv_data_clean.txt: ../../ref/snv/snv_ref_formatted_converted_original_seq.txt snv_data_clean.R
	@ echo "Cleaning up data..."
	@ $(Rscript) snv_data_clean.R


../../processed_data/snv/snv_func_annot.txt: ../../processed_data/snv/snv_data_clean.txt snv_add_annot.R
	@ echo "Adding functional annotation..."
	@ $(Rscript) snv_add_annot.R


../../processed_data/snv/snv_HAL_scores.txt: ../../processed_data/snv/snv_data_clean.txt snv_run_HAL.py
	@ echo "Calculating HAL scores..."
	@ python snv_run_HAL.py $< $@


../../processed_data/snv/snv_SPANR_scores.txt: ../../processed_data/snv/snv_data_clean.txt snv_run_SPANR.py
	@ echo "Calculating SPANR scores..."
	@ python snv_run_SPANR.py $< $@


../../processed_data/snv/snv_models_pr_curves_all.txt: ../../processed_data/snv/snv_func_annot.txt \
	../../processed_data/snv/snv_SPANR_scores_capped.txt \
	../../processed_data/snv/snv_HAL_scores.txt
	@ echo "Creating PR and ROC curves..."
	@ python snv_pr_curves.py


../../figs/snv/fig2C_snv_density.tiff: ../../processed_data/snv/snv_func_annot.txt fig2C_rel_position.R
	@ echo "Making Figure 2..."
	@ $(Rscript) fig2C_rel_position.R



