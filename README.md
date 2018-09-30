@author : Kimberly Insigne
kinsigne@ucla.edu

Feb 19, 2018

Please don't hesistate to email me with any questions.

Here's a basic description of where the data lies and how you can use it in your
own work. 

*To locate the raw read counts per bin:*
- There is an Excel file with all of the read counts per bin for each construct.
- The version of the SNV library used in the main paper is located at 
`processed_data/snv/snv_v2_all_alignments.csv`
- The data for the SRE library is located under `processed_data/sre/dhfr/dhfr_all_alignments.csv` and
`processed_data/sre/smn1/smn1_all_alignments.csv` for the two different intron backbones. 
- Raw sequencing data is not available at this time due to storage constraints but will
be made public upon publication.

*To located processed and cleaned data:*
- For the SNV library, the Excel file is processed by `process_scripts/snv/snv_data_clean.R`
and produces the file `processed_data/snv/snv_data_clean.txt`. **This is probably
the file you will want to work with**, it contains the calculated exon inclusion index
as well as information integrated from the reference (`ref/snv/snv_ref_formatted_converted_original_seq.txt`)
- `processed_data/snv/snv_simple_list.txt` contains a simplified list of the SDVs 
identified in the ExAC library and contains coordinates, the reference and alternate alleles,
and the delta exon inclusion index.
