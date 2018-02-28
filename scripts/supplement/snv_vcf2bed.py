# convert ExAC VCF annotation file to BED format. Add ID column with unique
# identifier and strand. Format:
# 1. chrom
# 2. chromStart
# 3. chromEnd
# 4. name
# 5. score (leave blank)
# 6. strand

import gzip
import argparse

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='compressed VCF input file')
	parser.add_argument('outfile', help='Name of output file')
	args = parser.parse_args()

	with gzip.open(args.infile, 'rb') as infile:
		index = 0
		with gzip.open(args.outfile, 'wb') as outfile:
			for line in infile:
				index += 1
				name = 'SNV_' + str(index)
				line = infile.readline()
				fields = line.strip().split('\t')
				# VCF is 1 based
				chrom, position = fields[:2]
				chrom = 'chr' + chrom
				# VCF is 1-based, convert to 0-based
				start = int(position) - 1
				end = position

				info = fields[7]
				# grab everything after 'CSQ=' field identifier
				csq = info.split('CSQ=')[1]
				# CSQ split by '|' , strand given in field 20
				csq_fields = csq.split('|')
				strand = csq_fields[20]

				if strand == '':
					strand = '.'

				new_fields = map(str, [chrom, start, end, name, '.', strand])
				outfile.write('\t'.join(new_fields) + '\n')