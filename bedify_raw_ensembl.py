#/usr/bin/python

# bedify_raw_ensembl.py
# Caleb Matthew Radens
# 2016_1_14

### This script parses an output from ensembl biomart and writes to file 
###  a BEDFILE formatted version of the data:
###  chr | start | end | all columns from input file
###  
### start is 'Gene Start (bp)'
### end is 'Gene End (bp)'
###
###
### USAGE:
### > python bedify_raw_ensembl.py input_file_name.gz output_file_name.gz
###
### input_file is sourced from:
###  "/project/voight_subrate/cradens/noncoding_seq_context/data/not_mine"
### output_file is sourced from:
###  "/project/voight_subrate/cradens/noncoding_seq_context/data/results"

import sys
import os
import gzip
from helper_functions import remove_all

print "Initiating bedify_raw_ensembl.py"
print "Argument List:", str(sys.argv[1:])

if (len(sys.argv) != 3):
	raise Exception("Expected two command arguments.")
if (str(sys.argv[1])[-3:] != ".gz") and (str(sys.argv[2])[-3:] != ".gz"):
	raise Exception("Expected command arguments to be 2 file names ending in '.gz'")

print "Passed script checks."

# Where data are stored
data_folder = "/project/voight_subrate/cradens/noncoding_seq_context/data/"
noncoding_raw = data_folder+"not_mine/"+str(sys.argv[1])
OUT_FILE = data_folder+"results/"+str(sys.argv[2])

# 'rb' means read
f_IN = gzip.open(noncoding_raw, 'rb')
# 'wb' means write
f_OUT = gzip.open(OUT_FILE, 'wb')

# Some expected columns from the enesmbl bedfile
chr = 'Chromosome Name'
exon_chr_start = 'Exon Chr Start (bp)'
exon_chr_end = 'Exon Chr End (bp)'
gene_start = 'Gene Start (bp)'
gene_end = 'Gene End (bp)'
transcript_start = 'Transcript Start (bp)'
transcript_end = 'Transcript End (bp)'
genomic_coding_start = 'Genomic coding start'
genomic_coding_end = 'Genomic coding end'
cds_start = 'CDS Start'
cds_end = 'CDS End'

i = 0
for line in f_IN:

	# Remove newline chars and split by tab
	split_line = line.rstrip('\r\n').split('\t')

	# First line of file is header info
	if i == 0:

		# The first 3 columns of all bedfiles need to look like:
		head = ['chr','start','end']

		# array.extend() adds individual elements of split_line to
		#  head. Compare that to .append() which would add the entire
		#  split_line as a single element to the end of head.
		head.extend(split_line)

		print "Input Header: "
		print head
		head = '\t'.join(head)

		# Add updated header to new file
		print>>f_OUT, head

		# Column index of the chr name
		chr_i = split_line.index(chr)

		# Column indeces for all the coordinates
		exon_chr_start_i = split_line.index(exon_chr_start)
		exon_chr_end_i = split_line.index(exon_chr_end)
		gene_start_i = split_line.index(gene_start)
		gene_end_i = split_line.index(gene_end)
		transcript_start_i = split_line.index(transcript_start)
		transcript_end_i = split_line.index(transcript_end)
		genomic_coding_start_i = split_line.index(genomic_coding_start)
		genomic_coding_end_i = split_line.index(genomic_coding_end)
		cds_start_i = split_line.index(cds_start)
		cds_end_i = split_line.index(cds_end)

	# Else we're looking at the data below the header now
	else:
		# Check what is at each column. Convert to int, or, if empty, convert to 0
		all_indeces = [gene_start_i, gene_end_i, transcript_start_i, transcript_end_i, exon_chr_start_i, exon_chr_end_i,  
		 genomic_coding_start_i, genomic_coding_end_i, cds_start_i, cds_end_i]

		 # Initiate empty vector of genome coordinates
		all_coords = []
		for element in all_indeces:
			# Make sure the column isn't empty
			if len(split_line[element]) > 0:
				all_coords.append(int(split_line[element]))
			# Else make its value 0
			else:
				all_coords.append(0)

		start_coord = all_coords[0]
		end_coord = all_coords[1]
		all_coords = all_coords[2:]

		# Make sure the start and end of the gene aren't 0...
		if (start_coord != 0) and (end_coord != 0):

			remove_all(all_coords, 0)
			# Make sure the range of start and end is wider than all the other coordinates asscociated with this gene
			if abs(min(all_coords)-max(all_coords)) > abs(start_coord-end_coord):
				print start_coord
				print end_coord
				raise Exception("Range of start-end coords doesn't make sense at line: "+str(i)+"\n"+str(split_line))
		else:
			print start_coord
			print end_coord
			raise Exception("start_coord or end_coord was 0 at line: "+str(i)+"\n"+str(split_line))

		chr_name = split_line[chr_i]

		# Pre-pend chr, start, end, to the bedfile
		new_line = [chr_name, str(start_coord), str(end_coord)]
		new_line.extend(split_line)

		# Join array elements, separated by tab
		new_line = '\t'.join(new_line)

		# Print new line to the outfile
		print>>f_OUT, new_line
	i = i + 1
	# if i == 100:
	# 	break
f_IN.close()
f_OUT.close()









