#!/usr/bin/python

### add_seq_context_col.py
### Caleb Matthew Radens
### 2015_2_2

### This script appends a column of hg19 reference genome sequence to each file in a folder of files
### 	that each has chr, start, and end positions as columns. Discards all other columns. Output
###	is saved in a sub-directory with '_seq_context' appended to the root file name.
###
###  Arguments:
###	input_directory/: where are the files located? 
###	 extant, non-empty directory
###	delim: what are the file deliminators? (default: ,)
###	 string
###	chrom: which column is the chromosome column?
###	 integer, valid column index
###	start: which column is the start position column?
###	 integer, valid column index
###	end: which column is the end position column?
###	 integer, valid column index
###	before_after: nucleotides before or after the start/end (needed for padding)
###	 integer, >=0
###
###  Assumptions:
###    All files are structured and deliminated identically
###    There are single line headers for each file
###  
###  Notes:
###    Any row with a non-acceptable chromosome (1:23 and X) will be skipped: NA instead of sequence
###
###  Depends:
###    On Varun's get_seq_context_interval module [see file path below]
###
###  Usage:
###    python folderize_by_gene.py input_dir/ delim chrom start end before_after

import sys
import os
# Useful for manipulating a list all at once
from operator import methodcaller, itemgetter
varun_scripts = "/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks"
sys.path.append(varun_scripts)
from find_context import get_seq_context_interval
caleb_scripts = "/project/voight_subrate/cradens/noncoding_seq_context/script/"
sys.path.append(caleb_scripts)
from helper_functions import remove_all
import csv

print "Initiating add_seq_context_col.py"
print "Argument List:", str(sys.argv[1:])

if (len(sys.argv)-1 != 6):
	raise Exception("Expected five command arguments.")
in_dir = str(sys.argv[1])
delim = str(sys.argv[2])
chrom_i = int(sys.argv[3])
start_i = int(sys.argv[4])
end_i = int(sys.argv[5])
padding = int(sys.argv[6])

if not (os.path.isdir(in_dir)):
	raise ValueError(in_dir+" not found. Is it a valid directory?")

out_dir = "/"+os.path.basename(in_dir)+"_seq_context/"
out_dir = os.path.join(in_dir, out_dir)

# Extract names of all files in in_dir
in_files = os.listdir(in_dir)

if len(in_files) <= 0:
	raise Exception("No files found in directory.")

if len(delim) <= 0:
	raise ValueError("delim must be of length > 0")

if chrom_i < 0 or start_i < 0 or end_i < 0:
	raise ValueError("Column indeces need to be an integers >= 0.")

if padding < 0:
	raise ValueError("Padding needs to be an integer >0")

ACCEPTED_CHROMOSOMES = ["1","2","3","4","5","6","7","8","9","10",
			"11","12","13","14","15","16","17","18",
			"19","20","21","22","X"]

for file_name in in_files:
	# Append base path firectory to file_name
	full_file_name = os.path.join(in_dir, file_name)
	# Open file for reading in binary format
	with open(full_file_name, 'rb') as file_handle:
		# Get all lines, removed of \n, from the file, as a list
		all_lines = file_handle.read().splitlines()
		# Use operator.itemgetter() to slit each line in the list at delim
		all_lines = map(methodcaller("split", delim), all_lines)
	# Pop the header from the list of lines
	header = all_lines.pop(0)
	header.append("sequence_context_interval")
	if len(all_lines) > 0:
		# Use operator.itemgetter() to get all index-specified elements from list of lines
		chroms = map(itemgetter(chrom_i), all_lines)
		starts = map(itemgetter(start_i), all_lines)
		ends = map(itemgetter(end_i), all_lines)
		sequences = list()
		# Add the sequence conext interval for each chr, start, end in file
		for chrom, start, end in zip(chroms, starts, ends):
			if chrom not in ACCEPTED_CHROMOSOMES:
				print "'"+chrom+"' isn't a valid chrom in file: "+file_name
				sequences.append("NA")
			else:
				sequences.append(get_seq_context_interval(chrom,start,end,padding))

		# Check if 'NA' is the only type of sequence data in the file (if so, skip it)
		sequences_list = list(sequences)
		remove_all(sequences_list,"NA")
		if len(sequences_list) == 0:
			print "Skipping: "+file_name+" because no valid chromosomes were in it."
			continue

		columnized = zip(chroms, starts, ends, sequences)

		# If sub-directory doesn't exist yet, make it:
		out_file_path = os.path.join(out_dir, file_name)
		print out_file_path
		if not os.path.exists(os.path.dirname(out_file_path)):
		    try:
			os.makedirs(os.path.dirname(out_file_path))
		    except OSError as exc: # Guard against race condition
			if exc.errno != errno.EEXIST:
			    raise

		with open(out_file_path,"wb+") as out:
		    csv_out=csv.writer(out)
		    csv_out.writerow(header)
		    for row in columnized:
			csv_out.writerow(row)
			
	else:
		print "File: "+file_name+" had a header, but was empty."

	



























