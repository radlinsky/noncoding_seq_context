#!/usr/bin/python

### add_seq_context_col.py
### Caleb Matthew Radens
### 2015_2_15

###	This script appends a column of hg19 reference genome sequence to all .BED files in a folder
###		that have chr, start, and end positions as columns. Discards all other columns. Output
###		is saved as a new file with '.seq_context' appended to the root file name.
###
### 	Arguments:
###			input_directory/: where are the files located? 
###				extant, non-empty directory
###			delim: what are the file deliminators? (default: ,)
###				string
###			skip: how many lines to skip before looking for genomic loci info?
###				integer >= 0
###			chrom: which column is the chromosome column?
###				integer, valid column index
###			start: which column is the start position column?
###				integer, valid column index
###			end: which column is the end position column?
###				integer, valid column index
###			before_after: nucleotides before or after the start/end (needed for padding)
###				integer, >=0
###
###		Assumptions:
###			All files are structured and deliminated identically
###			All files end in '.BED'
###			chromosome format: (# = 1:22) if # == 'X' or 'Y', line is ignored.
###				'chr#' or '#'
###  
###		Notes:
###			Any row with a non-acceptable chromosome (only 1:23, X/Y rejected) will be skipped
###			Only writes to file if there is an acceptable chromosome
###			Recursively gets all files in directory that contain '.BED' but lack '.seq_context'
###
###		Depends:
###			On Varun's get_seq_context_interval module [see file path below]
###
###		Usage:
###			python add_seq_context_col.py input_dir/ delim skip chrom start end before_after

import sys
import os
varun_scripts = "/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks"
sys.path.append(varun_scripts)
from find_context import get_seq_context_interval
caleb_scripts = "/project/voight_subrate/cradens/noncoding_seq_context/script/generally_useful"
sys.path.append(caleb_scripts)
import csv
import pdb

print "Initiating add_seq_context_col.py"
print "Argument List:", str(sys.argv[1:])

if (len(sys.argv)-1 != 7):
	raise Exception("Expected seven command arguments.")
in_dir = str(sys.argv[1])
delim = str(sys.argv[2])
skip = int(sys.argv[3])
chrom_i = int(sys.argv[4])
start_i = int(sys.argv[5])
end_i = int(sys.argv[6])
padding = int(sys.argv[7])

if not (os.path.isdir(in_dir)):
	raise ValueError(in_dir+" not found. Is it a valid directory?")

# Extract names of all files in in_dir
in_files = os.listdir(in_dir)

if len(in_files) <= 0:
	raise Exception("No files or folders found in directory.")

if len(delim) <= 0:
	raise ValueError("delim must be of length > 0.")

# If tab-delimited, need to make sure it will be python-interpretable:
if delim != "," and delim != "tab":
	raise ValueError("This script was only tested with ',' or 'tab', not '"+delim+"'")

# Convert delim to '\t' if it is tab
if delim == "tab":
	delim_check = '\t'
# Else delim is kept the same ',' in this case
else:
	delim_check = delim

if skip < 0:
	raise ValueError("skip needs to be integer >= 0.")

# from bash, >>> python \t sends a 't' to python. >>> python $'\t' sends \t.
#    this little check assumes user meant \t
#if "t" in delim and len(delim)==1:
#	delim = "\t"

if chrom_i < 0 or start_i < 0 or end_i < 0:
	raise ValueError("Column indeces need to be an integers >= 0.")

if padding < 0:
	raise ValueError("Padding needs to be an integer >0")

ACCEPTED_CHROMOSOMES = list()
ACCEPTED_CHROMOSOMES.extend(["1","2","3","4","5","6","7","8","9","10",
			"11","12","13","14","15","16","17","18",
			"19","20","21","22"])

in_files = list()
for root, subdirs, files in os.walk(in_dir):
	for f in files:
		if ".BED" in f and ".seq_context" not in f:
			in_files.append(os.path.join(root,f))
			
all_seq_context_files = list()
skipped_files = list()
for full_file_name in in_files:
	new_lines = list()
	# Get base of full_file_name
	file_name = os.path.basename(full_file_name)
	# Open file for reading in binary format
	with open(full_file_name, 'rb') as file_handle:
		# Skip lines
		i = 0
		for line in file_handle:
			if i < skip:
				i+=1
				continue
			
			if delim_check not in line:
				raise ValueError("Delim '"+delim_check+"' not found in line # "+str(i))
			# Add line to all_lines
			split_line = line.rstrip('\r\n').split(delim_check)
			chrom = split_line[chrom_i]
			
			# Extract # from chromosome. Expected format: 'chr#' or '#'
			chrom = [str(s) for s in chrom.split("chr") if s.isdigit()]
			if len(chrom) > 1:
				raise ValueError("Unexpected chromosome format: "+split_line[chrom_i])
			elif len(chrom) == 0:
				chrom = split_line[chrom_i]
			elif len(chrom) == 1:
				chrom = chrom[0]
			if chrom not in ACCEPTED_CHROMOSOMES:
				print "'"+chrom+"' isn't an accepted chromosome at line # "+str(i)+" in file: "+file_name
				i+=1
				continue
			start = int(split_line[start_i])
			end = int(split_line[end_i])
			if start >= end:
				#pdb.set_trace()
				raise ValueError("Start is >= End at line:\n"+line)
			sequence = get_seq_context_interval(chrom,start,end,padding)
			split_line.append(sequence)
			new_lines.append(split_line)
			i+=1
	if len(new_lines) == 0:
		print "Skipping file "+file_name+" because it didn't have any accepted chromosomes."
		skipped_files.append(full_file_name)
		continue
	
	# Add '.seq_context' before '.BED' and write to file
	new_file_name = full_file_name[:-4]+".seq_context"+full_file_name[-4:]
	with open(new_file_name,"wb+") as out:
		csv_out=csv.writer(out)
		csv_out.writerows(new_lines)
	all_seq_context_files.append(new_file_name)
	
print "=============="
print "Retrieved sequence context from "+str(len(all_seq_context_files))+" .BED files."
print "New files written:\n"+str(all_seq_context_files)

print "=============="
print "Skipped "+ str(len(skipped_files)) + " files because they didn't have acceptable chromosome names."
print "Skipped files were:\n"+str(skipped_files)