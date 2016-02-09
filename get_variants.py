#!/usr/bin/python

### get_variants.py
### Caleb Matthew Radens
### 2015_2_5

### Get SNV info from a given region of the genome from a user-specified population.
###
### Given a directory, recursively search all folders for '.seq_context' files, look up what
### 	variants are within the sequence context for a specified population from 1000G, then write
### 	a separate file with the variant information in the same directory as the .seq_context file.
###
### Assumptions:
###	Population variant files are permenantly located at:
###		/project/voight_subrate/avarun/Research/mutation_rate/whole_genome/
### 	Population variant file names are formatted: 'all_var_POP_chr_loc' where POP is variable
###	Population variant files have 7 columns without headers:
###		chr | loci | ref_allele | variant_allele | ancesteral_allele | rsid | MAF
###
### Depends on:
###  /project/voight_subrate/cradens/noncoding_seq_context/script/sequence_functions.py
###
###  Usage:
###    python get_variants.py input_dir/ POP

import sys
import os
import pdb
import csv
home_base = "/project/voight_subrate/cradens/noncoding_seq_context/"
# Append my script folder to the file path so I may import my modules
sys.path.append(home_base+"script/")
from sequence_functions import get_variant_dict

print "Initiating get_variants.py"
print "Argument List:", str(sys.argv[1:])

if (len(sys.argv)-1 != 2):
	raise Exception("Expected two command arguments.")
in_dir = str(sys.argv[1])
pop = str(sys.argv[2])

if not (os.path.isdir(in_dir)):
	raise ValueError(in_dir+" not found. Is it a valid directory?")

pop_dict = get_variant_dict(Pop = pop)
# pop_dict[chrom]-> chrom_dict[loci] -> SNV_data

bed_files = list()
for root, subdirs, files in os.walk(in_dir):
	# ID files that I know have sequence context in them
	for f in files:
		if ".seq_context" in f:
			bed_files.append(os.path.join(root,f))


ACCEPTED_CHROMOSOMES = ["1","2","3","4","5","6","7","8","9","10",
			"11","12","13","14","15","16","17","18",
			"19","20","21","22","X"]

for bed_file in bed_files:
	with open(bed_file, 'rb') as in_file:
		i = 0
		snv_data = list()
		for line in in_file:
			# Skip header
			if i == 0:
				i = i + 1
				continue
			# Split line at ','
			line = line.split(",")
			chrom = line[0]
			start = line[1]
			end = line[2]
			# If the dictionary has key == chromosome of interest
			if pop_dict.has_key(chrom):
				chrom_dict = pop_dict[chrom]
				# Dictionary keys and values are all strings
				#  xrange(0:3) = [0,1,2]
				for loci in xrange(int(start), int(end)+1):
					if chrom_dict.has_key(str(loci)):
						line_to_save = list()
						line_to_save.append(chrom)
						line_to_save.append(str(loci))
						line_to_save.extend(chrom_dict[str(loci)])
						snv_data.append(line_to_save)
			elif chrom not in ACCEPTED_CHROMOSOMES:
				raise StandardError("Chromosome '"+chrom+"' not an accepted chromosome.")
			i = i + 1
		# If varaints were found, write them to file
		if len(snv_data) != 0:
			directory = os.path.dirname(bed_file)
			new_file = os.path.join(directory,
						os.path.basename(directory)+"."+pop+".SNV")
			with open(new_file, 'wb+') as out_file:
				writer = csv.writer(out_file)
				writer.writerows(snv_data)
			print str(len(snv_data))+" variants found in file:\n"+bed_file
		else:
			print "No variants found in file:\n"+bed_file
