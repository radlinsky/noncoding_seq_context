#!/usr/bin/python

### get_variants.py
### Caleb Matthew Radens
### 2015_2_2

### Get the SNVs in a given region of the genome from all populations.
###
### Given a directory with folders that contain .BED files that each have 'seq_context' column,
### 	look up what variants are in the sequence context, and write a new file with the variants.
###
###
###
###
###
###
###
###
###
###
###
###
### Assumptions:
###	Population variant files are permenantly located at:
###		/project/voight_subrate/avarun/Research/mutation_rate/whole_genome/
### 	Population variant file names are formatted: 'all_var_POP_chr_loc' where POP is variable
###	Population variant files have 7 columns without headers:
###		chr | loci | ref_allele | variant_allele | ancesteral_allele | rsid | MAF
###
###
###
###
###
###
###
###
###
###
###  Usage:
###    python get_variants.py input_dir/ pop

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

bed_files = list()
for root, subdirs, files in os.walk(in_dir):
	for f in files:
		if ".seq_context" in f:
			bed_files.append(os.path.join(root,f))

	for bed_file in bed_files:
		with open(bed_file, 'rb') as in_file:
			i = 0
			for line in in_file:
				# Skip header
				if i == 0:
					i = i + 1
					continue
				chrom = line[0]
				start = line[1]
				end = line[2]
				snv_data = list()
				if pop_dict.has_key(chrom):
					chrom_dict = pop_dict[chrom]
					for loci in xrange(int(start), int(end)+1):
						if chrom_dict.has_key(str(loci)):
							line_to_save = list()
							line_to_save.append(chrom)
							line_to_save.append(str(loci))
							line_to_save.extend(chrom_dict[str(loci)])
							snv_data.append(line_to_save)
				else:
					raise StandardError("Chromosome '"+chrom+"' not in pop_dict.")
				if len(snv_data) == 0:
					print "No variants found for: "+new_file 
				else:
					directory = os.path.dirname(bed_file)
					new_file = os.path.join(directory,
								os.path.basename(directory)+"."+pop+".SNV")
					with open(new_file, 'wb+') as out_file:
						writer = csv.writer(out_file)
						writer.writerows(row_group_file)
			i = i + 1
			
		



# python get_variants.py /project/voight_subrate/cradens/noncoding_seq_context/data/results/ensemble_grch37/miRNA/ENSG00000221524/ EUR






