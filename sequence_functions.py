#!/usr/bin/python

### sequence_functions.py
### Caleb Matthew Radens
### 2015_2_3

### A bunch of helpful sequence-related functions I put together.
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
###

import sys
import os
# Useful for manipulating a list all at once
from operator import methodcaller, itemgetter


def get_variant_dict(Pop):
	"""Given a population [ex: 'EUR'], return a dictionary of dictionaries:
		dict{ chrom -> dict{ loci -> [ref_allele, variant, ancestral_allele, rsid, MAF] }}

		Assumptions:

		Population variant files are permenantly located at:
			/project/voight_subrate/avarun/Research/mutation_rate/whole_genome/
	 	Population variant file format: 'all_var_POP_chr_loc' where POP is variable
		Population variant files have 7 columns without headers:
			chr | loci | ref_allele | variant_allele | ancesteral_allele | rsid | MAF
	"""

	expected_pops = ['ASW', 'TSI', 'PUR',
	 		'CHB', 'GBR', 'EUR',
	 		'CLM', 'AFR', 'ASN',
			 'CHS', 'YRI', 'MXL',
			 'JPT', 'LWK', 'FIN',
			 'CEU', 'IBS']

	if not isinstance(Pop, str):
		raise ValueError("Pop needs to be a string.")

	if Pop not in expected_pops:
		raise ValueError("Expected valid POP name, but received: "+Pop)

	pop_dir = "/project/voight_subrate/avarun/Research/mutation_rate/whole_genome/"
	pop_files = list()
	for root, subdirs, files in os.walk(pop_dir):
		for f in files:
			if len(f) == 19 and "all_var" in f and "chr_loc" in f and Pop in f:
				pop_files.append(os.path.join(root,f))

	if len(pop_files) > 1:
		raise StandardError("More than 1 Pop file found in directory: "+str(pop_files))

	if len(pop_files) == 0:
		raise StandardError("Population not found in directory for some reason...")

	# map(function, list_or_things_to_apply_function_to)
	print "Population file to be anayzed:"+str(map(os.path.basename, pop_files))

	pop_dict = dict()
	with open(pop_files[0], 'rb') as file_handle:
		for line in file_handle:
			# Split line at the tab
			line = line.split()
			chrom = line[0]
			loci = line[1]
			variant_info = line[2:]
			# dictionary:
			# chr -> loci -> variant_info
			if len(pop_dict[chrom]) == 0:
				pop_dict[chrom]=dict({loci:variant_info})
			else:
				pop_dict[chrom][loci] = variant_info
