#!/usr/bin/python

### sequence_functions.py
### Caleb Matthew Radens
### 2015_2_8

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
import pdb

EXPECTED_POPS = ['ASW', 'TSI', 'PUR',
 		'CHB', 'GBR', 'EUR',
 		'CLM', 'AFR', 'ASN',
		 'CHS', 'YRI', 'MXL',
		 'JPT', 'LWK', 'FIN',
		 'CEU', 'IBS']

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
	
	if not isinstance(Pop, str):
		raise ValueError("Pop needs to be a string.")
	
	print "Initiating get_variant_dict(Pop = "+Pop+")"

	if Pop not in EXPECTED_POPS:
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
	print "Population file to be analyzed:"+str(map(os.path.basename, pop_files))

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
			if not pop_dict.has_key(chrom):
				pop_dict[chrom]=dict({loci:variant_info})
			else:
				pop_dict[chrom][loci] = variant_info
	print "get_variant_dict(Pop = "+Pop+") complete."
	return pop_dict

def n_variants(Directory, Pop):
	"""
	Given a directory with SNV data from a given Population,
		count how many SNVs are present.
	
	Arguments:
		Directory: /extant_dir/
			Must contain SNV data for the population of interest
			
		Pop: POP of interest to get data for
			Three letter string (e.g. EUR for European)
	Returns: integer
	"""
	#pdb.set_trace()
	if not isinstance(Directory, str):
		raise ValueError("Directory needs to be a String")
	
	if not isinstance(Pop, str):
		raise ValueError("Pop needs to be a String")
	
	if not (os.path.isdir(Directory)):
		raise ValueError("Directory not found.")
	
	if Pop not in EXPECTED_POPS:
		raise ValueError(Pop+" not found in expected pops:\n"+str(EXPECTED_POPS))
	
	# Get list of all files in directory
	files = os.listdir(Directory)
	
	# Counter of SNVs
	n_snvs = 0
	
	# Counter of matched files
	n_files = 0
	for f in files:
		
		# If .POP.SNV in file name, add it to the list of matched files
		if "."+Pop+".SNV" in f:
			n_files = n_files+1
			if n_files > 1:
				raise ValueError("More than one file matched "+"."+Pop+".SNV in Directory.")
			
			# Count how many lines are in the file (corresponds to # of SNVs)
			# In retrospect, this is a really stupid way to do this. Oh well, it works.
			with open(os.path.join(Directory,f), 'rb') as file_handle:
				content = file_handle.readlines()
				for line in content:
					n_snvs = n_snvs + 1
	return n_snvs

def get_pop_snv_files(Directory, Pop):
	"""
	Searches Directory recursively for ".Pop.SNV." files and adds them to a list.
		
	Arguments:
		Directory: extant directory
			str
		Pop: POP of interest to match file names with
			Three letter str (e.g. 'EUR' for European)
	
	Returns:
		list (see above)
	"""
	if not isinstance(Directory, str):
		raise ValueError("Directory needs to be a string.")
	
	if not isinstance(Pop, str):
		raise ValueError("Pop needs to be a string.")
	
	if not os.path.isdir(Directory):
		raise ValueError(Directory +" is not a valid Directory.")
	
	if Pop not in EXPECTED_POPS:
		raise ValueError(Pop+" is not a valid population. Format: 'EUR'")
	
	matched_files = list()
	
	for root, subdirs, files in os.walk(Directory):
		for f in files:
			if Pop+".SNV" in f:
				matched_files.append(os.path.join(root,f))
				
	return matched_files

def get_variants(SNV_files, Delim = ",", Skip = 0):
	"""
	Given a file or list of files with SNV data, return list
		of variants.
	
	Arguments:
		SNV_files: list of files with SNV data
			list of strings
		Delim: how are SNV files delimited?
			str: tab 'tab' or comma ',' were tested
		Skip: how many lines should be skipped in each SNV file before data begins?
			integer >= 0
			
	Assumptions:
		Each line in SNV file is formatted the same way:
		
	Returns: list of SNV_variants
		E.g. [[chr,locus,ref. allele,variant allele,ancestral allele,rsid,maf],[etc...],etc...]
	"""
	if not is_sequence(SNV_files):
		raise ValueError("SNV files needs to be a list of files.")
	
	if not isinstance(Delim, str) or (Delim != "," and Delim != "tab"):
		raise ValueError("Delim needs to be a string: ',' or 'tab'.")
	
	if Delim == "tab":
		delim = "\t"
	else:
		delim = ","
	
	if not isinstance(Skip, int) or Skip < 0:
		raise ValueError("Skip needs to be an integer >= 0.")
	
	SNV_variants = list()
	
	for f in SNV_files:
		if not isinstance(f, str) or not os.path.isfile(f):
			raise ValueError("SNV files needs to be a list of file path strings.")
		with open (f, 'rb') as handle:
			for line in handle:
				if Delim not in line:
					raise ValueError("You fool! The variant files aren't all delimited by: '"+delim+"'!")
				variant_data_list = line.rstrip("\r\n").split(delim)
				SNV_variants.append(variant_data_list)
				
	return SNV_variants
	

def unique_variants(SNV_list):
	"""
		Given list of SNV_variants, return list removed of duplicates.
		
		Arguments:
			SNV_list: list of variants with the following format:
				[[chr,locus,ref. allele,variant allele,ancestral allele,rsid,maf],[etc...],etc...]
				
		Assumptions:
			SNV_list is from 1000 Genomes, so it is sufficient to uniquely call a variant by its [chr, start, end].
				(or by chr, start, end)
		
		Returns: (set of unique (via [chr, str, end]) variants, list of duplicate variants)
	"""
	if not is_sequence(SNV_list):
		raise ValueError("SNV_list needs to be a list.")
	
	if len(SNV_list) == 0:
		return SNV_list
	
	if len(SNV_list[0]) != 7:
		raise ValueError("SNV data aren't formatted as expected.")
	
	unique = dict()
	discard = list()
	
	for snv in SNV_list:
		chrom, loci = snv[0:2]
		if not isinstance(chr, str) or not isinstance(loci, str):
			raise ValueError("Chrom and loci need to be strings.")
		if not chrom.isdigit():
			raise ValueError("Chrom needs to be a digit in string format '#'.")
		if not loci.isdigit():
			raise ValueError("Loci needs to be a digit in string format '###'.")
		
		chr_loci = chrom+loci
		
		# If unique doesn't have chr_loci in it yet:
		if not unique.has_key(chr_loci):
			unique[chr_loci] =  snv
		# Else chr_loci already found, add the snv to the discard pile
		else:
			discard.append(snv)
	
	return unique, discard


def is_sequence(arg):
	"""
		Checks if arg is a list. If it is, return True, else return False.
	"""
	return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))
	
def get_seq_context_files(Directory):
	"""
	Searches Directory recursively for ".seq_context.BED" files and adds them to a list.
		
	Arguments:
		Directory: extant directory
			str
	
	Assumptions:
		There is a single *.seq_context.BED file in a given directory.
	
	Returns:
		list (see above)
	"""
	if not isinstance(Directory, str):
		raise ValueError("Directory needs to be a string.")
	
	if not os.path.isdir(Directory):
		raise ValueError(Directory +" is not a valid Directory.")
	
	matched_files = list()
	matched_subdirs = set()
	
	for root, subdirs, files in os.walk(Directory):
		for f in files:
			if ".seq_context.BED" in f:
				if root in matched_subdirs:
					raise BaseException("More than one .seq_context.BED found in directory:\n"+root)
				matched_files.append(os.path.join(root,f))
				matched_subdirs.add(root)
				
	return matched_files
	
def get_seq_context(Directory, Delim, Col):
	"""
		Given a directory with a ".seq_context.BED" file, return the sequence context
		
		Arguments:
			Directory: extant
			Delim: how is file delimited?
				str: 'tab' or ','
			Col: integer >= 0
				int
		
		Assumptions:
			Sequence is in given column, in fasta format
			There is only one line in the seq_context.BED file
			
		Returns:
			fasta sequence (str)
	"""
	if not isinstance(Directory, str):
		raise ValueError("Directory needs to be a directory path string")
	
	if not os.path.isdir(Directory):
		raise ValueError(Directory+" not found. Is it a full valid path?")
	
	if not isinstance(Delim, str) or (Delim != "," and Delim != "tab"):
		raise ValueError("Delim needs to be a string: ',' or 'tab'.")
	
	if Delim == "tab":
		delim = "\t"
	else:
		delim = ","
	
	if not isinstance(Col, int) or Col < 0:
		raise ValueError("Col needs to be an integer >= 0.")
	
	files = os.listdir(Directory)
	matched = list()
	
	for f in files:
		if ".seq_context.BED" in f:
			matched.append(os.path.join(Directory,f))
			
	if len(matched) == 0:
		raise BaseException("No 'seq_context.BED files found in directory: "+Directory)
	
	if len(matched) >1:
		raise BaseException("More than 1 'seq_context.BED files found in directory: "+Directory)
	
	with open(matched[0], 'rb') as handle:
		line_count = 0
		for line in handle:
			if delim not in line:
				raise ValueError("You fool! The seq context file wasn't delimited by: '"+delim+"'!")
			
			split_line = line.rstrip("\n\r").split(delim)
			
			fasta = split_line[Col]
			
			line_count += 1
			
	if line_count != 1:
		raise BaseException("There were "+str(line_count)+" lines in file. Expected a single line. File:\n"+matched[0])
	
	return fasta
		
	
	
	