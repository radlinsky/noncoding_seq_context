#!/usr/bin/python

### gen_variant_summaries.py
### Caleb Matthew Radens
### 2015_2_8

### Searches a directory recursively for folders that contain .POP.SNV files and
###     generates a summary table with the numbers of SNVs per population. Writes
###     to file a summary table with each row corresponding to a different folder.
###     For example, the folders may each correspond to a Gene Transcript.
###
### Arguments:
###     in_dir
###         extant directory   
###     out_file: where should the summary file be written to?
###
### Assumptions:
###     in_dir has sub-directories that contain .SNV files
###         These sub-dirs are named after a sequence of interest (e.g. 'ENSG...###')
###     SNV files are not zipped
###     There are no headers in the .SNV files
###
###  Usage:
###    python in_dir out_file

import os
import sys
import pdb
import csv
home_base = "/project/voight_subrate/cradens/noncoding_seq_context/script/"
# Add Caleb's script directory to PATH
sys.path.append(home_base)
from sequence_functions import n_variants

pdb.set_trace()

print "Initiating gen_variant_summaries.py"
print "Argument List:", str(sys.argv[1:])

if (len(sys.argv)-1 != 2):
    raise Exception("Expected two command arguments.")
in_DIR = str(sys.argv[1])
out_FILE = str(sys.argv[2])

if not (os.path.isdir(in_DIR)):
    raise ValueError("Input directory not found.")

# Find all '.SNV" files in their respective sub directories
sub_dirs_with_snv_files = dict()
for root, subdirs, files in os.walk(in_DIR):  
    for f in files:
        # If .SNV file in the directory
        if ".SNV" in f:
            # If the dictionary doesn't have the sub-dir defined as a key yet
            if not sub_dirs_with_snv_files.has_key():
                # Add the directory as a key -> list()
                sub_dirs_with_snv_files[root] = list()
            # Add the matched file to the list in the sub directory
            sub_dirs_with_snv_files[root].append(os.path.join(root,f))

possible_pops = ['ASW', 'TSI', 'PUR',
         'CHB', 'GBR', 'EUR',
         'CLM', 'AFR', 'ASN',
         'CHS', 'YRI', 'MXL',
         'JPT', 'LWK', 'FIN',
         'CEU', 'IBS']


summary_data = list()
header = list()
header.append("ensembleID")
header.extend(possible_pops)


empty_pop_data = [0]*len(possible_pops)

for directory, snv_files in sub_dirs_with_snv_files:
    summary_row = list()    
    dir_basename = os.path.basename(directory)
    # Add ensembleID name to summary row
    summary_row.append(dir_basename)
    
    # Add 0s for each population
    summary_row.extend(empty_pop_data)
    
    for snv_file in snv_files:
        snv_file_basename = os.path.basename(snv_file)
        snv_i = snv_file_basename.index(".SNV")
        pop = snv_file_basename[snv_i-3:snv_i]
        # +1 because the first column is the ensemblID
        summary_row_pop_col = possible_pops.index(pop)+1
        
        # Get the number of SNVs
        n_pop_snvs = n_variants(directory, pop)
        
        # Update 0 -> number of SNVs
        summary_row[summary_row_pop_col] = n_pop_snvs
    
    # Add the ensembleID data row to the full list
    summary_data.append(summary_row)
    
    
with open(out_FILE, 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(summary_data)
    
    
    
    
    
    
    
    
    
    