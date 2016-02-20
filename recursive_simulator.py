#!/usr/bin/python

### recursive_simulator.py
### Caleb Matthew Radens
### 2015_2_20

###    This script recursively searches a directory for sub-directories that contain 
###        *.seq_context.BED files and uses Varun's intergenic polymorphism rate
###        simulator to generate 1000 simulations of mutations from the fasta sequence
###        inside the *.seq_contex.BED file.
###
###    Arguments:
###        in_dir
###            extant directory with sub-directories that contain *.seq_context.BED files
###        POP: population to simulate
###            str
###            'EUR', 'AFR', or 'ASN'
###        delim: how are *.seq_context.BED files delimted?
###            str
###            'tab' or ','
###        skip: how many lines to skip before data begin in *.seq_context.BED files?
###            int >= 0
###        col: which column of *.seq_context.BED is fasta sequence in?
###            int >= 0
###
###    Assumptions:
###        The sub-dirs are named after a sequence of interest (e.g. 'ENSG...###')
###        The *.seq_context.BED files are all formated the same way
###        There is a single *.seq_context.BED per sub-directory
###        There is a single fasta sequence in each *.seq_context.BED file
###
###  Usage:
###    python recursive_simulator.py in_dir POP delim skip col

import sys
#import pdb

home_base = "/project/voight_subrate/cradens/noncoding_seq_context/script/"
# Add Caleb's script directory to PATH
sys.path.append(home_base)
from sequence_functions import *

# Import Varun's fasta sequence extraction module
varun_scripts = "/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks"
sys.path.append(varun_scripts)
from find_simulated_variants import *

pdb.set_trace()

print "Initiating recursive_simulator.py"
print "Argument List:", str(sys.argv[1:])

if (len(sys.argv)-1 != 5):
    raise Exception("Expected five command arguments.")
in_DIR = str(sys.argv[1])
POP = str(sys.argv[2])
delim = str(sys.argv[3])
skip = str(sys.argv[4])
col = sys.argv[5]

if not (os.path.isdir(in_DIR)):
    raise ValueError("Input directory not found.")

ACCEPTED_POPS = ["EUR", "AFR", "ASN"]
if POP not in ACCEPTED_POPS:
    raise ValueError("Population code: '"+POP+"' not an accepted population. See script desc.")

if delim != "tab" and delim != ",":
    raise ValueError("Delim needs to be 'tab' or ','")

if not skip.isdigit():
    raise ValueError("Skip needs to be an integer >= 0.")
skip = int(skip)
if skip < 0:
    raise ValueError("Skip needs to be an integer >= 0.")

if not col.isdigit():
    raise ValueError("col needs to be an integer >= 0.")
col = int(col)
if col < 0:
    raise ValueError("col needs to be an integer >= 0.")

files = get_seq_context_files(Directory = in_DIR)
coverage = 0

for f in files:
    parent_dir_path = os.path.dirname(f)
    
    fasta = get_seq_context(Directory = parent_dir_path, Delim = delim, Col = col)
    if len(fasta) == 0:
        raise ValueError("Fasta sequence non-existent in file:\n"+f)
    if fasta == "NA":
        raise ValueError("Fasta sequence was 'NA' in file:\n"+f)
    A = fasta.count("A")
    T = fasta.count("T")
    C = fasta.count("C")
    G = fasta.count("G")
    if A+T+C+G != len(fasta):
        raise ValueError("Fasta sequence had characters other than A, T, C, and G in file:\n"+f)
    
    coverage += len(fasta)
    
    
    parent_dir_base = os.path.basename(parent_dir_path)
    new_file_path = os.path.join(parent_dir_path, parent_dir_base+"_"+POP+"_1000_sim")
    try:
        find_simulated_variants(fastaseq = fasta, pop = POP, filesave = new_file_path, nsim = 1000)
    except:
        print "Error caught when trying to use 'find_simulated_variants' .. initiating pdb.set_trace()"
        pdb.set_trace()
        
print "Finished recursive_simulator.py"
print "Number of seq_context.BED files found: "+str(len(files))
print "Total coverage of all fasta sequences:"
print "Bp: "+str(coverage)
print "Kb: "+str(coverage/1000)
print "Mb: "+str(coverage/1000000)