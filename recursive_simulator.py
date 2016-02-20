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
###    Note:
###        About the name, I was going to make this a recursive script, but I changed my
###            mind halfway through making it recursive. I'm too lazy to change the name. 
###
###  Usage:
###    python recursive_simulator.py in_dir POP delim skip col

import sys
from subprocess import Popen
#import pdb
import time

home_base = "/project/voight_subrate/cradens/noncoding_seq_context/script/"
# Add Caleb's script directory to PATH
sys.path.append(home_base)
from sequence_functions import *

# Import Varun's fasta sequence extraction module
varun_scripts = "/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks"
sys.path.append(varun_scripts)
from find_simulated_variants import *

#pdb.set_trace()

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

# Sub-routine
sub_routine = os.path.join(in_DIR,POP+"_never_gonna_give_u_up.py")
if os.path.isfile(sub_routine):
    print "FYI, we're overwriting something called "+sub_routine
with open(sub_routine, 'wb') as handle:
    handle.write("#!/usr/bin/python\n")
    handle.write("import sys, os, pdb\nfrom subprocess import call\n")
    handle.write("home_base = '/project/voight_subrate/cradens/noncoding_seq_context/script/'\n")
    handle.write("sys.path.append(home_base)\n")
    handle.write("from sequence_functions import *\n")
    handle.write("varun_scripts = '/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks'\n")
    handle.write("sys.path.append(varun_scripts)\n")
    handle.write("from find_simulated_variants import *\n")
    handle.write("in_DIR = str(sys.argv[1])\n")
    handle.write("delim = str(sys.argv[2])\n")
    handle.write("if delim == 'tab':\n\tdelim='\t'\n")
    handle.write("Column_index = int(sys.argv[3])\n")
    handle.write("Pop = str(sys.argv[4])\n")
    handle.write("fasta = get_seq_context(Directory = in_DIR, Delim = delim, Col = Column_index)\n")
    handle.write("if len(fasta) == 0:\n")
    handle.write("\traise ValueError('Fasta sequence not found in directory: '+in_DIR)\n")
    handle.write("if fasta == 'NA':\n")
    handle.write("\traise ValueError('Fasta sequence was NA in directory: '+in_DIR)\n")
    handle.write("A = fasta.count('A')\n")
    handle.write("T = fasta.count('T')\n")
    handle.write("C = fasta.count('C')\n")
    handle.write("G = fasta.count('G')\n")
    handle.write("if A+T+C+G != len(fasta):\n")
    handle.write("\traise ValueError('Fasta sequence had characters other than A, T, C, and G in dir: '+in_DIR)\n")
    handle.write("parent_dir_base = os.path.basename(in_DIR)\n")
    handle.write("new_file_path = os.path.join(in_DIR, parent_dir_base+'_'+Pop+'_1000_sim')\n")
    handle.write("find_simulated_variants(fastaseq = fasta, pop = Pop, filesave = new_file_path, nsim = 1000)\n")
    handle.write("parent_sqrd = os.path.dirname(in_DIR)\n")
    handle.write("finish_file = os.path.join(parent_sqrd,parent_dir_base+'_'+Pop+'_never_gonna_give_u_up_FINISHED')\n")
    handle.write("with open(finish_file, 'wb'):\n")
    handle.write("\tpass\n")
    

files = get_seq_context_files(Directory = in_DIR)
n_files = len(files)

for f in files:
    parent_dir_path = os.path.dirname(f)
    command = "bsub -e "+sub_routine[:-3]+".err "
    command += "-o "+sub_routine[:-3]+".out -q voight_normal "
    command += "python "+sub_routine+" "+parent_dir_path+" "+delim+" "+str(col)+" "+POP
    proc = Popen([command],shell=True, stdin=None, stdout=None, stderr=None, close_fds=True)

#pdb.set_trace()

all_done = False
# Wait until all the simulations have been run
while not all_done:
    finish_files = list()
    current_files = os.listdir(in_DIR)
    for f in current_files:
        if POP+"_never_gonna_give_u_up_FINISHED" in f:
            finish_files.append(os.path.join(in_DIR, f))
    if len(finish_files) == n_files:
        all_done = True
    else:
        time.sleep(5)

#pdb.set_trace()

for f in finish_files:
    os.remove(f)

coverage = 0
print "=============="
print "error file looks like:"
with open(sub_routine[:-3]+".err") as handle:
    for line in handle:
        print line.rstrip("\n\r")
        
#pdb.set_trace()

print "=============="
print "out file looks like:"
with open(sub_routine[:-3]+".out") as handle:
    for line in handle:
        print line.rstrip("\n\r")

#pdb.set_trace()

# Remove sub-routine and the error / out files
os.remove(sub_routine)
os.remove(sub_routine[:-3]+".err")
os.remove(sub_routine[:-3]+".out")
        
print "Finished recursive_simulator.py"
print "Number of seq_context.BED files found: "+str(len(files))