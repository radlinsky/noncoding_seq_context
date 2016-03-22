#!/usr/bin/python

### gen_sims.py
### Caleb Matthew Radens
### 2015_2_15

###    This script runs Varun's sequence context simulator and saves results to file.
###
###     Arguments:
###            in_DIR/: where is the seq_context file located? 
###                extant, non-empty directory, that contains a single '.seq_context.BED' file
###            delim: what are the seq_context file deliminators? (default: ,)
###                string
###            Column_index: which column has the fasta sequence?
###                integer >= 0
###            Pop: What population should the simulator use?
###                string, "AFR", "ASN", or "EUR"
###            Line_s: which row should fasta sequence start being pulled from to simulate on?
###                integer, valid row index
###            Line_e: which row to end pulling fasta sequence from?
###                integer, valid row index
###        
###        Note:
###            Also writes to file the results+"_FINISHED" so you can check to make sure all the rows
###                you expected to generate simulations for actually got simulated.
###
###        Depends:
###            On Varun's find_simulated_variants module [see file path below]
###
###        Usage:
###            python gen_sims.py in_DIR/ delim Column_index Pop Line_s Line_e


#!/usr/bin/python
import sys, os, pdb
from subprocess import call
home_base = '/project/voight_subrate/cradens/noncoding_seq_context/script/'
sys.path.append(home_base)
from sequence_functions import *
varun_scripts = '/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks'
sys.path.append(varun_scripts)
from find_simulated_variants import *

in_DIR = str(sys.argv[1])
delim = str(sys.argv[2])
Column_index = int(sys.argv[3])
Pop = str(sys.argv[4])
Line_s= int(sys.argv[5])
Line_e= int(sys.argv[6])

fasta = get_seq_context(Directory = in_DIR, Delim = delim, Col = Column_index, Start = Line_s, End = Line_e)

n_lines = Line_e-Line_s

for Line in xrange(0,n_lines+1):
    seq = fasta[Line]
    new_file_path = os.path.join(in_DIR,Pop+'_1000_sim'+"_"+str(Line_s+Line))
    find_simulated_variants(fastaseq = seq, pop = Pop, filesave = new_file_path, nsim = 1000)
    finish_file = os.path.join(in_DIR,Pop+'_FINISHED_'+str(Line_s+Line))
    print "Finish file written to:"
    print finish_file
    with open(finish_file, 'wb'):
        pass
