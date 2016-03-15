#!/usr/bin/python
import sys, os, pdb
from subprocess import call
home_base = '/project/voight_subrate/cradens/noncoding_seq_context/script/'
sys.path.append(home_base)
from sequence_functions import *
varun_scripts = '/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks'
sys.path.append(varun_scripts)
from find_simulated_variants import *
#pdb.set_trace()
in_DIR = str(sys.argv[1])
delim = str(sys.argv[2])
Column_index = int(sys.argv[3])
Pop = str(sys.argv[4])
Line_s= int(sys.argv[5])
Line_e= int(sys.argv[6])
fasta = get_seq_context(Directory = in_DIR, Delim = delim, Col = Column_index, Start = Line_s, End = Line_e)

for Line in xrange(Line_s,Line_e+1):
    seq = fasta[Line]
    new_file_path = os.path.join(in_DIR,Pop+'_1000_sim'+"_"+str(Line))
    find_simulated_variants(fastaseq = seq, pop = Pop, filesave = new_file_path, nsim = 1000)
    finish_file = os.path.join(in_DIR,Pop+'_FINISHED_'+str(Line))
    print "Finish file written to:"
    print finish_file
    with open(finish_file, 'wb'):
        pass
