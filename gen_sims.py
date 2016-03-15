#!/usr/bin/python
import sys, os, pdb
from subprocess import call
home_base = '/project/voight_subrate/cradens/noncoding_seq_context/script/'
sys.path.append(home_base)
from sequence_functions import *
varun_scripts = '/project/voight_subrate/avarun/Research/mutation_rate/scripts_for_folks'
sys.path.append(varun_scripts)
from find_simulated_variants import *
pdb.set_trace()
in_DIR = str(sys.argv[1])
delim = str(sys.argv[2])
Column_index = int(sys.argv[3])
Pop = str(sys.argv[4])
Line= int(sys.argv[5])
fasta = get_seq_context(Directory = in_DIR, Delim = delim, Col = Column_index, Lines = Line)
if len(fasta) == 0:
    raise ValueError('Fasta sequence not found in directory: '+in_DIR)
if fasta == 'NA':
    raise ValueError('Fasta sequence was NA in directory: '+in_DIR)
A = fasta.count('A')
T = fasta.count('T')
C = fasta.count('C')
G = fasta.count('G')
if A+T+C+G != len(fasta):
    raise ValueError('Fasta sequence had characters other than A, T, C, and G in dir: '+in_DIR)
new_file_path = os.path.join(in_DIR,Pop+'_1000_sim'+"_"+str(Line))
find_simulated_variants(fastaseq = fasta, pop = Pop, filesave = new_file_path, nsim = 1000)
finish_file = os.path.join(in_DIR,Pop+'_FINISHED_'+str(Line))
with open(finish_file, 'wb'):
    pass
