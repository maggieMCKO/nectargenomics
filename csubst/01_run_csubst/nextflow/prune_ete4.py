
import sys
from ete4 import Tree
# import ete3 # same problem
from csubst import sequence

treefile = sys.argv[1]
align = sys.argv[2]
outtreefile = sys.argv[3]
print('treefile: ', treefile)
print('align: ', align)
print('outtreefile: ', outtreefile)

# test
# treefile='/home/mpg08/mko/Nectar/analysis/csubst/00_inputs/new_sp_tree.nw'
# align='/home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst_2024_proto/nf3_cusbst1.4_env/results_treeswift/rm_HLcinPul1_MBD5_rna-XM_015289876.2.fa'
# outtreefile = '/home/mpg08/mko/Nectar/analysis/csubst/01_run_csubst_2024_proto/nf3_cusbst1.4_env//pruned_tree_MBD5_rna-XM_015289876.2.nw'

# read tree
t = Tree(open(treefile))
# len(t)
# print(t)

# read align
seqs = sequence.read_fasta(path=align)
seq_names = list(seqs.keys())


# prune tree
t.prune(seq_names)
# len(t)

# export pruned tree
t.write( outtreefile )