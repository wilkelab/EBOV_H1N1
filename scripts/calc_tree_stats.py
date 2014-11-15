## SJS. Created 11/15/14. Grab various stats from trees
from dendropy import *
import numpy as np
import os

head_dirs = ['EBOV_all_outbreaks/', 'EBOV_2014/', 'H1N1/', 'all_ebolavirus/']


outfile = '../phylogenies/tree_stats.txt'
outf = open(outfile, 'w')
outf.write('type\tfile\ttree_length\tmean_root-to-tip\tmean_patristic\n')


for dir in head_dirs:
    treefiles = os.listdir('../phylogenies/' + dir)
    type = dir.replace('/','')
    for file in treefiles:
        # only tree files
        if file.endswith('.tree'):
            name = file.split('.')[0]
            print name
            t = Tree.get_from_path('../phylogenies/' + dir + file, 'newick')
        
            # tree length
            tree_length = t.length()
        
            # root-to-tip
            treetips = t.leaf_nodes()
            rtt = []
            for tip in treetips:
                rtt.append( tip.distance_from_root() )
        
            # patristic
            pd = []
            dist = treecalc.PatristicDistanceMatrix(tree=t)
            for i, t1 in enumerate(t.taxon_set):
                for t2 in t.taxon_set[i:]:
                        d = dist(t1,t2)
                        pd.append( float(d) )
            
        
            outf.write(type + '\t' + name + '\t' + str(round(tree_length, 5)) + '\t' + str(round(np.mean(rtt), 5)) + '\t' + str(round(np.mean(pd), 5)) + '\n')
outf.close()
        
        