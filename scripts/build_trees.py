## SJS. Created 11/15/14. Build nucleotide trees w/ FastTree. 

import subprocess
import os

head_dirs = ['EBOV_2014/', 'H1N1/', 'all_ebolavirus/']
for dir in head_dirs:
    alnfiles = os.listdir('../alignments/' + dir)
    for file in alnfiles:
        alnfile = '../alignments/' + dir + file
        # only process nucleotide .aln files
        if alnfile.endswith('.aln') and 'aa' not in alnfile:
            treefile = '../phylogenies/' + dir + file.split('.')[0] + '.tree'
            print alnfile
            call_fasttree = 'FastTree -quiet -nt -gtr -nosupport ' + alnfile  + ' > ' + treefile
            run = subprocess.call(call_fasttree, shell=True)
            if run != 0:
                print "tree wouldn't build for " + alnfile
                print