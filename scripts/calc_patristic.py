import dendropy
from dendropy import treecalc
import numpy as np

all_data = []

tree = dendropy.Tree.get_from_path("H1_noDup.tree", "newick")
pdm = treecalc.PatristicDistanceMatrix(tree)
for i, t1 in enumerate(tree.taxon_set):
    for t2 in tree.taxon_set[i+1:]:
        all_data.append(pdm(t1, t2))
        
print(np.mean(all_data))
