import sys, os

sys.path.append("/gpfs/ysm/project/kleinstein/mw957/repos/spectral-tree-inference/spectraltree")

import pickle
import numpy as np
import utils
import generation
import reconstruct_tree
import dendropy
import scipy

from dendropy.model.discrete import simulate_discrete_chars, Jc69, Hky85
from dendropy.calculate.treecompare import symmetric_difference

tree_path = "/home/mw957/project/repos/spec_tree/data/skygrid_J2.newick"
fasta_path = "/home/mw957/project/repos/spec_tree/data/H3N2_NewYork.fasta"

H3N2_tree = dendropy.Tree.get(path=tree_path, schema="newick")
N = 1000
data_HKY = simulate_discrete_chars(N, H3N2_tree, Hky85(kappa = 2))
raxml_HKY = reconstruct_tree.RAxML()
raxml_HKY_tree = raxml_HKY(data_HKY, raxml_args="-T 2 --HKY85 -c 1")

RF,F1 = reconstruct_tree.compare_trees(raxml_HKY_tree, H3N2_tree)
print("")
print("NJ: ")
print("RF = ",RF)
print("F1% = ",F1)
print("")

with open("/gpfs/ysm/scratch60/morgan_levine/mw957/raxml_hky.pickle", "wb") as f:
    pickle.dump(raxml_HKY_tree, f)
