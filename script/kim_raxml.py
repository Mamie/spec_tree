import sys, os

sys.path.append("/gpfs/ysm/project/kleinstein/mw957/repos/spectral-tree-inference/spectraltree")

import numpy as np
import utils
import generation
import reconstruct_tree
import dendropy
import scipy
import time
from itertools import product
import matplotlib.pyplot as plt
import pandas as pd

from dendropy.model.discrete import simulate_discrete_chars, Jc69, Hky85
from dendropy.calculate.treecompare import symmetric_difference

def run_method(method, tree, seqs, threshold = None):
    ch_list = list()
    taxons = [x.label for x in tree.taxon_namespace]
    for t in taxons:
        ch_list.append([x.symbol for x in seqs[t]])
    ch_arr = np.array(ch_list)
    ch_arr[ch_arr == "U"] = "T"

    ch_dendro = dendropy.DnaCharacterMatrix()
    ch_dendro.taxon_namespace = tree.taxon_namespace
    for t, taxon in enumerate(taxons):
        ch_dendro.new_sequence(tree.taxon_namespace[t], ch_arr[t,:].tolist())
 
    if method == "RaXML":
        raxml_HKY = reconstruct_tree.RAxML()
        start_time = time.time()
        tree_rec = raxml_HKY(ch_dendro, raxml_args="-T 2 --HKY85 -c 1")      
    if method == "SNJ":
        snj = reconstruct_tree.SpectralNeighborJoining(reconstruct_tree.HKY_similarity_matrix_missing_data)
        start_time = time.time()
        tree_rec = snj(ch_arr, tree.taxon_namespace)
    if method == "NJ":
        nj = reconstruct_tree.NeighborJoining(reconstruct_tree.HKY_similarity_matrix_missing_data)
        start_time = time.time()
        tree_rec = nj(ch_arr, tree.taxon_namespace)
    if method == "STR + NJ":
        spectral_method = reconstruct_tree.SpectralTreeReconstruction(reconstruct_tree.NeighborJoining, reconstruct_tree.HKY_similarity_matrix_missing_data)
        start_time = time.time()
        tree_rec = spectral_method.deep_spectral_tree_reonstruction(ch_arr, reconstruct_tree.HKY_similarity_matrix_missing_data, 
                                                            taxon_namespace = tree.taxon_namespace, 
                                                            threshhold = threshold, min_split = 5)
    if method == "STR + SNJ":
        spectral_method = reconstruct_tree.SpectralTreeReconstruction(reconstruct_tree.SpectralNeighborJoining, reconstruct_tree.HKY_similarity_matrix_missing_data)
        start_time = time.time()
        tree_rec = spectral_method.deep_spectral_tree_reonstruction(ch_arr, reconstruct_tree.HKY_similarity_matrix_missing_data, 
                                                            taxon_namespace = tree.taxon_namespace, 
                                                            threshhold = threshold, min_split = 5)
    if method == "STR + RaXML":
        spectral_method = reconstruct_tree.SpectralTreeReconstruction(reconstruct_tree.RAxML,
                                                              reconstruct_tree.HKY_similarity_matrix_missing_data)
        start_time = time.time()
        tree_rec = spectral_method.deep_spectral_tree_reonstruction(ch_arr, reconstruct_tree.HKY_similarity_matrix_missing_data, 
                                                            taxon_namespace = tree.taxon_namespace, 
                                                            threshhold = threshold,
                                                            raxml_args = "-T 2 --HKY85 -c 1", min_split = 5)
    runtime = time.time() - start_time
    RF,F1 = reconstruct_tree.compare_trees(tree_rec, tree)
    print(method)
    if threshold is not None: print(threshold)
    print("--- %s seconds ---" % runtime)
    print("RF = ",RF)
    print("F1% = ",F1) 
    return([method, str(threshold), runtime, RF, F1])

tree_path = "/gpfs/ysm/project/kleinstein/mw957/repos/spec_tree/data/steps-1024_0.nex"
kim_tree = dendropy.Tree.get(path=tree_path, schema="nexus")
kim_seqs = dendropy.RnaCharacterMatrix.get(path=tree_path, schema="nexus")

methods = ["RaXML", "SNJ", "NJ", "STR + NJ", "STR + NJ", "STR + NJ", "STR + SNJ", "STR + SNJ", "STR + SNJ", "STR + RaXML", "STR + RaXML", "STR + RaXML"]
thresholds = [None, None, None, 32, 64, 128, 32, 64, 128, 32, 64, 128]

ms = []
ts = []
rts = []
rfs = []
f1s = []

for j in range(len(methods)):
    method = methods[j]
    threshold = thresholds[j]
    print(method, threshold)
    res = run_method(method, kim_tree, kim_seqs, threshold = threshold)
    ms.append(res[0])
    ts.append(res[1])
    rts.append(res[2])
    rfs.append(res[3])
    f1s.append(res[4])

perf_metrics = pd.DataFrame({'method': ms, 'threshold': ts, 'runtime': rts, 'RF': rfs, "F1": f1s})
perf_metrics.to_csv("/gpfs/ysm/project/kleinstein/mw957/repos/spec_tree/script/kim_1024.csv")

