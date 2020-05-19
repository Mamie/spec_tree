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

def run_method(method, tree, threshold = None):
    data_HKY = simulate_discrete_chars(1000, tree, Hky85(kappa = 2), mutation_rate=0.1)
    ch_list = list()
    for t in data_HKY.taxon_namespace:
        ch_list.append([x.symbol for x in data_HKY[t]])
    ch_arr = np.array(ch_list)
    
    if method == "raxml":
        raxml_HKY = reconstruct_tree.RAxML()
        start_time = time.time()
        tree_rec = raxml_HKY(data_HKY, raxml_args="-T 2 --HKY85 -c 1")      
    if method == "snj":
        snj = reconstruct_tree.SpectralNeighborJoining(reconstruct_tree.HKY_similarity_matrix)
        start_time = time.time()
        tree_rec = snj(ch_arr, tree.taxon_namespace)
    if method == "nj":
        nj = reconstruct_tree.NeighborJoining(reconstruct_tree.HKY_similarity_matrix)
        start_time = time.time()
        tree_rec = nj(ch_arr, tree.taxon_namespace)
    if method == "nj sp deep":
        spectral_method = reconstruct_tree.SpectralTreeReconstruction(reconstruct_tree.NeighborJoining, reconstruct_tree.HKY_similarity_matrix)
        start_time = time.time()
        tree_rec = spectral_method.deep_spectral_tree_reonstruction(ch_arr, reconstruct_tree.HKY_similarity_matrix, 
                                                            taxon_namespace = tree.taxon_namespace, 
                                                            threshhold = threshold, min_split = 5)
    if method == "raxml sp deep":
        spectral_method = reconstruct_tree.SpectralTreeReconstruction(reconstruct_tree.RAxML,
                                                              reconstruct_tree.HKY_similarity_matrix)
        start_time = time.time()
        tree_rec = spectral_method.deep_spectral_tree_reonstruction(ch_arr, reconstruct_tree.HKY_similarity_matrix, 
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

tree_path = "/home/mw957/project/repos/spec_tree/data/H1N1_NY_Skygrid_cutoff6.newick"

H1N1_tree = dendropy.Tree.get(path=tree_path, schema="newick")
n_runs = 10

methods = ["raxml", "snj", "nj", "nj sp deep", "nj sp deep", "nj sp deep", "nj sp deep", "raxml sp deep", "raxml sp deep"]
thresholds = [None, None, None, 35, 64, 128, 256, 64, 128]

ms = []
ts = []
rts = []
rfs = []
f1s = []

for i in range(n_runs):
    for j in range(len(methods)):
        method = methods[j]
        threshold = thresholds[j]
        print(method, threshold)
        res = run_method(method, H1N1_tree, threshold = threshold)
        ms.append(res[0])
        ts.append(res[1])
        rts.append(res[2])
        rfs.append(res[3])
        f1s.append(res[4])

perf_metrics = pd.DataFrame({'method': ms, 'threshold': ts, 'runtime': rts, 'RF': rfs, "F1": f1s})
perf_metrics.to_csv("/gpfs/ysm/project/kleinstein/mw957/repos/spec_tree/script/h1n1_angle.csv")

