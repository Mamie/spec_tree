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
from sklearn.decomposition import TruncatedSVD

def HKY_similarity_matrix(observations, classes=None, verbose = False):
    m, N = observations.shape
    if classes is None:
        classes = np.unique(observations)
    k = len(classes)
    # From Tamura, K., and M. Nei. 1993
    # for each pair of sequences, 
    # 1. estimate the average base frequency for pairs of sequences
    # 2. compute purine transition proportion P1 (A <-> G)
    # 3. compute pyrimidine transition proportion P2 (T <-> C)
    # 3. compute transversion proportion Q (A <-> C, A <-> T, G <-> C, G <-> T)

    if verbose: print("Computing the average base frequency for each pair of sequences...")
    g = {}
    for x in classes:
        obs_x = observations == x
        g[x] = np.array([np.mean(np.hstack([a, b])) for a, b in product(obs_x, repeat = 2)]).reshape((m, m))
    
    g["R"] = g["A"] + g["G"]
    g["Y"] = g["T"] + g["C"]
    
    # compute transition and transversion proportion
    if verbose: print("Computing transition and transversion proportion for each pair of sequences...")
    P = {}
    for i, x in enumerate(classes):
        other_classes = np.delete(classes, i)
        for y in other_classes:
            P_x_y = np.array([np.mean(np.logical_and(a == x, b == y)) for a, b in product(observations, repeat = 2)]).reshape((m, m))
            P[x + y] = P_x_y
            
    P_1 = P['AG'] + P["GA"]
    P_2 = P['CT'] + P['TC']
    Q = P['AC'] + P['CA'] + P['AT'] + P['TA'] +\
        P['GC'] + P['CG'] + P['GT'] + P['TG']

    # compute the similarity (formula 7)
    if verbose: print("Computing similarity matrix")
    R = (1 - g["R"]/(2 * g["A"] * g["G"]) * P_1 - 1 / (2 * g["R"]) * Q)
    Y = (1 - g["Y"]/(2 * g["T"] * g["C"]) * P_2 - 1 / (2 * g["Y"]) * Q)
    T = (1 - 1/(2 * g["R"] * g["Y"]) * Q)
    S = np.sign(R) * (np.abs(R))**(2 * g["A"] * g["G"] / g["R"])
    S += np.sign(Y) * (np.abs(Y))**(2 * g["T"] * g["C"] / g["Y"])
    S += np.sign(T) * (np.abs(T))**(2 * (g["R"] * g["Y"] - g["A"] * g["G"] * g["Y"] / g["R"] - g["T"] * g["C"] * g["R"] / g["Y"]))

    return S

def check_is_bipartition(tree, bool_partition):
    bipartitions = [str(x)[::-1] for x in tree.encode_bipartitions()]
    partition_1 = "".join(list(bool_partition.astype('int').astype('str')))
    partition_2 = "".join(list((1 - bool_partition).astype('int').astype('str')))
    is_bipartition = (partition_1 in bipartitions) or (partition_2 in bipartitions)
    return is_bipartition


SVD2_OBJ = TruncatedSVD(n_components=2, n_iter=7)
def svd2(mat):
    if (mat.shape[0] == 1) | (mat.shape[1] == 1):
        return 0
    elif (mat.shape[0] == 2) | (mat.shape[1] == 2):
        return np.linalg.svd(mat,False,False)[1]
    else:
        return SVD2_OBJ.fit(mat).singular_values_[1]

def partition_taxa(v,similarity,num_gaps):
    m = len(v)
    v_sort = np.sort(v)
    gaps = v_sort[1:m]-v_sort[0:m-1]
    ind_partition = np.argpartition(gaps, -num_gaps)[-num_gaps:]
    smin = 1000
    for p_idx in ind_partition:        
        threshold = (v_sort[p_idx]+v_sort[p_idx+1])/2
        if (p_idx>0) & (p_idx<m-2):
            bool_bipartition = v<threshold
            s_sliced = similarity[bool_bipartition,:]
            s_sliced = s_sliced[:,~bool_bipartition]
            s2 = svd2(s_sliced)
            if s2<smin:
                partition_min = bool_bipartition
                smin = s2
        elif p_idx == 0: partition_min = v <= v_sort[0]
        elif p_idx == m - 2: partition_min = v < v_sort[m-1]
    return partition_min

tree_path = "/home/mw957/project/repos/spec_tree/data/skygrid_J2.newick"
fasta_path = "/home/mw957/project/repos/spec_tree/data/H3N2_NewYork.fasta"

H3N2_tree = dendropy.Tree.get(path=tree_path, schema="newick")

B = 50
N = [50, 100, 400, 600, 800, 1000] 
mean_is_bipartition = []

for n in N:
    print(n)
    bipartitions = []
    start_time = time.time()
    for b in range(B):
        print(b)
        data_HKY = simulate_discrete_chars(n, H3N2_tree, Hky85(kappa = 2), mutation_rate=0.1)
        ch_list = list()
        for t in data_HKY.taxon_namespace: 
            ch_list.append([x.symbol for x in data_HKY[t]])
        ch_arr = np.array(ch_list)
        
        HKY_sim = HKY_similarity_matrix(ch_arr)
        _, eigvec = np.linalg.eigh(HKY_sim)
        partition = partition_taxa(eigvec[:,-2], HKY_sim, 1)
        is_bipartition = check_is_bipartition(H3N2_tree, partition)
        bipartitions.append(int(is_bipartition))
    runtime = time.time() - start_time
    print("--- %s seconds ---" % runtime)

    mean_is_bipartition.append(np.mean(bipartitions))

metrics = pd.DataFrame({"N": N, "valid partition freq": mean_is_bipartition})
metrics.to_csv("/gpfs/ysm/project/kleinstein/mw957/repos/spec_tree/script/split_metrics.csv")

