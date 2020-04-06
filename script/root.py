import sys, os

sys.path.append("/gpfs/ysm/project/kleinstein/mw957/repos/spectral-tree-inference/spectraltree")

import copy
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

tree_path = "/home/mw957/project/repos/spec_tree/data/skygrid_J2.newick"
fasta_path = "/home/mw957/project/repos/spec_tree/data/H3N2_NewYork.fasta"

H3N2_tree = dendropy.Tree.get(path=tree_path, schema="newick")
all_bipartitions = np.array([str(x)[::-1] for x in H3N2_tree.encode_bipartitions()][0:-1])

taxon_namespace_label = np.array([x.label for x in H3N2_tree.taxon_namespace])

def to_bool(partition_str):
    return np.array(list(partition_str)) == '1'

def min_partition_size(bipartition_encoding):
    n_ones = np.sum(np.array(list(bipartition_encoding)) == '1')
    n_zeros = np.sum(np.array(list(bipartition_encoding)) == '0')
    return(min(n_ones, n_zeros))

min_bipar = np.array([min_partition_size(x) for x in all_bipartitions])
filtered_bipar = all_bipartitions[np.where(min_bipar > 50)[0]]


N = [50, 100, 400, 600, 800, 1000]

Ns = []
par1s = []
par2s = []
RFs = []
F1s = []
rts = []

for n in N:
    print(n)
    data_HKY = simulate_discrete_chars(n, H3N2_tree, Hky85(kappa = 2), mutation_rate=0.1)
    ch_list = list()
    for t in data_HKY.taxon_namespace: 
        ch_list.append([x.symbol for x in data_HKY[t]])
    ch_arr = np.array(ch_list)
    HKY_sim = HKY_similarity_matrix(ch_arr)
    
    for partition in filtered_bipar:
        partition = to_bool(partition)
        par1_size = np.sum(partition)
        par2_size = np.sum(np.logical_not(partition))
        print("Partition size: ", par1_size, " vs ", par2_size)
        left_namespace = list(taxon_namespace_label[np.where(partition)[0]])
        left_taxa = dendropy.TaxonNamespace([taxon for taxon in H3N2_tree.taxon_namespace
            if taxon.label in left_namespace])

        T_left = copy.deepcopy(H3N2_tree).extract_tree_with_taxa_labels(labels = left_namespace)
        T_left.purge_taxon_namespace()
        s = T_left.as_string(schema = "newick")
        T_left = dendropy.Tree.get(data=s, schema="newick", taxon_namespace = left_taxa)
        right_namespace = list(taxon_namespace_label[np.where(np.logical_not(partition))[0]])
        right_taxa = dendropy.TaxonNamespace([taxon for taxon in H3N2_tree.taxon_namespace
            if taxon.label in right_namespace])
        T_right = copy.deepcopy(H3N2_tree).extract_tree_with_taxa_labels(labels = right_namespace)
        T_right.purge_taxon_namespace()
        s = T_right.as_string(schema = "newick")
        T_right = dendropy.Tree.get(data=s,
        schema="newick", taxon_namespace = right_taxa)
        
        start_time = time.time()
        joined_tree = reconstruct_tree.join_trees_with_spectral_root_finding(
            HKY_sim, T_left, T_right, taxon_namespace = H3N2_tree.taxon_namespace)
        runtime = time.time() - start_time
        
        RF,F1 = reconstruct_tree.compare_trees(joined_tree, H3N2_tree)
        
        Ns.append(n)
        par1s.append(par1_size)
        par2s.append(par2_size)
        RFs.append(RF)
        F1s.append(F1)
        rts.append(runtime)
        
perf_metrics = pd.DataFrame({'seqlength': Ns, 'par1_size': par1s, 'par2_size': par2s, 
                             'RF': RFs, "F1": F1s, "runtime": rts})
perf_metrics.to_csv("/gpfs/ysm/project/kleinstein/mw957/repos/spec_tree/script/rooting_metrics_normalized_2.csv")
