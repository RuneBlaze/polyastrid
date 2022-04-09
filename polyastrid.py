from copy import copy
from math import fsum
from joblib import Parallel, delayed
from random import random
import treeswift as ts
import numpy as np
from functools import lru_cache
import json
import argparse
import asterid as ad
import asterid as astrid
import time
from polybase import *
import os

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

with open(os.path.join(__location__, "probs.json")) as fh: correction_map_ = json.load(fh)

correction_map = {}
for k in correction_map_:
    correction_map[int(k)] = correction_map_[k]
for i in range(3):
    correction_map[i] = 0

# if branch support are percentages, convert them to probabilities

def degree(node):
    return (1 if not node.is_root() else 0) + node.num_children()

def dfs(l):
    # state: (leaf, distances)
    stack = [(l, [])]
    visited = set([l])
    results = {}
    while stack:
        u, lengths = stack.pop()
        if not u.is_leaf():
            # then we can go down
            for c in u.children:
                if c in visited:
                    continue
                if c.is_leaf():
                    results[str(c)] = fsum(lengths) + 1
                else:
                    stack.append((c, lengths + [1 + correction_map[degree(c)]]))
                visited.add(c)
        if not u.is_root():
            # going up
            if u.parent not in visited:
                if u.parent.is_root() and u.parent.num_children() == 2:
                    stack.append((u.parent, lengths))
                else:
                    stack.append((u.parent, lengths + [1 + correction_map[degree(u.parent)]]))
                visited.add(u.parent)
    return results

def get_distance(tree, ts2int):
    D = np.zeros((len(ts2int), len(ts2int)))
    for n in tree.traverse_leaves():
        distances = dfs(n)
        for u in distances:
            d = distances[u]
            D[ts2int[n.label], ts2int[u]] = d
            D[ts2int[u], ts2int[n.label]] = d
    return D

def get_distance_from_newick(tree, ts2int):
    return get_distance(ts.read_tree_newick(tree), ts2int)

def all_pairs_matrix(tree_, ts2int):
    tree = ts.read_tree_newick(tree_)
    N = len(ts2int)
    D = np.zeros((N, N))
    leaf_dists = dict()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            leaf_dists[node] = [[node,0]]
        else:
            calculated_root = False
            for c in node.children:
                if c.edge_length is not None:
                    if node.is_root() and node.num_children() == 2:
                        if calculated_root:
                            continue
                        calculated_root = True
                    for i in range(len(leaf_dists[c])):
                        l = leaf_dists[c][i][1]
                        leaf_dists[c][i][1] = l + 1 + correction_map[degree(node)]
            for c1 in range(0,len(node.children)-1):
                leaves_c1 = leaf_dists[node.children[c1]]
                for c2 in range(c1+1,len(node.children)):
                    leaves_c2 = leaf_dists[node.children[c2]]
                    for i in range(len(leaves_c1)):
                        for j in range(len(leaves_c2)):
                            u, d1 = leaves_c1[i]
                            v, d2 = leaves_c2[j]
                            d = d1 + d2 + correction_map[degree(node)]
                            u_key = u.label
                            v_key = v.label
                            estimated = d
                            D[ts2int[u_key], ts2int[v_key]] = estimated
                            D[ts2int[v_key], ts2int[u_key]] = estimated
            leaf_dists[node] = leaf_dists[node.children[0]]; del leaf_dists[node.children[0]]
            for i in range(1,len(node.children)):
                leaf_dists[node] += leaf_dists[node.children[i]]; del leaf_dists[node.children[i]]
    return D

def build_D(trees, n_jobs = -1):
    taxons = ad.get_ts(trees)
    # tsw_trees = [ for t in trees]
    ts2int = get_ts_mapping(ts.read_tree_newick(trees[0]))
    DMs = Parallel(n_jobs=n_jobs)(delayed(get_distance_from_newick)(trees[k], ts2int) for k in range(len(trees)))
    Ds = all_matrices(taxons, trees)
    for k in range(len(trees)):
        DM = DMs[k]
        D = Ds[k]
        for i, j in taxon_pairs(taxons):
            iname, jname = taxons[i], taxons[j]
            if i == j:
                D[i, j] = 0
                D.setmask((i, j), 1)
                continue
            D[i, j] = DM[ts2int[iname], ts2int[jname]]
            D.setmask((i, j), 1)
    return taxons, matrix_elementwise(taxons, Ds, np.mean)

def build_D2(trees, n_jobs = -1):
    taxons = ad.get_ts(trees)
    # tsw_trees = [ for t in trees]
    ts2int = get_ts_mapping(ts.read_tree_newick(trees[0]))
    DMs = Parallel(n_jobs=n_jobs)(delayed(all_pairs_matrix)(trees[k], ts2int) for k in range(len(trees)))
    Ds = all_matrices(taxons, trees)
    for k in range(len(trees)):
        DM = DMs[k]
        D = Ds[k]
        for i, j in taxon_pairs(taxons):
            iname, jname = taxons[i], taxons[j]
            if i == j:
                D[i, j] = 0
                D.setmask((i, j), 1)
                continue
            D[i, j] = DM[ts2int[iname], ts2int[jname]]
            D.setmask((i, j), 1)
    return taxons, matrix_elementwise(taxons, Ds, np.mean)

def get_ts_mapping(tree):
    result = {}
    leaves = list(tree.traverse_leaves())
    leaves = sorted(leaves, key=lambda x: x.label)
    for i, t in enumerate(leaves):
        result[t.label] = i
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required = True)
    parser.add_argument('-m', '--montecarlo', type=int, default = 0)
    parser.add_argument('-w', '--weighting', type=str, default = 's')
    parser.add_argument('--renormalize', action='store_true')
    parser.add_argument('-c', '--cores', type=int, default = -1)
    parser.add_argument('-o', '--output', type=str, default = "-")

    args = parser.parse_args()
    trees = open(args.input, "r").readlines()
    ts_trees = [ts.read_tree_newick(t) for t in trees]
    normalize(ts_trees, args.renormalize)
    if args.montecarlo > 0:
        start = time.time()
        trees = explode(ts_trees, args.montecarlo, args.weighting)
        end = time.time()
        print(f"Exploded all trees in {end - start} seconds")
    start = time.time()
    taxa, D = build_D(trees, args.cores)
    end = time.time()
    print(f"Built the distance matrix in {end - start} seconds")
    T = run_iterations(taxa, D, "s")
    if args.output == "-":
        print(T)
    else:
        with open(args.output, "w+") as fh:
            fh.write(T)