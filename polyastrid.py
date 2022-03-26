from copy import copy
from math import fsum
from random import random
import treeswift as ts
import numpy as np
from functools import lru_cache
import json
import argparse
import asterid as ad
import asterid as astrid
import os

__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

with open(os.path.join(__location__,"probs.json")) as fh: correction_map_ = json.load(fh)

correction_map = {}
for k in correction_map_:
    correction_map[int(k)] = correction_map_[k]
for i in range(3):
    correction_map[i] = 0

# if branch support are percentages, convert them to probabilities
def normalize(trees):
    maximum_value = 1
    for t in trees:
        for n in t.traverse_internal():
            if n.label:
                maximum_value = max(maximum_value, float(n.label))
    if maximum_value > 1:
        for t in trees:
            for n in t.traverse_internal():
                if n.label:
                    n.label = float(n.label) / 100
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

def all_matrices(ts, trees):
    return [ad.DistanceMatrix(ts, t) for t in trees]

def taxon_pairs(ts):
    for i in range(len(ts)):
        for j in range(i, len(ts)):
            yield i, j

def matrix_elementwise(ts, Ds, f):
    R = ad.DistanceMatrix(ts)
    for i, j in taxon_pairs(ts):
        if i == j:
            R[i, j] = 0
            R.setmask((i, j), len(Ds))
            continue
        l = []
        for D in Ds:
            if D.has(i, j):
                l.append(D[i, j])
        if l:
            R[i, j] = f(l)
            R.setmask((i, j), len(l))
        else:
            R[i, j] = 0
            R.setmask((i, j), 0)
    return R

def build_D(trees):
    taxons = ad.get_ts(trees)
    Ds = all_matrices(taxons, trees)
    tsw_trees = [ts.read_tree_newick(t) for t in trees]
    ts2int = get_ts_mapping(tsw_trees[0])
    for k in range(len(trees)):
        DM = get_distance(tsw_trees[k], ts2int)
        # print(DM.sum())
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

def monte_carlo_contract(tree_):
    tree = copy(tree_)
    to_contract = []
    for n in tree.traverse_internal():
        if n.is_root():
            continue
        if random() < (1 - float(n.label)):
            to_contract.append(n)
    for n in to_contract:
        n.contract()
    return tree

def run_iterations(ts, D, methods):
    fns = {
        "u": lambda ts, D: astrid.upgma_star(ts, D) + ";",
        "f": lambda ts, D: astrid.fastme_balme(ts, D, 0, 0),
        "n": lambda ts, D: astrid.fastme_balme(ts, D, 1, 0),
        "s": lambda ts, D: astrid.fastme_balme(ts, D, 1, 1),
        "j": lambda ts, D: astrid.fastme_nj(ts, D, 0, 0),
        "N": lambda ts, D: astrid.fastme_nj(ts, D, 1, 0),
        "S": lambda ts, D: astrid.fastme_nj(ts, D, 1, 1),
    }
    t = None
    for m in methods:
        if m not in fns:
            raise f"{m} not found as a method!"
        f = fns[m]
        t = f(ts, D)
        D.fill_in_transient(ts, t)
    return t

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required = True)
    parser.add_argument('-o', '--output', type=str, default = "-")

    args = parser.parse_args()
    trees = open(args.input, "r").readlines()
    taxa, D = build_D(trees)
    T = run_iterations(taxa, D, "s")
    if args.output == "-":
        print(T)
    else:
        with open(args.output, "w+") as fh:
            fh.write(T)
# trees = []
# with open("avian.tre") as fh:
#     for l in fh:
#         trees.append(ts.read_tree_newick(l))
# normalize(trees)
# for i in range(10):
#     tree = monte_carlo_contract(trees[0])
#     print(tree.newick())
#     ts = get_ts_mapping(tree)
#     D = get_distance(tree, ts)
#     print(ts)
#     print(D)
#     print(D.sum()/2)