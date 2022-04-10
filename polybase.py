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
import sys
from math import exp

def normalize(trees, mabayes = False):
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
    if mabayes:
        for t in trees:
            for n in t.traverse_internal():
                if n.label:
                    n.label = (float(n.label) - 0.333) * 3 / 2
    # for t in trees:
    #     if t.root.num_children() == 2:
    #         support = -1
    #         for c in t.root.children:
    #             if c.is_leaf():
    #                 support = 1
    #             elif c.label:
    #                 support = max(support, float(c.label))
    #         for c in t.root.children:
    #             c.label = support

def all_matrices(ts, trees):
    return [ad.DistanceMatrix(ts, t) for t in trees]

def taxon_pairs(ts):
    for i in range(len(ts)):
        for j in range(i, len(ts)):
            yield i, j

def is_fake_node(node):
    return node.num_children() == 2 and node.is_root()

def monte_carlo_contract(tree_, mode):
    tree = copy(tree_)
    to_contract = []
    considered_one_part = False
    for n in tree.traverse_internal():
        if n.is_root() or n.label == '' or not n.label:
            continue
        if considered_one_part and is_fake_node(n.parent):
            continue
        if is_fake_node(n.parent):
            considered_one_part = True
        if random() > calc_weight(n, mode):
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

def to_newick(tree):
    T = tree.newick()
    if '[&R]' in T:
        T = T.replace('[&R]', '', 1)
    T = T.strip()
    return T

def explode(trees, factor, mode = 's'):
    results = []
    for t in trees:
        for _ in range(factor):
            results.append(to_newick(monte_carlo_contract(t, mode)))
    return results

def calc_support(node):
    if node.is_leaf():
        return 1
    if is_fake_node(node.parent) and not node.label:
        for c in node.parent.children:
            if c.label and not c.is_leaf():
                return float(c.label)
            elif c.is_leaf():
                return 1
    if not node.label:
        return 0 # turns out aBayes does attribute things so we might need to impute missing data
    return float(node.label)

def calc_length(node):
    if node.parent.is_root() and node.parent.num_children() == 2:
        return sum(c.edge_length for c in node.parent.children)
    return node.edge_length

def calc_weight(node, mode = "s"):
    if mode == 'i':
        return 1
    if mode == "s":
        return calc_support(node)
    if mode == "l":
        return 1 - exp(-calc_length(node))
    if mode == "h":
        return calc_weight(node, "s") * calc_weight(node, "l")
    if mode == "h2":
        s = calc_support(node)
        return 1 - (1 - s) * exp(-calc_length(node))
    if mode == "hu":
        s = calc_support(node)
        return max(0, s - exp(-calc_length(node)))