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

def all_matrices(ts, trees):
    return [ad.DistanceMatrix(ts, t) for t in trees]

def taxon_pairs(ts):
    for i in range(len(ts)):
        for j in range(i, len(ts)):
            yield i, j


def monte_carlo_contract(tree_):
    tree = copy(tree_)
    to_contract = []
    for n in tree.traverse_internal():
        if n.is_root() or n.label == '':
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
        T = T.replace('[&R]', '')
    T = T.strip()
    return T

def explode(trees, factor):
    results = []
    for t in trees:
        for _ in range(factor):
            results.append(to_newick(monte_carlo_contract(t)))
    return results