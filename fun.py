import timeit
import treeswift as ts
from polybase import taxon_pairs, calc_support
from xastrid import build_D, build_D2, import_trees, dfs
trees = import_trees("scratch/quartet.tre")
t = ts.read_tree_newick(trees[0])
lookup = t.label_to_node()
print(calc_support(lookup['1']))
print(dfs(lookup['1']))