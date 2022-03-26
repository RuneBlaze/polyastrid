import asterid as astrid
import treeswift as ts
import argparse
from polybase import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required = True)
    parser.add_argument('-m', '--montecarlo', type=int, default = 100)
    parser.add_argument('-o', '--output', type=str, default = "-")
    args = parser.parse_args()

    with open(args.input, "r") as fh: trees = fh.readlines()
    tsw_trees = [ts.read_tree_newick(t) for t in trees]
    normalize(tsw_trees)
    trees = explode(tsw_trees, args.montecarlo)
    taxa = astrid.get_ts(trees)
    D = astrid.mk_distance_matrix(taxa, trees)
    T = run_iterations(taxa, D, "s")
    if args.output == "-":
        print(T)
    else:
        with open(args.output, "w+") as fh:
            fh.write(T)