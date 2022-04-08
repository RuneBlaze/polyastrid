from math import fsum, exp
from polybase import *

def postprocess(length, weight, mode):
    if mode == "i":
        return weight
    if mode == "l":
        return (1 - exp(-length)) * weight

def dfs(l, mode = "s", postprocessing = "i"):
    # state: (leaf, (total length, total weight))
    stack = [(l, ([], []))]
    visited = set([l])
    results = {}
    while stack:
        u, (lengths, weights) = stack.pop()
        if not u.is_leaf():
            # going down
            for c in u.children:
                if c in visited:
                    continue
                if c.is_leaf():
                    tt_lengths = fsum(lengths + [calc_length(c)])
                    tt_weights = fsum(weights + [calc_weight(c, mode)])
                    results[str(c)] = postprocess(tt_lengths, tt_weights, postprocessing)
                else:
                    stack.append((c, (lengths + [calc_length(c)], weights + [calc_weight(c, mode)])))
                visited.add(c)
        if not u.is_root():
            # going up
            if u.parent not in visited:
                if u.parent.is_root() and u.parent.num_children() == 2:
                    stack.append((u.parent, (lengths, weights)))
                else:
                    stack.append((u.parent, (lengths + [calc_length(u)], weights + [calc_weight(u, mode)])))
                visited.add(u.parent)
    return results

def get_distance(tree, ts2int, mode = "s", postprocessing = "i"):
    D = np.zeros((len(ts2int), len(ts2int)))
    for n in tree.traverse_leaves():
        distances = dfs(n, mode, postprocessing)
        for u in distances:
            d = distances[u]
            D[ts2int[n.label], ts2int[u]] = d
            D[ts2int[u], ts2int[n.label]] = d
    return D

def get_ts_mapping(tree):
    result = {}
    leaves = list(tree.traverse_leaves())
    leaves = sorted(leaves, key=lambda x: x.label)
    for i, t in enumerate(leaves):
        result[t.label] = i
    return result

def build_D(trees, mode, postprocessing):
    taxons = ad.get_ts(trees)
    Ds = all_matrices(taxons, trees)
    tsw_trees = [ts.read_tree_newick(t) for t in trees]
    ts2int = get_ts_mapping(tsw_trees[0])
    for k in range(len(trees)):
        DM = get_distance(tsw_trees[k], ts2int, mode, postprocessing)
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required = True)
    parser.add_argument('-w', '--weighting', type=str, default = 's')
    parser.add_argument('-p', '--postprocessing', type=str, default = 'i')
    parser.add_argument('-o', '--output', type=str, default = "-")
    args = parser.parse_args()
    trees = open(args.input, "r").readlines()
    ts_trees = [ts.read_tree_newick(t) for t in trees]
    normalize(ts_trees)
    taxa, D = build_D(trees, args.weighting, args.postprocessing)
    T = run_iterations(taxa, D, "s")
    if args.output == "-":
        print(T)
    else:
        with open(args.output, "w+") as fh:
            fh.write(T)