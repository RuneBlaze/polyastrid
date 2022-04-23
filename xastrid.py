from math import fsum
from tkinter import X
from numpy import exp, log
from polybase import *

def normalization_factor(tree, ts2int, strategy):
    N = len([t for t in tree.traverse_leaves()])
    if strategy == "none":
        return 1
    if strategy.endswith("path"):
        D, _ = all_pairs_matrix(tree, ts2int, "L", "i")
        if strategy == "longest_path":
            return D.max() / (0.25 * N)
        if strategy == "longest_internal_path":
            int2ts = {v: k for k, v in ts2int.items()}
            trans = tree.label_to_node()
            i = np.argmax(D)
            x, y = np.unravel_index(i, D.shape)
            u = int2ts[x]
            v = int2ts[y]
            n1 = trans[u]
            n2 = trans[v]
            r = D.max() - calc_length(n1) - calc_length(n2)
            return r / (0.25 * N)
    if strategy == "longest_internal_edge":
        return max(calc_length(i) for i in tree.traverse_internal() if not i.is_root())
    if strategy == "longest_edge":
        return max(calc_length(i) for i in tree.traverse_postorder() if not i.is_root())
    if strategy == "total_length":
        return sum(calc_length(i) for i in tree.traverse_postorder() if not i.is_root()) / N
    if strategy == "total_internal_length":
        return sum(calc_length(i) for i in tree.traverse_internal() if not i.is_root()) / (0.5 * N)

def normalize_tree(tree, ts2int, strategy, normalize_leaf):
    f = normalization_factor(tree, ts2int, strategy) * 0.5
    if f <= 0:
        return
    for n in tree.traverse_internal():
        if n.is_root():
            continue
        if n.edge_length:
            n.edge_length = float(n.edge_length) / f
    for n in tree.traverse_leaves():
        if normalize_leaf:
            n.edge_length = 1
        else:
            n.edge_length = float(n.edge_length) / f

def postprocess(prodsupp, weight, mode):
    if mode == "i":
        return weight
    if mode == "m":
        return prodsupp * weight


# def pseudoglass(tree, ts2int):
#     N = len(ts2int)
#     D = np.zeros((N, N)) # the average coalescence time matrix
#     leaf_dists = dict()
#     for node in tree.traverse_postorder():
#         if node.is_leaf():
#             leaf_dists[node] = [[node,(0, 0, 0)]]
#         else:
#             calculated_root = False
#             for c in node.children:
#                 if c.edge_length is not None:
#                     if node.is_root() and node.num_children() == 2:
#                         if calculated_root:
#                             continue
#                         calculated_root = True
#                     for i in range(len(leaf_dists[c])):
#                         dist, u, v = leaf_dists[c][i][1]
#                         u_acc = 0
#                         leaf_dists[c][i][1] = (dist + 1, u + u_acc, v + calc_length(c))
#             for c1 in range(0,len(node.children)-1):
#                 leaves_c1 = leaf_dists[node.children[c1]]
#                 for c2 in range(c1+1,len(node.children)):
#                     leaves_c2 = leaf_dists[node.children[c2]]
#                     for i in range(len(leaves_c1)):
#                         for j in range(len(leaves_c2)):
#                             u, (ub, ul, ud) = leaves_c1[i]
#                             v, (vb, vl, vd) = leaves_c2[j]
#                             l, d = (ul + vl, ud + vd)
#                             u_key = u.label
#                             v_key = v.label
#                             estimated = d
#                             weight = 1
#                             if postprocessing == "as":
#                                 weight = l / (ub + vb)
#                             elif postprocessing == "was":
#                                 if estimated <= 0:
#                                     weight = 0
#                                 else:
#                                     weight = l / estimated
#                             elif postprocessing == "naive":
#                                 weight = exp(-l)
#                             tu = ts2int[u_key]
#                             tv = ts2int[v_key]
#                             D[tu, tv] = estimated
#                             D[tv, tu] = estimated
#                             W[tu, tv] = weight
#                             W[tv, tu] = weight
#             leaf_dists[node] = leaf_dists[node.children[0]]; del leaf_dists[node.children[0]]
#             for i in range(1,len(node.children)):
#                 leaf_dists[node] += leaf_dists[node.children[i]]; del leaf_dists[node.children[i]]
#     return D, W

def all_pairs_matrix(tree, ts2int, mode, postprocessing, ct_mat = None, ct_normalizer = None):
    N = len(ts2int)
    D = np.zeros((N, N))
    W = np.zeros((N, N))
    leaf_dists = dict()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            leaf_dists[node] = [[node,(0, 0, 0)]]
        else:
            calculated_root = False
            for c in node.children:
                if c.edge_length is not None:
                    if node.is_root() and node.num_children() == 2:
                        if calculated_root:
                            continue
                        calculated_root = True
                    for i in range(len(leaf_dists[c])):
                        dist, u, v = leaf_dists[c][i][1]
                        u_acc = 0
                        if postprocessing == "as":
                            u_acc = calc_support(c)
                        elif postprocessing == "was":
                            u_acc = calc_support(c) * calc_length(c)
                        else:
                            u_acc = calc_length(c)
                        leaf_dists[c][i][1] = (dist + 1, u + u_acc, v + calc_weight(c, mode))
            for c1 in range(0,len(node.children)-1):
                leaves_c1 = leaf_dists[node.children[c1]]
                for c2 in range(c1+1,len(node.children)):
                    leaves_c2 = leaf_dists[node.children[c2]]
                    for i in range(len(leaves_c1)):
                        for j in range(len(leaves_c2)):
                            u, (ub, ul, ud) = leaves_c1[i]
                            v, (vb, vl, vd) = leaves_c2[j]
                            l, d = (ul + vl, ud + vd)
                            u_key = u.label
                            v_key = v.label
                            estimated = d
                            weight = 1
                            tu = ts2int[u_key]
                            tv = ts2int[v_key]
                            if postprocessing == "as":
                                weight = l / (ub + vb)
                            elif postprocessing == "was":
                                if estimated <= 0:
                                    weight = 0
                                else:
                                    weight = l / estimated
                            elif postprocessing == "naive":
                                weight = exp(-l)
                            elif postprocessing == "naive2":
                                n = ct_normalizer#[tu, tv]
                                if n <= 0:
                                    n = 1
                                weight = exp(-l/n)
                                # if weight == 0:
                                #     print(l, ct_normalizer[tu, tv])
                                # print(weight, exp(-l))
                            elif postprocessing == "msc":
                                n = ct_normalizer#[tu, tv]
                                if n <= 0:
                                    n = 1
                                cu = (l - ct_mat[tu, tv]) / n
                                cu = max(cu, 0)
                                weight = exp(-cu)
                            elif postprocessing == "msc2":
                                cu = (l - ct_mat[tu, tv])
                                cu = max(cu, 0)
                                weight = exp(-cu)
                            D[tu, tv] = estimated
                            D[tv, tu] = estimated
                            W[tu, tv] = weight
                            W[tv, tu] = weight
            leaf_dists[node] = leaf_dists[node.children[0]]; del leaf_dists[node.children[0]]
            for i in range(1,len(node.children)):
                leaf_dists[node] += leaf_dists[node.children[i]]; del leaf_dists[node.children[i]]
    return D, W

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

def median_of_means(arr, buckets):
    """compute the median of means of arr with number of buckets"""
    np.random.shuffle(arr)
    buckets = np.array_split(arr, buckets)
    # Compute the mean within each bucket
    b_means = [np.mean(x) for x in buckets]
    # Compute the median-of-means
    median = np.median(np.array(b_means))
    return median

def build_D2(trees, mode, postprocessing, norm_strategy, keep_leaves):
    taxons = ad.get_ts(trees)
    Ds = all_matrices(taxons, trees)
    Ws = all_matrices(taxons, trees)
    tsw_trees = [ts.read_tree_newick(t) for t in trees]
    ts2int = get_ts_mapping(tsw_trees[0])
    for t in tsw_trees:
        normalize_tree(t, ts2int, norm_strategy, keep_leaves)
    # initiate the GLASS-like matrices
    N = len(ts2int)
    minDis = None
    avgDis = None
    rateNorm = None
    rateMed = None
    if postprocessing in ['msc', 'msc2', 'naive2']:
        minDis = np.zeros((N, N))
        avgDis = np.zeros((N, N))
        first_element = True
        for tree in tsw_trees:
            dis_mat, w = all_pairs_matrix(tree, ts2int, 'L', 'i')
            avgDis += dis_mat
            # weights += w
            if first_element:
                minDis = dis_mat
                first_element = False
            else:
                minDis = np.minimum(minDis, dis_mat)
        avgDis /= len(tsw_trees)
        rateNorm = (avgDis - minDis) / 2
        rates = []
        for i in range(N):
            for j in range(i+1,N):
                rates.append(rateNorm[i,j])
        rateMed = median_of_means(np.asarray(rates), 10)
        print("rateMed:", rateMed)
    for k in range(len(trees)):
        DM, WM = all_pairs_matrix(tsw_trees[k], ts2int, mode, postprocessing, minDis, rateMed)
        D = Ds[k]
        W = Ws[k]
        for i, j in taxon_pairs(taxons):
            iname, jname = taxons[i], taxons[j]
            if i == j:
                D[i, j] = 0
                D.setmask((i, j), 1)
                continue
            D[i, j] = DM[ts2int[iname], ts2int[jname]]
            W[i, j] = WM[ts2int[iname], ts2int[jname]]
            D.setmask((i, j), 1)
            W.setmask((i, j), 1)
    return taxons, matrix_weightedaverage(taxons, Ds, Ws)

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

def import_trees(path):
    trees = open(path, "r").readlines()
    ts_trees = [ts.read_tree_newick(t) for t in trees]
    normalize(ts_trees)
    return [to_newick(t) for t in ts_trees]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required = True)
    parser.add_argument('-w', '--weighting', type=str, default = 's')
    parser.add_argument('-p', '--postprocessing', type=str, default = 'i')
    parser.add_argument('--renormalize', action='store_true')
    parser.add_argument('-s', "--strategy", type=str, default = 'none')
    parser.add_argument('--eachtree', action='store_true')
    parser.add_argument('-o', '--output', type=str, default = "-")
    args = parser.parse_args()
    trees = open(args.input, "r").readlines()
    ts_trees = [ts.read_tree_newick(t) for t in trees]
    if args.renormalize:
        normalize(ts_trees)
    taxa, D = build_D2([to_newick(t) for t in ts_trees], args.weighting, args.postprocessing, args.strategy, args.renormalize)
    T = run_iterations(taxa, D, "s")
    if args.output == "-":
        print(T)
    else:
        with open(args.output, "w+") as fh:
            fh.write(T)