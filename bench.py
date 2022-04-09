import timeit
from polybase import taxon_pairs
from xastrid import build_D, build_D2, import_trees

trees = import_trees("scratch/avian.tre")

def sum_distance(ts, D):
    return sum(D[i, j] for i, j in taxon_pairs(ts))

def experiment1():
    ts, D = build_D(trees, "s", "i")
    # print(D[ts["1"], ts["4"]],D[ts["1"], ts["2"]],D[ts["1"], ts["3"]])
    # print(D.str())
    s1 = sum_distance(ts, D)
    print(s1)

def experiment2():
    ts, D = build_D2(trees, "s", "i")
    # print(D[ts["1"], ts["4"]],D[ts["1"], ts["2"]],D[ts["1"], ts["3"]])
    # print(D.str())
    s2 = sum_distance(ts, D)
    print(s2)

if __name__ == '__main__':
    t = timeit.Timer('experiment1()', setup="from __main__ import experiment1")
    print(f"{t.timeit(number=3)} seconds to run experiment1")
    t2 = timeit.Timer('experiment2()', setup="from __main__ import experiment2")
    print(f"{t2.timeit(number=3)} seconds to run experiment2")