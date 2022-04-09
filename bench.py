import timeit
from polybase import taxon_pairs
from xastrid import build_D, build_D2, import_trees

trees = import_trees("scratch/avian.tre")
def sum_distance(ts, D):
    return sum(D[i, j] for i, j in taxon_pairs(ts))

def experiment1():
    s1 = sum_distance(*build_D(trees, "s", "i"))
    print(s1)
    # s2 = sum_distance(*build_D2(trees, "s", "i"))
    # print(s1, s2, abs(s1 - s2))

def experiment2():
    s2 = sum_distance(*build_D2(trees, "s", "i"))
    print(s2)

if __name__ == '__main__':
    t = timeit.Timer('experiment1()', setup="from __main__ import experiment1")
    print(f"{t.timeit(number=3)} seconds to run experiment1")
    t2 = timeit.Timer('experiment2()', setup="from __main__ import experiment2")
    print(f"{t2.timeit(number=3)} seconds to run experiment2")