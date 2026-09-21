# k=4 check of the reduction: chi(J(8,4)) = 5  <=>  the 70 four-subsets of [8] partition into 5 Steiner systems S(3,4,8).
# Enumerate all S(3,4,8) on [8] (14 blocks, every 3-set in exactly one block) by backtracking, then find the maximum
# number of pairwise disjoint systems.  Expected: 30 systems, max disjoint < 5.
import itertools, sys
from itertools import combinations
V=range(8); blocks=[frozenset(b) for b in combinations(V,4)]; triples=[frozenset(t) for t in combinations(V,3)]
tri_of={b:[frozenset(t) for t in combinations(b,3)] for b in blocks}
systems=[]
def bt(chosen, covered):
    if len(chosen)==14: systems.append(frozenset(chosen)); return
    # smallest uncovered triple must be covered by some block
    t=next(t for t in triples if t not in covered)
    for b in blocks:
        if t<=b and not any(x in covered for x in tri_of[b]):
            bt(chosen+[b], covered|set(tri_of[b]))
bt([],set())
print("S(3,4,8) systems on [8]:",len(systems))
# every system is complement-closed?
print("all complement-closed:", all(all(frozenset(V)-b in s for b in s) for s in systems))
# max number of pairwise disjoint systems (clique search on the disjointness graph)
n=len(systems); adj=[[i!=j and not (systems[i]&systems[j]) for j in range(n)] for i in range(n)]
best=0
def clique(cur, cand):
    global best
    best=max(best,len(cur))
    for i in cand:
        clique(cur+[i], [j for j in cand if j>i and adj[i][j]])
clique([], list(range(n)))
print("max pairwise disjoint S(3,4,8):",best,"(need 5 for a 5-colouring of J(8,4))")
