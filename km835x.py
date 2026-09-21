"""Exact-cover version of the Kramer-Mesner search (replaces the SAT step of km835g.py).
Rows = G-orbits of (k-1)-sets, columns = G-orbits of k-sets.  Sparse KM data: for each row i the column ids of the
v-k+1 supersets rep_i + {x}.  A design = set of columns such that every row appears exactly once.
Preprocessing:  (1) a column that hits some row twice is impossible;  (2) when v = 2k every S(k-1,k,2k) is closed
under complements (Hoffman-bound equality), so columns j and comp(j) are merged into one variable, impossible if their
row sets overlap;  (3) unit propagation: a row with one candidate forces it, a row with none proves nonexistence.
Then Algorithm X (exact cover) on the residue.  Usage: python km835x.py <group> <q> <k>"""
import sys, time
import numpy as np
from numba import njit
from km835 import binom_table, rank_of
from km835g import group, orbits

@njit(cache=True)
def km_sparse(reps_t, ids_k, v, k, B):
    R = np.empty((reps_t.shape[0], v - k + 1), np.int64)
    for i in range(reps_t.shape[0]):
        T = reps_t[i]; c = 0
        for x in range(v):
            if not (T >> x) & 1: R[i, c] = ids_k[rank_of(T | (1 << x), B)]; c += 1
    return R

def solve(R, ncols, comp, log):
    nrows = R.shape[0]
    rows_of = [set() for _ in range(ncols)]; bad = np.zeros(ncols, bool)
    for i in range(nrows):
        seen = set()
        for j in R[i]:
            j = int(j)
            if j in seen: bad[j] = True
            seen.add(j); rows_of[j].add(i)
    log(f"  columns hitting a row twice: {int(bad.sum())}")
    # merge complement pairs into variables
    var_rows = {}; nvar_pairs = 0
    for j in range(ncols):
        if bad[j]: continue
        c = int(comp[j]) if comp is not None else j
        if c < j: continue                       # handled from its partner
        if c == j: var_rows[j] = rows_of[j]
        else:
            if bad[c]: continue
            if rows_of[j] & rows_of[c]: continue  # the pair would cover a row twice
            var_rows[j] = rows_of[j] | rows_of[c]; nvar_pairs += 1
    log(f"  variables after complement merging: {len(var_rows)} ({nvar_pairs} pairs)")
    cand = {i: set() for i in range(nrows)}
    for vv, rs in var_rows.items():
        for i in rs: cand[i].add(vv)
    # unit propagation
    chosen = []; alive_rows = set(range(nrows)); alive_vars = set(var_rows)
    changed = True
    while changed:
        changed = False
        for i in list(alive_rows):
            cs = cand[i] & alive_vars
            if not cs: log(f"  row {i} has no candidate -> NO invariant Steiner system"); return None
            if len(cs) == 1:
                vv = next(iter(cs)); chosen.append(vv); changed = True
                covered = var_rows[vv] & alive_rows
                for r in covered:
                    for u in cand[r] & alive_vars: alive_vars.discard(u)     # every var touching a covered row dies
                alive_rows -= covered
                break
    log(f"  after propagation: {len(chosen)} forced variables, {len(alive_rows)} rows and {len(alive_vars)} variables left")
    # Algorithm X on the residue
    X = {i: cand[i] & alive_vars for i in alive_rows}
    Y = {vv: var_rows[vv] & alive_rows for vv in alive_vars}
    nodes = [0]
    def select(r):
        cols = []
        for i in Y[r]:
            for u in X[i]:
                for i2 in Y[u]:
                    if i2 != i: X[i2].remove(u)
            cols.append(X.pop(i))
        return cols
    def deselect(r, cols):
        for i in reversed(list(Y[r])):
            X[i] = cols.pop()
            for u in X[i]:
                for i2 in Y[u]:
                    if i2 != i: X[i2].add(u)
    def rec(sol):
        nodes[0] += 1
        if not X: return list(sol)
        i = min(X, key=lambda r: len(X[r]))
        for u in list(X[i]):
            sol.append(u); cols = select(u)
            res = rec(sol)
            deselect(u, cols); sol.pop()
            if res is not None: return res
        return None
    t = time.time(); res = rec([]); log(f"  Algorithm X: {nodes[0]} nodes, {time.time()-t:.1f}s")
    if res is None: log("  NO invariant Steiner system"); return None
    sol = chosen + res
    cols = []
    for vv in sol:
        cols.append(vv)
        if comp is not None and int(comp[vv]) != vv: cols.append(int(comp[vv]))
    return sorted(cols)

def main(name, q, k):
    lines = []
    def log(s): print(s, flush=True); lines.append(s)
    P, v = group(name, q); B = binom_table()
    log(f"S({k-1},{k},{v}) invariant under {name}({q}), |G| = {P.shape[0]}  [exact cover]")
    ids_k, reps_k, sizes_k = orbits(k, P, v, B, log)
    ids_t, reps_t, sizes_t = orbits(k - 1, P, v, B, log); del ids_t
    R = km_sparse(reps_t, ids_k, v, k, B)
    comp = None
    if 2 * k == v:
        full = (1 << v) - 1; comp = np.array([ids_k[rank_of(full ^ int(r), B)] for r in reps_k])
        assert (comp[comp] == np.arange(len(reps_k))).all()
    del ids_k
    np.savez(f"km_{name}{q}_k{k}.npz", R=R, reps_k=reps_k, sizes_k=sizes_k, reps_t=reps_t, comp=comp if comp is not None else np.zeros(0))
    import os
    if os.environ.get('KM_NOSOLVE'):
        log('  saved npz, solve skipped'); open(f"kmx_{name}{q}_k{k}.log", "w").write("\n".join(lines) + "\n"); return
    cols = solve(R, len(reps_k), comp, log)
    if cols is not None:
        blocks = int(sizes_k[cols].sum()); expect = int(B[v, k-1] // k)
        # independent check: every row covered exactly once
        cnt = np.zeros(R.shape[0], np.int64); S = set(cols)
        for i in range(R.shape[0]): cnt[i] = sum(1 for j in R[i] if int(j) in S)
        log(f"  SOLUTION: {len(cols)} orbits, {blocks} blocks (expected {expect}), rows covered once: {bool((cnt == 1).all())}")
        assert blocks == expect and (cnt == 1).all()
        np.save(f"sol_{name}{q}_k{k}.npy", np.array([int(reps_k[j]) for j in cols], np.int64))
    open(f"kmx_{name}{q}_k{k}.log", "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
