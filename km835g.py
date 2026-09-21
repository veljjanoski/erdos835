"""Kramer-Mesner search for Steiner systems S(k-1,k,v) invariant under a given permutation group G on v points,
G given by generators and enumerated by closure (|G| up to ~1e5, stored as permutations).  Same method as km835.py
(orbits by marking, ranked tables), but the final system A x = 1 is solved with a SAT solver (exactly-one per row).
Groups: pgl2 q / psl2 q  (projective line, v = q+1);  agl1 / agammal1 (F_32 = F_2[t]/(t^5+t^2+1), v = 32).
Usage: python km835g.py <group> <q> <k>     e.g.  psl2 11 6  (S(5,6,12) must be found),  pgl2 31 16,  agl1 32 16."""
import sys, time
import numpy as np
from numba import njit, prange
from km835 import binom_table, rank_of, next_unassigned

def closure(gens, v):
    ident = tuple(range(v)); seen = {ident}; frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                c = tuple(h[g[x]] for x in range(v))
                if c not in seen: seen.add(c); nxt.append(c)
        frontier = nxt
    return np.array(sorted(seen), np.uint8)

def group(name, q):
    if name in ("pgl2", "psl2"):
        v = q + 1; INF = q
        def mob(f): return tuple(f(x) for x in range(v))
        def inv(x): return pow(x, q - 2, q)
        g = next(a for a in range(2, q) if all(pow(a, (q-1)//p, q) != 1 for p in {2, 3, 5, 7, 11, 13} if (q-1) % p == 0))
        add = mob(lambda x: INF if x == INF else (x + 1) % q)
        if name == "pgl2":
            mul = mob(lambda x: INF if x == INF else (g * x) % q)
            invm = mob(lambda x: INF if x == 0 else (0 if x == INF else inv(x)))
        else:
            mul = mob(lambda x: INF if x == INF else (g * g * x) % q)
            invm = mob(lambda x: INF if x == 0 else (0 if x == INF else (-inv(x)) % q))
        return closure([add, mul, invm], v), v
    if name in ("agl1", "agammal1"):
        v = 32
        def mult(x):                       # multiply by t in F_2[t]/(t^5+t^2+1)
            y = x << 1
            if y & 32: y ^= 0b100101
            return y
        mt = tuple(mult(x) for x in range(v)); tr = tuple(x ^ 1 for x in range(v))
        gens = [mt, tr]
        if name == "agammal1":
            def sq(x):                      # Frobenius x -> x^2 via repeated multiplication
                r = 0; a = x; b = x
                while b:
                    if b & 1: r ^= a
                    a = mult(a); b >>= 1
                return r
            gens.append(tuple(sq(x) for x in range(v)))
        return closure(gens, v), v
    raise ValueError(name)

@njit(parallel=True, cache=True)
def mark_orbit_g(ids, S, oid, P, v, B):
    N = P.shape[0]
    bits = np.empty(v, np.int64); w = 0
    for x in range(v):
        if (S >> x) & 1: bits[w] = x; w += 1
    stabs = np.zeros(N, np.int64)
    for a in prange(N):
        U = 0
        for q in range(w): U |= 1 << (P[a, bits[q]] ^ 0)
        ids[rank_of(U, B)] = oid
        if U == S: stabs[a] = 1
    return stabs.sum()

@njit(cache=True)
def next_unassigned32(ids, mask, w, v, B):
    limit = 1 << v
    while mask < limit:
        if ids[rank_of(mask, B)] == 0xFFFFFFFF: return mask
        c = mask & (-mask); r = mask + c
        mask = (((r ^ mask) >> 2) // c) | r
    return -1

def orbits(w, P, v, B, log):
    n_sets = int(B[v, w]); ids = np.full(n_sets, 0xFFFFFFFF, np.uint32)
    reps = []; sizes = []; G = P.shape[0]; mask = (1 << w) - 1; t = time.time()
    while True:
        mask = next_unassigned32(ids, mask, w, v, B)
        if mask < 0: break
        stab = mark_orbit_g(ids, mask, len(reps), P, v, B)
        assert G % stab == 0
        reps.append(mask); sizes.append(G // stab)
    assert sum(sizes) == n_sets, (sum(sizes), n_sets)
    log(f"  weight {w}: {len(reps)} orbits, {n_sets} sets, {time.time()-t:.1f}s")
    return ids, np.array(reps, np.int64), np.array(sizes, np.int64)

@njit(cache=True)
def km_rows(reps_t, ids_k, v, B, ncols):
    A = np.zeros((reps_t.shape[0], ncols), np.int32)
    for i in range(reps_t.shape[0]):
        T = reps_t[i]
        for x in range(v):
            if not (T >> x) & 1: A[i, ids_k[rank_of(T | (1 << x), B)]] += 1
    return A

def solve_sat(A, log, complement=None):
    from pysat.solvers import Cadical153
    rows, cols = A.shape
    forced0 = (A >= 2).any(0)
    log(f"  KM matrix {A.shape}, columns forced to 0: {int(forced0.sum())}, max entry {A.max()}")
    s = Cadical153(); ncl = 0; nv = cols
    for i in range(rows):
        lits = [int(j) + 1 for j in np.nonzero(A[i] == 1)[0] if not forced0[j]]
        if not lits: log(f"  row {i} has no admissible column -> no solution"); return None
        s.add_clause(lits); ncl += 1
        # at most one (Sinz sequential counter): aux r_1..r_{n-1}
        n = len(lits)
        if n > 1:
            r = [nv + a + 1 for a in range(n - 1)]; nv += n - 1
            s.add_clause([-lits[0], r[0]]); ncl += 1
            for a in range(1, n - 1):
                s.add_clause([-lits[a], r[a]]); s.add_clause([-r[a-1], r[a]]); s.add_clause([-lits[a], -r[a-1]]); ncl += 3
            s.add_clause([-lits[n-1], -r[n-2]]); ncl += 1
    for j in np.nonzero(forced0)[0]: s.add_clause([-(int(j) + 1)]); ncl += 1
    if complement is not None:
        for j in range(cols):
            c = int(complement[j])
            if c != j: s.add_clause([-(j + 1), c + 1]); s.add_clause([j + 1, -(c + 1)]); ncl += 2
    t = time.time(); ok = s.solve(); log(f"  SAT: {ncl} clauses, {'SAT' if ok else 'UNSAT'} in {time.time()-t:.1f}s")
    return [j for j in range(cols) if s.get_model()[j] > 0] if ok else None

def main(name, q, k):
    lines = []
    def log(s): print(s, flush=True); lines.append(s)
    P, v = group(name, q); B = binom_table()
    log(f"S({k-1},{k},{v}) invariant under {name}({q}), |G| = {P.shape[0]}")
    ids_k, reps_k, sizes_k = orbits(k, P, v, B, log)
    ids_t, reps_t, sizes_t = orbits(k - 1, P, v, B, log)
    A = km_rows(reps_t, ids_k, v, B, len(reps_k)); assert (A.sum(1) == v - k + 1).all()
    comp = None
    if 2 * k == v:
        full = (1 << v) - 1; comp = np.array([ids_k[rank_of(full ^ int(r), B)] for r in reps_k])
    sol = solve_sat(A, log)
    if sol is not None:
        blocks = int(sizes_k[sol].sum()); expect = int(B[v, k-1] // k)
        cc = (all(int(comp[j]) in set(sol) for j in sol)) if comp is not None else None
        log(f"  SOLUTION: {len(sol)} orbits, total blocks {blocks} (expected {expect}), complement-closed {cc}")
        assert blocks == expect
        np.save(f"sol_{name}{q}_k{k}.npy", np.array([int(reps_k[j]) for j in sol], np.int64))
    else:
        log("  no invariant Steiner system")
    open(f"kmg_{name}{q}_k{k}.log", "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
