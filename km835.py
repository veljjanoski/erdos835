"""Kramer-Mesner search for Steiner systems S(k-1, k, 2^n) invariant under the affine group AGL(n,2) acting on the
points F_2^n (Erdős #835: a 17-colouring of J(32,16) needs 17 disjoint S(15,16,32); this looks for ONE such system
with prescribed symmetry).  A G-invariant design is a union of G-orbits on k-sets; it is a Steiner system iff every
(k-1)-set lies in exactly one chosen block.  Kramer-Mesner: rows = G-orbits on (k-1)-sets, columns = G-orbits on
k-sets, A[i][j] = number of blocks in orbit j containing a fixed representative of orbit i; need 0/1 x with A x = 1.
Orbits are enumerated by marking: every set of the given weight is ranked (combinatorial number system) into a table.
Usage: python km835.py n k        e.g. 3 4 (S(3,4,8): must be found), 4 4 (S(3,4,16): must be found), 5 16 (target)."""
import sys, time
import numpy as np
from numba import njit, prange

def binom_table(N=33):
    B = np.zeros((N, N), np.int64)
    for i in range(N):
        B[i, 0] = 1
        for j in range(1, i + 1): B[i, j] = B[i-1, j-1] + B[i-1, j]
    return B

@njit(cache=True)
def is_invertible(rows, n):
    r = rows.copy(); rank = 0
    for col in range(n):
        piv = -1
        for i in range(rank, n):
            if (r[i] >> col) & 1: piv = i; break
        if piv < 0: return False
        t = r[piv]; r[piv] = r[rank]; r[rank] = t
        for i in range(n):
            if i != rank and (r[i] >> col) & 1: r[i] ^= r[rank]
        rank += 1
    return True

@njit(parallel=True, cache=True)
def gl_flags(n):
    total = 1 << (n * n); flags = np.zeros(total, np.bool_)
    for m in prange(total):
        rows = np.empty(n, np.int64)
        for i in range(n): rows[i] = (m >> (n * i)) & ((1 << n) - 1)
        flags[m] = is_invertible(rows, n)
    return flags

@njit(parallel=True, cache=True)
def gl_perms(idx, n):
    N = idx.shape[0]; v = 1 << n; P = np.empty((N, v), np.uint8)
    for a in prange(N):
        m = idx[a]
        for x in range(v):
            y = 0
            for i in range(n):
                r = (m >> (n * i)) & (v - 1)
                c = 0; z = r & x
                while z: c ^= 1; z &= z - 1
                y |= c << i
            P[a, x] = y
    return P

@njit(cache=True)
def rank_of(mask, B):
    r = 0; j = 0; p = 0
    while mask:
        if mask & 1:
            j += 1; r += B[p, j]
        mask >>= 1; p += 1
    return r

@njit(cache=True)
def next_unassigned(ids, mask, w, v, B):
    """Gosper iteration from mask (inclusive) over weight-w masks below 2^v; return first unassigned mask or -1."""
    limit = 1 << v
    while mask < limit:
        if ids[rank_of(mask, B)] == 65535: return mask
        c = mask & (-mask); r = mask + c
        mask = (((r ^ mask) >> 2) // c) | r
    return -1

@njit(parallel=True, cache=True)
def mark_orbit(ids, S, oid, P, v, B):
    """Mark the AGL-orbit of S (image = {M x xor b}) with oid; return |stabilizer| (elements fixing S)."""
    N = P.shape[0]
    bits = np.empty(v, np.int64); w = 0
    for x in range(v):
        if (S >> x) & 1: bits[w] = x; w += 1
    stabs = np.zeros(N, np.int64)
    for a in prange(N):
        loc = 0
        for b in range(v):
            U = 0
            for q in range(w): U |= 1 << (P[a, bits[q]] ^ b)
            ids[rank_of(U, B)] = oid
            if U == S: loc += 1
        stabs[a] = loc
    return stabs.sum()

def orbits(w, P, v, B, log):
    n_sets = int(B[v, w]); ids = np.full(n_sets, 65535, np.uint16)
    reps = []; sizes = []; G = P.shape[0] * v
    mask = (1 << w) - 1; t = time.time()
    while True:
        mask = next_unassigned(ids, mask, w, v, B)
        if mask < 0: break
        oid = len(reps); stab = mark_orbit(ids, mask, oid, P, v, B)
        assert G % stab == 0
        reps.append(mask); sizes.append(G // stab)
        if len(reps) >= 65535: raise RuntimeError("too many orbits for uint16")
    assert sum(sizes) == n_sets, (sum(sizes), n_sets)     # orbits partition all weight-w sets
    log(f"  weight {w}: {len(reps)} orbits, {n_sets} sets, {time.time()-t:.1f}s")
    return ids, np.array(reps, np.int64), np.array(sizes, np.int64)

def km_matrix(reps_t, ids_k, v, B):
    A = np.zeros((len(reps_t), int(ids_k.max()) + 1), np.int64)
    for i, T in enumerate(reps_t):
        for x in range(v):
            if not (T >> x) & 1: A[i, ids_k[rank_of(T | (1 << x), B)]] += 1
    return A

def exact_cover(A, forced_zero):
    """All 0/1 x with A x = 1 (A entries 0/1 after forcing).  Algorithm X on small matrices; returns list of solutions."""
    rows, cols = A.shape; sols = []
    cols_ok = [j for j in range(cols) if not forced_zero[j]]
    def rec(uncovered, chosen):
        if not uncovered: sols.append(sorted(chosen)); return
        i = min(uncovered, key=lambda r: sum(1 for j in cols_ok if A[r, j]))
        for j in cols_ok:
            if A[i, j] and all(r in uncovered for r in range(rows) if A[r, j]):
                rec(uncovered - {r for r in range(rows) if A[r, j]}, chosen + [j])
    rec(frozenset(range(rows)), [])
    return sols

def main(n, k):
    v = 1 << n; B = binom_table()
    lines = []
    def log(s): print(s, flush=True); lines.append(s)
    log(f"S({k-1},{k},{v}) invariant under AGL({n},2)")
    t = time.time(); flags = gl_flags(n); idx = np.nonzero(flags)[0].astype(np.int64); P = gl_perms(idx, n)
    log(f"  |GL({n},2)| = {len(idx)}, |AGL| = {len(idx)*v}, {time.time()-t:.1f}s")
    ids_k, reps_k, sizes_k = orbits(k, P, v, B, log)
    ids_t, reps_t, sizes_t = orbits(k - 1, P, v, B, log)
    A = km_matrix(reps_t, ids_k, v, B)
    # each row: sum_j A[i,j] * (blocks in orbit j) ... consistency: sum_j A[i,j] = v-(k-1) supersets of the rep
    assert (A.sum(1) == v - k + 1).all()
    full = (1 << v) - 1
    comp = np.array([ids_k[rank_of(full ^ int(r), B)] for r in reps_k])       # complement permutes k-orbits (k = v/2)
    forced0 = (A >= 2).any(0)
    log(f"  KM matrix {A.shape}, columns forced to 0 (some entry >= 2): {int(forced0.sum())}, max entry {A.max()}")
    sols = exact_cover(A, forced0)
    log(f"  exact-cover solutions: {len(sols)}")
    for s in sols:
        blocks = int(sizes_k[s].sum()); expect = int(B[v, k-1] // k)
        cc = all(comp[j] in s for j in s)
        log(f"    orbits {s} sizes {[int(sizes_k[j]) for j in s]} total blocks {blocks} (expected {expect}) complement-closed {cc}")
        assert blocks == expect
    np.savez(f"km_n{n}_k{k}.npz", A=A, reps_k=reps_k, sizes_k=sizes_k, reps_t=reps_t, sizes_t=sizes_t, comp=comp)
    open(f"km_n{n}_k{k}.log", "w").write("\n".join(lines) + "\n")
    return sols

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
