"""LP / MIP feasibility of the Kramer-Mesner exact-cover system saved by km835x.py (npz: R, comp, sizes_k).
Variables = complement-merged column orbits (as in km835x.py).  Constraints: every row covered exactly once.
LP relaxation infeasible  =>  no invariant Steiner system (rigorous up to LP tolerance, then re-checked with the
dual Farkas certificate in exact arithmetic is NOT done here: report as 'LP infeasible' only).
Usage: python km835lp.py <npz> [mip_time_limit_s]"""
import sys, time
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

f = sys.argv[1]; tl = float(sys.argv[2]) if len(sys.argv) > 2 else 600
d = np.load(f); R = d['R']; comp = d['comp']; sizes = d['sizes_k']
nrows, m = R.shape; ncols = len(sizes)
rows_of = [[] for _ in range(ncols)]; bad = np.zeros(ncols, bool)
for i in range(nrows):
    seen = set()
    for j in R[i]:
        j = int(j)
        if j in seen: bad[j] = True
        seen.add(j); rows_of[j].append(i)
vars_ = []; var_cols = []
for j in range(ncols):
    if bad[j]: continue
    c = int(comp[j]) if len(comp) else j
    if c < j: continue
    if c == j: vars_.append([j])
    else:
        if bad[c] or set(rows_of[j]) & set(rows_of[c]): continue
        vars_.append([j, c])
nv = len(vars_); print(f"{f}: {nrows} rows, {ncols} columns, {nv} variables", flush=True)
ri, ci = [], []
for t, cols in enumerate(vars_):
    for j in cols:
        for i in rows_of[j]: ri.append(i); ci.append(t)
A = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(nrows, nv))
empty = np.nonzero(np.diff(A.indptr) == 0)[0]
if len(empty): print(f"rows with no candidate: {len(empty)} -> infeasible"); sys.exit()
t = time.time()
res = linprog(np.zeros(nv), A_eq=A, b_eq=np.ones(nrows), bounds=(0, 1), method='highs')
print(f"LP relaxation: {res.message}  status {res.status}  ({time.time()-t:.1f}s)", flush=True)
if res.status == 2: print("LP INFEASIBLE -> no invariant Steiner system"); sys.exit()
x = res.x; print(f"  LP solution: {int((x > 1e-6).sum())} nonzero, {int((x > 1-1e-6).sum())} at 1, max fractional {x[(x>1e-6)&(x<1-1e-6)].max() if ((x>1e-6)&(x<1-1e-6)).any() else 0:.3f}")
t = time.time()
res = milp(np.zeros(nv), constraints=LinearConstraint(A, 1, 1), integrality=np.ones(nv), bounds=Bounds(0, 1), options=dict(time_limit=tl, disp=False))
print(f"MIP: {res.message}  status {res.status}  ({time.time()-t:.1f}s)", flush=True)
if res.status == 0:
    sel = [t for t in range(nv) if res.x[t] > 0.5]; cols = [j for t in sel for j in vars_[t]]
    print(f"  SOLUTION: {len(cols)} orbits, blocks {int(sizes[cols].sum())}")
    np.save(f.replace('km_', 'solmip_'), np.array(cols))
elif res.status == 2: print("MIP INFEASIBLE -> no invariant Steiner system")
