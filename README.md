# Erdős problem #835, case k = 16: Steiner systems S(15,16,32) with prescribed symmetry

**Problem.** Is there a k > 2 such that the k-subsets of {1, …, 2k} can be coloured with k+1 colours so that
every (k+1)-subset contains k-subsets of all k+1 colours? Equivalently: is the chromatic number of the Johnson
graph J(2k,k) equal to k+1? (https://www.erdosproblems.com/835). Known: false for 3 ≤ k ≤ 8 by computation;
false whenever k+1 is not prime (Ma–Tang); false for k = 10 and k = 12, because a colour class would be an
S(9,10,20) or an S(11,12,24) (see the reduction below), and deriving it five or seven times gives an S(4,5,15) or
an S(4,5,17), neither of which exists (Mendelsohn–Hung, Utilitas Math. 1 (1972); Östergård–Pottonen, J. Combin.
Theory Ser. A 115 (2008); remark by athvedt on the problem page). So the smallest open case is k = 16.

**Reduction (already noted on the problem page).** In a (k+1)-colouring, the k+1 sets T ∪ {x} through a fixed
(k−1)-set T form a clique, so each colour class contains exactly one k-set through every (k−1)-set: each colour
class is a Steiner system S(k−1, k, 2k). Conversely such systems are independent sets. Hence a 17-colouring of
J(32,16) exists iff the 16-subsets of a 32-set can be partitioned into 17 Steiner systems S(15,16,32)
(35,357,670 blocks each). No Steiner system with t ≥ 6 is known; deriving eleven times turns S(15,16,32) into an
S(4,5,21), whose existence is open.

**What this repository does.** It asks the weaker question whether a single S(15,16,32) exists with a prescribed
automorphism group G ≤ Sym(32), by the Kramer–Mesner method: a G-invariant system is a union of G-orbits of 16-sets,
and it is a Steiner system iff every G-orbit of 15-sets is covered exactly once, which is a 0/1 linear system
(rows = orbits of 15-sets, columns = orbits of 16-sets, entry = number of blocks of the column orbit containing
a fixed representative of the row orbit).

## Results

| group G on 32 points | order | orbits on 16-sets / 15-sets | result |
|---|---|---|---|
| AGL(5,2) | 319,979,520 | 38 / 30 | **no invariant S(15,16,32)** (exhaustive exact cover, `km835.py`) |
| PGL(2,31) | 29,760 | 20,636 / 19,234 | **no invariant S(15,16,32)**: some orbit of 15-sets has no admissible column at all, i.e. every orbit of 16-sets containing a superset of its representative also contains two sets sharing 15 points (`km835g.py`, `km835x.py`) |
| PSL(2,31) | 14,880 | 40,843 / 38,039 | undecided: 21,213 variables after merging complementary orbits; CaDiCaL ran 6 CPU-hours without a verdict, Algorithm X and HiGHS LP/MIP did not finish |
| AΓL(1,32) | 4,960 | 121,282 / 114,065 | undecided: 60,714 variables, system saved (`km_agammal132_k16.npz`) |
| AGL(1,32) | 992 | 606,330 / 570,285 | undecided: 304,174 variables, system regenerable in 3 minutes (`KM_NOSOLVE=1 python km835x.py agl1 32 16`) |

The saved LP runs do not decide the three undecided groups. For AGL(1,32) the interior-point method on the LP
relaxation reached primal infeasibility 3.55·10⁻¹⁵ before HiGHS ran out of memory and returned no solution
(`lp_agl1_ipm.out`); this suggests that the relaxation is feasible, but no feasible point was saved or checked.
The saved PSL(2,31) log ends after 4 interior-point iterations with no result (`lp_psl2_ipm.out`), and the
AΓL(1,32) log holds only the problem size (`lp_agammal1.out`). The results above say nothing about
S(15,16,32) without such symmetry, and nothing about the 17 disjoint copies that a colouring would need.

Two side facts used or checked here:

* Every S(15,16,32) is closed under complements (Hoffman-bound equality, argument by KentaKitamura on the problem
  page); this is used to merge each orbit with its complementary orbit (`km835x.py`, `km835lp.py`).
* The block intersection numbers of a hypothetical S(15,16,32) are determined by the parameters and are
  nonnegative integers (n_0, …, n_14 = 1, 0, 960, 17920, 196560, 1118208, 3779776, 7687680, 9755460, 7687680,
  3779776, 1118208, 196560, 17920, 960; `intersect.py`), so this classical test gives no nonexistence proof. In the
  family S(t,t+1,2t+2) the divisibility conditions hold exactly when t+2 is prime, matching Ma–Tang.

## Validation

The same code finds the known systems and rejects the known non-existent ones:
S(3,4,8) under AGL(3,2) and under PSL(2,7) (the 14 planes), S(3,4,16) under AGL(4,2) (the 140 planes),
S(5,6,12) under PSL(2,11) (the Witt design, 132 blocks); nothing for S(2,3,8), S(1,2,4), S(2,3,6), for S(5,6,12)
under PGL(2,11) (not an automorphism group of the Witt design) and for S(7,8,16) under AGL(4,2). Orbit sizes are
computed from stabiliser counts and must sum exactly to the number of sets; every row of the matrix must sum to
v−k+1. `k4_check.py` verifies the colouring reduction at k = 4 directly: 30 systems S(3,4,8), all
complement-closed, at most 2 pairwise disjoint (χ(J(8,4)) = 6).

## Files

`km835.py` (affine groups AGL(n,2), orbits by marking with numba, Algorithm X), `km835g.py` (any permutation group
given by generators, SAT via python-sat), `km835x.py` (exact cover with complement merging and propagation, saves
the system as `.npz`), `km835lp.py`, `km835lp_ipm.py` (HiGHS LP relaxation and MIP), `intersect.py`,
`k4_check.py`; logs `km_n5_k16.log`, `kmg_*.log`, `kmx_*.log`, `run_*.out`, `lp_*.out`; systems `km_*.npz`
(the 92 MB AGL(1,32) file is not included).

Requirements: numpy, numba, scipy, python-sat. Reproduce the two proven negatives with
`python km835.py 5 16` (3 minutes) and `python km835x.py pgl2 31 16` (2 minutes).
