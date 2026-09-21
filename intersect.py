"""Block intersection numbers of a hypothetical Steiner system S(t,t+1,v), v=2t+2.
n_j = number of blocks (other than a fixed block B) meeting B in exactly j points, j=0..t-1 (n_t=0, two blocks share
<= t-1 points).  Counting blocks through each i-subset of B (i=0..t-1) gives sum_j C(j,i) n_j = C(k,i)(lambda_i - 1),
lambda_i = C(v-i, t-i)/C(k-i, t-i).  Square triangular system -> unique rational solution; must be nonnegative integers."""
from fractions import Fraction as F
from math import comb
def numbers(t):
    k = t + 1; v = 2 * t + 2
    lam = [F(comb(v - i, t - i), comb(k - i, t - i)) for i in range(t + 1)]
    n = [F(0)] * k                      # n_0..n_{k-1}; n_t = 0
    for i in range(t - 1, -1, -1):      # solve from the top: equation i involves n_i..n_{t-1}
        s = sum(comb(j, i) * n[j] for j in range(i + 1, t))
        n[i] = (comb(k, i) * (lam[i] - 1) - s) / comb(i, i)
    return lam, n[:t]
for t in range(1, 22):
    lam, n = numbers(t)
    ok_lam = all(x.denominator == 1 for x in lam)
    ok_n = all(x.denominator == 1 and x >= 0 for x in n)
    flag = "OK" if ok_lam and ok_n else ("lambda non-integer" if not ok_lam else "INTERSECTION NUMBERS FAIL")
    print(f"t={t:2d} S({t},{t+1},{2*t+2}): {flag}", end="")
    if t <= 7 or not (ok_lam and ok_n): print("  n =", [str(x) for x in n])
    else: print()
