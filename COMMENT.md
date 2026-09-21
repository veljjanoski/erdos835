Draft comment for https://www.erdosproblems.com/835 (LaTeX, for Daniel's review; not posted)

---

For the smallest open case $k=16$ we looked for a single Steiner system $S(15,16,32)$ with a prescribed automorphism group $G$, by the Kramer–Mesner method: a $G$-invariant system is a union of $G$-orbits of $16$-sets, and the $G$-orbits of $15$-sets give the linear conditions.

Result: there is no $S(15,16,32)$ invariant under $\mathrm{AGL}(5,2)$ ($38$ orbits of $16$-sets, $30$ of $15$-sets) or under $\mathrm{PGL}(2,31)$ ($20636$ and $19234$ orbits; here some orbit of $15$-sets admits no orbit of $16$-sets at all). For $\mathrm{PSL}(2,31)$, $\mathrm{A\Gamma L}(1,32)$ and $\mathrm{AGL}(1,32)$ the systems have $21213$, $60714$ and $304174$ variables after identifying each orbit with its complementary orbit, and our solvers (SAT, exact cover, LP/MIP) did not decide them.

The code was checked on $S(3,4,8)$, $S(3,4,16)$ and $S(5,6,12)$ with their known groups. We also checked that the block intersection numbers of a hypothetical $S(15,16,32)$ are nonnegative integers, so this classical test gives no obstruction. None of this says anything about an $S(15,16,32)$ without such symmetry, or about the $17$ disjoint copies a colouring would need.

Code, logs and the undecided systems: https://github.com/veljjanoski/erdos835

AI-usage disclosure: Claude (Anthropic) was used as a coding assistant.
