# Restricted 2-bases of cyclic groups: p(k) ≥ (7/18 − o(1))k²

Status: referee ACCEPT WITH FIXES (notes/reviews/sigma2_seven_eighteenths.md); fixes applied 2026-11-24. Written to be read on its own by a
specialist in additive bases. It is the case d = 2 of the seed lemma in proofs/blowup_lim_sup.md
(refereed), specialised to one seed, with every constant explicit and a direct proof that does not
need the general d. Exact checks: certificates/sigma2_seven_eighteenths/ (see RECORD.txt there).

## Definitions

- For a finite set A in an abelian group, A ∔ A = {a + b : a, b ∈ A, a ≠ b}.
- For k ≥ 3, p(k) is the largest N ≥ 1 such that some k-subset A ⊂ ℤ_N has A ∔ A = ℤ_N.
  This is Graham and Sloane's n_γ(k) (SIAM J. Algebraic Discrete Methods 1 (1980), 382–404),
  Haanpää's n′_ss(k) (J. Integer Seq. 7 (2004), Art. 04.2.6) and OEIS A004131. Since every
  element of ℤ_N must be one of at most C(k,2) sums, p(k) ≤ C(k,2).
- Q = lcm(1, 2, …, 11) = 27720 = 2³·3²·5·7·11. Note 56 = 2³·7 divides Q.

## Theorem

For every k ≥ 360Q = 9 979 200,

```
p(k) ≥ (7/18)·(k − 288Q)²,        hence   liminf_{k→∞} p(k)/k² ≥ 7/18.
```

Context (from the sources above, not proved here): Graham and Sloane prove (5/18)(k−1)² < n_γ(k)
(their Theorem 1(2)). For the unrestricted problem (sums a + b with a = b allowed), Jia and Shen
(SIAM J. Discrete Math. 31 (2017), 796–804) prove (1/3 − o(1))k². A set with A ∔ A = ℤ_N also has
A + A = ℤ_N, so the theorem gives 7/18 for the unrestricted problem too. We found no published
asymptotic bound above 1/3 for either problem (notes/novelty_check.md); a specialist should
confirm this.

## The construction in one line

In ℤ_56 × ℤ_m × ℤ_{m+Q} ≅ ℤ_{56m(m+Q)} (gcd(m, Q) = 1), take
S = ⋃_{t=1}^{12} {l_t} × {i·(1, t) : 0 ≤ i < m + 23Q}, where L = {l₁, …, l₁₂} is the seed below.

## Proof

**Lemma 1 (the seed).** L = {0, 1, 12, 15, 22, 29, 32, 43, 44, 48, 50, 52} ⊂ ℤ_56 has L ∔ L = ℤ_56.

*Proof.* A finite check; certificates/sigma2_seven_eighteenths/sigma2_check.py prints one representation r = l_s + l_t (s < t) for each
r ∈ ℤ_56. This set is the example for k = 12 in Haanpää's paper, where n′_ss(12) = 56. ∎

Number the elements l₁, …, l₁₂ as listed. For each r ∈ ℤ_56 fix indices s_r < t_r with
l_{s_r} + l_{t_r} = r.

**Lemma 2 (two progressions cover).** Let m ≥ 1 with gcd(m, Q) = 1. Put a = m, b = m + Q,
M = ab and T = m + 23Q. Then gcd(a, b) = 1, so ℤ_M ≅ ℤ_a × ℤ_b. For t ∈ {1, …, 12} let g_t ∈ ℤ_M
be the element with g_t ≡ 1 (mod a) and g_t ≡ t (mod b), and X_t = {i·g_t : i = 0, 1, …, T − 1}.
Then X_s + X_t = ℤ_M for all 1 ≤ s < t ≤ 12.

*Proof.* gcd(b, Q) = gcd(m, Q) = 1 and gcd(a, b) = gcd(m, Q) = 1. Fix s < t and put
δ = t − s ∈ {1, …, 11}, so δ | Q and gcd(δ, b) = 1.

- *The map.* Let φ : ℤ² → ℤ_M, φ(x, y) = x·g_s + y·g_t. In ℤ_a × ℤ_b it reads
  φ(x, y) = (x + y mod a, sx + ty mod b).
- *φ is onto.* Given (u, v) ∈ ℤ_a × ℤ_b, pick y ∈ ℤ with δy ≡ v − su (mod b) (δ is invertible
  mod b) and put x = u − y. Then x + y = u and sx + ty = su + δy ≡ v. So the kernel K of φ is a
  subgroup of ℤ² of index M.
- *A basis of K.* Let

  ```
  G = [ m − Qs/δ     −Qt/δ    ]
      [   Qs/δ      m + Qt/δ  ]
  ```

  Its entries are integers since δ | Q. For the first column (x, y): x + y = m ≡ 0 (mod a), and
  sx + ty = sm − Qs²/δ + Qst/δ = sm + Qs = s(m + Q) ≡ 0 (mod b). For the second column:
  x + y = m, and sx + ty = −Qst/δ + tm + Qt²/δ = t(m + Q). So Gℤ² ⊂ K. Moreover
  det G = (m − Qs/δ)(m + Qt/δ) + Q²st/δ² = m² + mQ(t − s)/δ = m(m + Q) = M, so Gℤ² has index M
  as well, and K = Gℤ².
- *A box meets every coset of K.* Every coset of K in ℤ² contains a point z = Gw with
  w ∈ [0, 1)² (write any point of the coset as Gw′ with w′ ∈ ℚ² and subtract G⌊w′⌋). Row j of G
  (j = 1, 2) has entries γ_{j1}, γ_{j2}, not both zero, and z_j = γ_{j1}w₁ + γ_{j2}w₂ lies in the
  interval [α_j, β_j] with α_j = Σ_i min(γ_{ji}, 0) and β_j = Σ_i max(γ_{ji}, 0), integers.
  - If row j has a negative entry, then α_j is not attained (that entry times w_i < 1 exceeds the
    entry), so z_j ∈ {α_j + 1, …, β_j}; put c_j = α_j + 1.
  - Otherwise α_j = 0 and some entry is positive, so β_j is not attained and
    z_j ∈ {α_j, …, β_j − 1}; put c_j = α_j.
  - Either way z_j ∈ c_j + {0, …, β_j − α_j − 1}, and β_j − α_j = |γ_{j1}| + |γ_{j2}|. Both row sums
    are at most m + Q(s + t)/δ ≤ m + 23Q = T (|m − Qs/δ| ≤ m + Qs/δ, s + t ≤ 23, δ ≥ 1).
  - So every coset of K meets c + {0, …, T−1}², where c = (c₁, c₂) depends only on G.
- *Conclusion.* Hence φ(c + {0, …, T−1}²) = ℤ_M, and since φ is a homomorphism,
  φ({0, …, T−1}²) = ℤ_M − φ(c) = ℤ_M. That is X_s + X_t = ℤ_M. ∎

**Lemma 3 (the lift).** With m, M, T, X_t as in Lemma 2, N = 56M and
S = {(l_t, x) : 1 ≤ t ≤ 12, x ∈ X_t} ⊂ ℤ_56 × ℤ_M ≅ ℤ_N, we have |S| ≤ 12T and S ∔ S = ℤ_N.

*Proof.* gcd(M, 56) = 1 because 56 | Q and gcd(M, Q) = 1, so ℤ_56 × ℤ_M ≅ ℤ_N. Clearly
|S| ≤ 12T. Let (r, h) ∈ ℤ_56 × ℤ_M and s = s_r < t = t_r. By Lemma 2, h = x + y with x ∈ X_s,
y ∈ X_t. Then (r, h) = (l_s, x) + (l_t, y), and the two summands are distinct elements of S because
l_s ≠ l_t. ∎

**Proof of the theorem.**

- *Monotonicity.* If A ⊂ ℤ_N has A ∔ A = ℤ_N and |A| ≤ k ≤ N, any k-set A′ ⊃ A in ℤ_N has
  A′ ∔ A′ = ℤ_N. So by Lemma 3, p(k) ≥ 56m(m + Q) whenever gcd(m, Q) = 1 and
  12(m + 23Q) ≤ k ≤ 56m(m + Q).
- *Choice of m.* Let k ≥ 360Q and let m be the largest integer with m ≤ k/12 − 23Q and
  gcd(m, Q) = 1. The Q consecutive integers ⌊k/12 − 23Q⌋ − Q + 1, …, ⌊k/12 − 23Q⌋ contain one that
  is ≡ 1 (mod Q), hence coprime to Q, so m ≥ ⌊k/12 − 23Q⌋ − Q + 1 > k/12 − 24Q ≥ 6Q.
  Then 12(m + 23Q) ≤ k, and k < 12(m + 24Q) ≤ 12(m + 4m) = 60m < 56m(m + Q) since m + Q ≥ 2.
- *The bound.* p(k) ≥ 56m(m + Q) ≥ 56m² ≥ 56(k/12 − 24Q)² = (56/144)(k − 288Q)²
  = (7/18)(k − 288Q)². Dividing by k² and letting k → ∞ gives liminf p(k)/k² ≥ 7/18. ∎

## Remarks

- **Why 12.** The same proof works for any seed: a k₀-subset of ℤ_P with restricted 2-sum cover
  gives liminf p(k)/k² ≥ P/k₀² (proofs/blowup_lim_sup.md, part (b), any d). So
  lim p(k)/k² exists and equals sup_k p(k)/k² (part (c) there). Among the exact values in
  Haanpää's table (k = 3, …, 14: 3, 6, 9, 13, 17, 24, 30, 36, 42, 56, 61, 72) the largest ratio
  is 56/144 = 7/18 at k = 12; k = 4 gives 6/16 = 3/8 < 7/18, and k = 13, 14 give 61/169 ≈ 0.361
  and 72/196 ≈ 0.367.
- **What a specialist should check.** (i) That Lemma 2 is right: two arithmetic progressions of
  length ≈ √M, with steps (1, s) and (1, t) in ℤ_m × ℤ_{m+Q}, cover ℤ_M; this is the only
  non-routine step. (ii) That no published bound above 1/3 exists, for the restricted or the
  unrestricted cyclic problem.
- **The constants are not optimised.** Only the ratio 56/144 matters asymptotically; T could be
  reduced by choosing the t's and Q more carefully.
- **Checks.** certificates/sigma2_seven_eighteenths/sigma2_check.py verifies, in exact integer
  and rational arithmetic: Lemma 1; Q = 27720 with 56 | Q; for all 66 pairs s < t the integrality,
  kernel identities, det G = m(m + Q) and the row-sum bound of Lemma 2, as polynomial identities in
  m; by brute force, Lemma 2 for small analogues (4 or 5 indices, Q′ = 6, 12, M up to 165 613) and
  Lemma 3 for a lifted toy seed (Q′ = 36, N up to 1 577 133); and Lemma 2 at the real parameters
  (Q = 27720, m = 6Q + 1 and m = 1 000 003, M up to about 10¹²) for all 66 pairs, by the arc test:
  X_s + X_t = ℤ_M iff the T arcs [i·u, i·u + T) mod M, u = g_s g_t⁻¹, 0 ≤ i < T, cover ℤ_M,
  i.e. the sorted start points have all circular gaps ≤ T.
