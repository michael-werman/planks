Floating-point and high-precision checks for Paper II (not part of any proof).

- `h_polygon.py`: h(z) = (2/(3|K|^2)) ∫_0^π L_z(θ)^3 dθ for the unit square and the equilateral triangle of side 1, minimized
  over a 41x41 grid (floating point). Output (2026-11-27): square min h ≈ 1.18647 at the edge midpoint (0, 0.5), c ≈ 6.7427;
  triangle min h ≈ 1.73119 at the edge midpoint (0.5, 0), c ≈ 4.6211.
- `square_exact.py`: mpmath (30 digits) quadrature check of Corollary 2.6 (unit square): A(1/2) = √5/4 + (11/16) log φ,
  m = √5/3 + (11/12) log φ ≈ 1.186467, c = 96/(4√5 + 11 log φ) ≈ 6.742709, h along an edge, and A'(a) = 3 g(1/a).
- `h_triangle_scan.py`: floating-point scan (boundary step 1/200 per edge, 40x40 interior grid) of h for two triangles
  (2026-11-27): equilateral of side 1: min h ≈ 1.73119 at an edge midpoint, c ≈ 4.6211, interior grid min ≈ 2.387;
  vertices (0,0),(1,0),(0.05,0.3): min h ≈ 3.13481 on the long edge at x ≈ 0.335 (not its midpoint), interior grid min ≈ 4.80.
  The scan uses K=600 midpoint nodes in θ, then 20000 at the best grid point.
- `h_triangle_refine.py` (2026-11-27, floating point, ~11 s): golden-section search along the two candidate edges of the
  triangle (0,0),(1,0),(0.05,0.3) with K θ-nodes. K=2·10^5 and 10^6: bottom-edge minimum h ≈ 3.134735 at (0.3329, 0),
  c ≈ 2.55205; a second local minimum on the edge from (1,0) to (0.05,0.3), h ≈ 3.13618 at (0.364, 0.201), larger by
  ≈ 0.046%. With K=2000 the order of the two flips, so the location (not the value of c) is resolution-sensitive.
