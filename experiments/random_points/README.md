Floating-point experiments for Paper II (not part of any proof).

- `density_disc.py`: Monte Carlo estimate of P(dist(z, line U1U2) < t)/t in the unit disc at t = 0.002 from 4·10^7 pairs,
  at |z| = 0, 0.5, 0.9, 1, against h(z) = (2/(3π²))∫_0^π L_z(θ)³ dθ. Output (2026-11-25):
  0.0: 1.7013 vs 1.6977; 0.5: 1.3952 vs 1.3948; 0.9: 0.8502 vs 0.8485; 1.0: 0.7213 vs 0.7205;
  limit constant 8/h(boundary) = 9π²/8 = 11.1033.
- `rho_disc.py`: the covering radius of n uniform points in the unit disc (grid lower bound and Lipschitz upper
  bound), medians of n²ρ/log n over 20–40 trials per n (n = 8, 12, 20, 30, 40); output in `rho_disc.out`.
