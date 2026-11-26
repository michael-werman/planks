"""Figures for the paper "Planks through points" (drawn from exact constructions; the data of
figure 3 are floating-point best-known values from the paper's tables).
Usage: python3 make_figures.py   (writes ../paper1/figures/fig1_square_hexagon.png, fig2_sums.png, fig3_constants.png,
       fig4_hexagon_sector.png)
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations
from math import sqrt, pi, sin, cos
import os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper1", "figures")
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
INK, STRIP, PT, HI = "#2b2b2b", "#4c78a8", "#d62728", "#f58518"

def line_through(p, q):
    d = q - p; d = d/np.linalg.norm(d); nrm = np.array([-d[1], d[0]])
    return nrm, float(nrm @ p)          # {x : <nrm, x> = c}

def draw_panel(ax, P, eps, title):
    th = np.linspace(0, 2*pi, 400)
    ax.plot(np.cos(th), np.sin(th), color=INK, lw=1.2)
    X, Y = np.meshgrid(np.linspace(-1, 1, 700), np.linspace(-1, 1, 700))
    inside = X**2 + Y**2 <= 1
    lines = {}
    for p, q in combinations(range(len(P)), 2):
        nrm, c = line_through(P[p], P[q])
        if nrm[0] < -1e-12 or (abs(nrm[0]) < 1e-12 and nrm[1] < 0): nrm, c = -nrm, -c
        lines[(round(nrm[0], 9), round(nrm[1], 9), round(c, 9))] = (nrm, c)
    cover = np.zeros_like(X)
    for nrm, c in lines.values():
        cover += (np.abs(nrm[0]*X + nrm[1]*Y - c) <= eps)
    ax.contourf(X, Y, np.where(inside, cover, np.nan), levels=[0.5, 1.5, 2.5, 3.5, 20],
                colors=["#c6dbef", "#9ecae1", "#6baed6", "#4292c6"], alpha=0.55)
    for nrm, c in lines.values():
        d = np.array([nrm[1], -nrm[0]]); base = c*nrm; t = np.linspace(-2, 2, 2)
        seg = base[None] + t[:, None]*d[None]
        ax.plot(seg[:, 0], seg[:, 1], color=STRIP, lw=0.8)
    ax.scatter(P[:, 0], P[:, 1], s=22, color=PT, zorder=5)
    ax.scatter([0], [1], s=40, marker="*", color=HI, zorder=6)
    ax.annotate("(0, 1)", (0, 1), xytext=(0.12, 1.06), fontsize=9, color=HI)
    ax.set_xlim(-1.08, 1.08); ax.set_ylim(-1.08, 1.12); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10)

fig, axs = plt.subplots(1, 2, figsize=(8, 4.2))
a = 1/sqrt(10 + 4*sqrt(2))
draw_panel(axs[0], np.array([[a, a], [-a, a], [-a, -a], [a, -a]]), a,
           f"n = 4: centred square, ρ = a ≈ {a:.4f}")
r = 1/sqrt(14 + 4*sqrt(3))
H = np.array([[r*cos(k*pi/3), r*sin(k*pi/3)] for k in range(6)])
draw_panel(axs[1], H, r/2, f"n = 6: regular hexagon, ρ = r/2 ≈ {r/2:.4f}")
fig.tight_layout(); fig.savefig(os.path.join(out, "fig1_square_hexagon.png"), dpi=200); plt.close(fig)

# Figure 2: chord directions are sums
N, S = 6, [0, 1, 2, 4]
fig, axs = plt.subplots(1, 2, figsize=(8, 4.2))
ax = axs[0]; th = np.linspace(0, 2*pi, 400)
ax.plot(np.cos(th), np.sin(th), color="#bbbbbb", lw=1)
pts = {s: np.array([cos(2*pi*s/N), sin(2*pi*s/N)]) for s in S}
cols = plt.cm.tab10(np.arange(6))
for s, t in combinations(S, 2):
    k = (s + t) % N
    ax.plot(*zip(pts[s], pts[t]), color=cols[k], lw=2)
    mid = (pts[s] + pts[t])/2
    ax.annotate(f"{s}+{t}≡{k}", mid, fontsize=8, color=cols[k], ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
for s, p in pts.items():
    ax.scatter(*p, s=40, color=INK, zorder=5); ax.annotate(f"s={s}", p*1.2, ha="center", va="center", fontsize=9)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
ax.set_title("Points at angles 2πs/6, s ∈ {0,1,2,4}:\nchord {s,t} has direction π(s+t)/6 + π/2", fontsize=10)
ax = axs[1]
ax.plot(np.cos(th), np.sin(th), color=INK, lw=1.2)
for s, t in combinations(S, 2):
    k = (s + t) % N; ang = pi*k/N + pi/2
    d = np.array([cos(ang), sin(ang)])
    ax.plot([-d[0], d[0]], [-d[1], d[1]], color=cols[k], lw=2)
    ax.annotate(f"{k}", 1.2*d, color=cols[k], ha="center", va="center", fontsize=9)
gap_mid = pi/2 + pi/12; m = np.array([cos(gap_mid), sin(gap_mid)])
ax.scatter(*m, s=40, marker="*", color=HI, zorder=6)
ax.annotate("ρ = sin(π/12) here", m, xytext=(-1.45, 1.38), fontsize=9, color=HI,
            arrowprops=dict(arrowstyle="-", color=HI, lw=0.8))
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.4, 1.5)
ax.set_title("Shrunk to the centre: a pencil with all gaps π/6,\nso ρ → sin(π/12)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(os.path.join(out, "fig2_sums.png"), dpi=200, bbox_inches="tight", pad_inches=0.15); plt.close(fig)

# Figure 3: best known n^2 rho against the asymptotic constants (floating point values from the paper)
best = {3: 9*0.5, 4: 16/sqrt(10+4*sqrt(2)), 5: 4.261, 6: 36/sqrt(56+16*sqrt(3)),
        7: 4.051, 8: 4.044, 9: 4.045, 10: 4.185, 11: 4.095, 12: 4.274}
kind = {3: "triangle", 4: "square", 6: "hexagon"}
fig, ax = plt.subplots(figsize=(8, 4.2))
ns = sorted(best)
ax.plot(ns, [best[n] for n in ns], "o-", color=PT, lw=1.2, label="best known configuration (n²ρ)")
for n in ns:
    if n in kind: ax.annotate(kind[n], (n, best[n]), xytext=(n - 0.2, best[n] + (-0.17 if n == 6 else 0.08)), fontsize=8)
refs = [(9*pi/7, "9π/7 ≈ 4.04  upper constant (Thm B)", "#54a24b"),
        (pi, "π  central clusters, lower (Prop. 3)", "#b279a2"),
        (32/15, "32/15 ≈ 2.133  lower constant (Thm A)", STRIP),
        (2.0, "2  Bang's plank bound", "#888888")]
for y, lab, c in refs:
    ax.hlines(y, 2.7, 12.25, color=c, lw=1, ls="--"); ax.text(12.4, y, lab, va="center", fontsize=8, color=c)
ax.set_xlim(2.7, 17.2); ax.set_xticks(range(3, 13)); ax.set_ylim(1.6, 4.8); ax.set_xlabel("n"); ax.set_ylabel("n² · ρ")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(out, "fig3_constants.png"), dpi=200, bbox_inches="tight", pad_inches=0.15); plt.close(fig)

# Figure 4: the sector argument for the hexagon (Theorem D, n = 6)
r = 1/sqrt(14 + 4*sqrt(3)); h = r/2
fig, ax = plt.subplots(figsize=(6.2, 4.6))
X, Y = np.meshgrid(np.linspace(-0.05, 1.32, 1000), np.linspace(-0.12, 0.80, 700))
L1 = Y                                  # the long diagonal y = 0
L2 = Y - r*sqrt(3)/2                    # the side y = r*sqrt(3)/2
L3 = X/2 - sqrt(3)/2*Y - r/2            # the short diagonal x/2 - (sqrt3/2) y = r/2
for L, lab in ((L1, None), (L2, None), (L3, None)):
    ax.contourf(X, Y, (np.abs(L) <= h).astype(float), levels=[0.5, 1.5], colors=["#9ecae1"], alpha=0.55)
sector = (Y >= 0) & (Y <= X*np.tan(pi/6))
free = sector & (L1 > h) & (L2 > h) & (L3 > h)
ax.contourf(X, Y, free.astype(float), levels=[0.5, 1.5], colors=[HI], alpha=0.35)
th = np.linspace(-0.12, 0.80, 300); ax.plot(np.cos(th), np.sin(th), color=INK, lw=1.3)
ax.text(cos(0.62)+0.015, sin(0.62), "unit circle", fontsize=9, color=INK)
for ang in (0, pi/6): ax.plot([0, 1.3*cos(ang)], [0, 1.3*sin(ang)], color=INK, lw=0.8, ls="--")
ax.text(1.06, 0.585, "arg z = 30°", fontsize=9, rotation=30, color=INK)
xs = np.array([-0.05, 1.32])
ax.plot(xs, 0*xs, color=STRIP, lw=1); ax.plot(xs, 0*xs + r*sqrt(3)/2, color=STRIP, lw=1)
ax.plot(xs, (xs/2 - r/2)/(sqrt(3)/2), color=STRIP, lw=1)
ax.text(-0.04, 0.012, "long diagonal  y = 0", fontsize=8.5, color=STRIP)
ax.text(-0.04, r*sqrt(3)/2 + 0.012, "side  y = r√3/2", fontsize=8.5, color=STRIP)
ax.text(0.47, 0.035, "short diagonal\nx/2 − (√3/2)y = r/2", fontsize=8.5, color=STRIP)
V = np.array([[r*cos(k*pi/3), r*sin(k*pi/3)] for k in range(6)])
ax.scatter(V[:, 0], V[:, 1], s=18, color=PT, zorder=5)
cx, cy = (7 + sqrt(3))*r/2, (1 + sqrt(3))*r/2
ax.scatter([cx], [cy], s=40, marker="*", color=HI, zorder=6)
ax.annotate("corner ((7+√3)r/2, (1+√3)r/2),\non the unit circle", (cx, cy), xytext=(0.5, 0.6),
            fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))
ax.set_xlim(-0.05, 1.32); ax.set_ylim(-0.12, 0.80); ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout(); fig.savefig(os.path.join(out, "fig4_hexagon_sector.png"), dpi=200, bbox_inches="tight", pad_inches=0.1); plt.close(fig)
print("wrote", sorted(os.listdir(out)))
