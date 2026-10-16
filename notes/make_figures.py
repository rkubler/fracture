"""Figures for 'Stress intensity at the crack tip' lecture notes.
Every figure is computed from the same closed-form solutions used in the HTML pages."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10,
                     "axes.grid": True, "grid.alpha": .3, "savefig.bbox": "tight"})
BLUE, ORANGE, RED, GREY = "#1d58b8", "#d9561a", "#9b1230", "#566475"
CMAP = "jet"

def save(fig, name, png=False):
    fig.savefig(os.path.join(OUT, name + ".pdf"))
    if png: fig.savefig(os.path.join(OUT, name + ".png"), dpi=130)
    plt.close(fig)

def contour(ax, X, Y, F, vmin, vmax, n=13):
    """Filled rainbow contours (blue = min, red = max) with thin isolines."""
    lv = np.linspace(vmin, vmax, n)
    Fc = np.clip(F, vmin, vmax)
    cs = ax.contourf(X, Y, Fc, levels=lv, cmap=CMAP, extend="neither")
    ax.contour(X, Y, Fc, levels=lv[1:-1], colors="k", linewidths=.35, alpha=.55)
    ax.set_aspect("equal"); ax.grid(False)
    return cs

# ---------- elliptical hole, Kolosov-Muskhelishvili ----------
def ellipse_syy(x, y, a, b, s=1.0):
    b = min(b, 0.9998 * a)
    R = (a + b) / 2; m = (a - b) / (a + b); c = (1 + m) * (1 + m * m) / m
    z = x + 1j * y
    q = np.sqrt(z * z - 4 * R * R * m + 0j)
    z1, z2 = (z + q) / (2 * R), (z - q) / (2 * R)
    ze = np.where(np.abs(z1) >= np.abs(z2), z1, z2)
    om1 = R * (1 - m / ze**2); om2 = 2 * R * m / ze**3
    p1 = s * R / 4 * (1 + (2 + m) / ze**2); p2 = s * R / 4 * (-2 * (2 + m) / ze**3)
    d = ze**2 - m
    ps1 = -(s * R / 2) * (-1 + 1 / (m * ze**2) + c * (d - 2 * ze**2) / d**2)
    Phi = p1 / om1; dPhi = (p2 * om1 - p1 * om2) / om1**3; Psi = ps1 / om1
    S = 4 * Phi.real; D = 2 * (np.conj(z) * dPhi + Psi)
    syy = (S + D.real) / 2
    syy[(x / a)**2 + (y / b)**2 <= 1] = np.nan
    return syy

X, Y = np.meshgrid(np.linspace(-2.4, 2.4, 481), np.linspace(-1.5, 1.5, 301))
fig, axs = plt.subplots(1, 3, figsize=(10, 2.6), constrained_layout=True)
for ax, ba in zip(axs, [1.0, 0.3, 0.05]):
    f = ellipse_syy(X, Y, 1.0, ba)
    im = contour(ax, X, Y, f, -1, 5)
    t = np.linspace(0, 2 * np.pi, 200); ax.fill(np.cos(t), ba * np.sin(t), "w"); ax.plot(np.cos(t), ba * np.sin(t), "k", lw=.8)
    ax.set_title(rf"$b/a={ba}$,  $K_t={1+2/ba:.0f}$"); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
fig.colorbar(im, ax=axs, shrink=.9, label=r"$\sigma_{yy}/\sigma_\infty$ (clipped to $[-1,5]$)")
save(fig, "ellipse_fields", png=True)

aro = np.logspace(0, 4, 200)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.2), constrained_layout=True)
a1.loglog(aro, 1 + 2 * np.sqrt(aro), color=ORANGE, lw=2)
a1.set_xlabel(r"$a/\rho$"); a1.set_ylabel(r"$K_t$"); a1.set_title(r"$K_t = 1+2\sqrt{a/\rho}$ grows without limit")
a2.semilogx(aro, 1 + 0.5 / np.sqrt(aro), color=BLUE, lw=2); a2.axhline(1, ls="--", color=GREY)
a2.set_xlabel(r"$a/\rho$"); a2.set_ylabel(r"$\frac{\sqrt{\pi}}{2}\sigma_{\max}\sqrt{\rho}\,/\,\sigma_\infty\sqrt{\pi a}$")
a2.set_title(r"tends to 1: the limit is $K_I=\sigma_\infty\sqrt{\pi a}$"); a2.set_ylim(.95, 1.55)
save(fig, "kt_limit")

# ---------- Westergaard centre crack ----------
def westergaard_syy(x, y, a=1.0):
    z = x + 1j * y
    sq = np.sqrt(z - a + 0j) * np.sqrt(z + a + 0j)
    Z = z / sq; Zp = -a * a / ((z * z - a * a) * sq)
    return Z.real + y * Zp.imag
Yw = np.linspace(-1.5, 1.5, 300) + 0.005
X2, Y2 = np.meshgrid(np.linspace(-2.4, 2.4, 481), Yw)
def westergaard_all(x, y, a=1.0):
    z = x + 1j * y
    sq = np.sqrt(z - a + 0j) * np.sqrt(z + a + 0j)
    Z = z / sq; Zp = -a * a / ((z * z - a * a) * sq)
    return Z.real - y * Zp.imag - 1.0, Z.real + y * Zp.imag, -y * Zp.real   # sxx (uniaxial), syy, txy
sxx, syy, txy = westergaard_all(X2, Y2)
fig, axs = plt.subplots(1, 3, figsize=(10.5, 2.7), constrained_layout=True)
for ax, F, lab, rng in zip(axs, [syy, sxx, txy], [r"\sigma_{yy}", r"\sigma_{xx}", r"\tau_{xy}"], [(-1, 5), (-2, 4), (-2, 2)]):
    cs = contour(ax, X2, Y2, F, *rng); ax.plot([-1, 1], [0, 0], "w", lw=2.5); ax.plot([-1, 1], [0, 0], "k", lw=1)
    ax.set_title(rf"${lab}/\sigma_\infty$"); ax.set_xticks([-2, -1, 0, 1, 2]); ax.set_yticks([-1, 0, 1]); ax.set_xlabel(r"$x/a$")
    fig.colorbar(cs, ax=ax, shrink=.8)
axs[0].set_ylabel(r"$y/a$")
save(fig, "westergaard_fields", png=True)
fig, a2 = plt.subplots(figsize=(5.5, 3.2), constrained_layout=True)
xs = 1 + np.logspace(-3, np.log10(1.4), 300)
a2.plot(xs - 1, xs / np.sqrt(xs**2 - 1), color=ORANGE, lw=2, label="exact")
r = xs - 1; a2.plot(r, 1 / np.sqrt(2 * r), "--", color=BLUE, lw=1.5, label=r"$K/\sqrt{2\pi r}$ term")
a2.axhline(1, ls=":", color=GREY); a2.set_ylim(0, 8); a2.set_xlim(0, 1.4)
a2.set_xlabel(r"$(x-a)/a = r/a$"); a2.set_ylabel(r"$\sigma_{yy}/\sigma_\infty$ on $y=0$"); a2.legend()
save(fig, "westergaard")

# ---------- Williams near-tip field ----------
def fI(t):
    c, s = np.cos(t / 2), np.sin(t / 2)
    return dict(xx=c * (1 - s * np.sin(1.5 * t)), yy=c * (1 + s * np.sin(1.5 * t)), xy=c * s * np.cos(1.5 * t))
X3, Y3 = np.meshgrid(np.linspace(-1.6, 1.6, 481), np.linspace(-1, 1, 301) + 0.002)
Rr, T = np.hypot(X3, Y3), np.arctan2(Y3, X3)
fig, axs = plt.subplots(1, 3, figsize=(10.5, 2.6), constrained_layout=True)
for ax, k, lab, rng in zip(axs, ["yy", "xx", "xy"], [r"\sigma_{yy}", r"\sigma_{xx}", r"\tau_{xy}"], [(0, 3), (0, 3), (-1.5, 1.5)]):
    v = fI(T)[k] / np.sqrt(Rr)
    cs = contour(ax, X3, Y3, v, *rng)
    ax.plot([-1.6, 0], [0, 0], "w", lw=2.5); ax.plot([-1.6, 0], [0, 0], "k", lw=1)
    ax.set_title(rf"${lab}\,\sqrt{{2\pi\,r_0}}/K_I$,  $r_0=1$ mm"); ax.set_xlabel(r"$x$ (mm)"); ax.set_yticks([-1, 0, 1])
    fig.colorbar(cs, ax=ax, shrink=.8)
axs[0].set_ylabel(r"$y$ (mm)")
save(fig, "williams_fields", png=True)

K = 25.0; rr = np.logspace(-2, np.log10(1.6), 100)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3), constrained_layout=True)
a1.loglog(rr, K / np.sqrt(2 * np.pi * rr * 1e-3), color=ORANGE, lw=2)
a1.set_xlabel(r"$r$ (mm)"); a1.set_ylabel(r"$\sigma_{yy}(r,0)$ (MPa)"); a1.set_title(r"slope $-1/2$ for $K_I=25$ MPa$\sqrt{\rm m}$")
a2.semilogx(rr, np.full_like(rr, K), color=BLUE, lw=2); a2.set_ylim(0, 40)
a2.set_xlabel(r"$r$ (mm)"); a2.set_ylabel(r"$\sigma_{yy}\sqrt{2\pi r}$ (MPa$\sqrt{\rm m}$)"); a2.set_title(r"constant $=K_I$ at every $r$")
save(fig, "singularity_slope")

# ---------- geometry factors ----------
Yc = lambda r: 1 / np.sqrt(np.cos(np.pi * r))
Ye = lambda r: 1.12 - 0.231 * r + 10.55 * r**2 - 21.72 * r**3 + 30.39 * r**4
fig, ax = plt.subplots(figsize=(6, 3.4))
r1 = np.linspace(0, .35, 100); r2 = np.linspace(0, .6, 100)
ax.plot(r1, Yc(r1), color=BLUE, lw=2, label=r"centre crack, $\sqrt{\sec(\pi a/W)}$")
ax.plot(r2, Ye(r2), color=ORANGE, lw=2, label="single edge crack (Tada)")
ax.axhline(1, ls="--", color=GREY, label="infinite plate, $Y=1$"); ax.axhline(1.122, ls=":", color=RED, label="semi-infinite edge, $Y=1.122$")
ax.plot(.1, Yc(.1), "ko"); ax.annotate(r"example: $a/W=0.1$, $Y=1.025$", (.1, Yc(.1)), (.15, 1.6), arrowprops=dict(arrowstyle="->"))
ax.set_xlabel(r"$a/W$"); ax.set_ylabel(r"$Y$"); ax.set_ylim(.8, 4.2); ax.legend(fontsize=8)
save(fig, "geometry_factors")

# ---------- extrapolation ----------
sigma, a, = 120.0, 20.0
Kref = sigma * np.sqrt(np.pi * a * 1e-3)
n = 18; rs = 0.02 * (2.0 / 0.02) ** (np.arange(n) / (n - 1))
noise = (np.sin(np.arange(n) * 12.9898) * 43758.5453) % 1
noise = np.fmod(np.sin(np.arange(n) * 12.9898) * 43758.5453, 1)
data = {"str": Kref * (1 + rs / a) / np.sqrt(1 + rs / (2 * a)) * (1 - 0.32 * np.exp(-rs / 0.07)) * (1 + 0.004 * noise),
        "dis": Kref * np.sqrt(1 - rs / (2 * a)) * (1 - 0.07 * np.exp(-rs / 0.07)) * (1 + 0.004 * noise)}
res = {}
fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), constrained_layout=True)
for ax, k, title in zip(axs, ["str", "dis"], ["stress method, nodes at $\\theta=0$", "displacement method, flank nodes $\\theta=\\pm\\pi$"]):
    d = data[k]; ax.axhline(Kref, ls=":", color=GREY, label=rf"$\sigma\sqrt{{\pi a}}={Kref:.2f}$")
    for rmin, col, ls in [(0.3, BLUE, "-"), (0.0, ORANGE, "--")]:
        sel = rs >= rmin; p = np.polyfit(rs[sel], d[sel], 1); res[(k, rmin)] = p[1]
        xx = np.linspace(0, 2.1, 10); ax.plot(xx, np.polyval(p, xx), ls, color=col, lw=1.5, label=rf"fit $r\geq{rmin}$ mm: $K={p[1]:.2f}$")
    ax.plot(rs[rs >= .3], d[rs >= .3], "o", color=BLUE, ms=4); ax.plot(rs[rs < .3], d[rs < .3], "o", mfc="none", color=ORANGE, ms=4)
    ax.set_xlabel(r"$r$ (mm)"); ax.set_ylabel(r"$K_{\rm app}$ (MPa$\sqrt{\rm m}$)"); ax.set_title(title); ax.legend(fontsize=7.5)
save(fig, "extrapolation")

# ---------- fracture: K vs a ----------
W, sig = 200.0, 150.0
mats = [("Al 7075-T6", 29, 505), ("AISI 4340", 50, 1500), ("Ti-6Al-4V", 66, 880)]
aa = np.linspace(0.1, 0.35 * W, 400); Ka = Yc(aa / W) * sig * np.sqrt(np.pi * aa * 1e-3)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(aa, Ka, color=BLUE, lw=2, label=r"$K_I=\sqrt{\sec(\pi a/W)}\,\sigma\sqrt{\pi a}$")
acs = {}
for (nme, kic, sy), col in zip(mats, [ORANGE, RED, GREY]):
    ax.axhline(kic, ls="--", color=col, lw=1.2)
    Kf = lambda x: Yc(x / W) * sig * np.sqrt(np.pi * x * 1e-3)
    ac = None
    if Kf(0.35 * W) >= kic:
        lo, hi = 1e-3, 0.35 * W
        for _ in range(80):
            mid = (lo + hi) / 2; lo, hi = (lo, mid) if Kf(mid) >= kic else (mid, hi)
        ac = hi
    acs[nme] = ac
    if ac: ax.plot(ac, kic, "o", color=col); ax.annotate(rf"{nme}: $a_c={ac:.1f}$ mm", (ac, kic), (ac + 2, kic + 4), color=col, fontsize=8)
ax.axvline(20, ls=":", color="k"); ax.text(21, 5, "current $a=20$ mm", fontsize=8)
ax.set_xlabel(r"half crack length $a$ (mm)"); ax.set_ylabel(r"$K_I$ (MPa$\sqrt{\rm m}$)"); ax.set_ylim(0, 110); ax.legend(fontsize=8, loc="upper left")
ax.set_title(r"centre crack, $W=200$ mm, $\sigma=150$ MPa")
save(fig, "fracture_curve")

# ---------- angular functions ----------
t = np.linspace(-np.pi, np.pi, 721); c, s = np.cos, np.sin; h = t / 2
F = {"I": {"cart": [(r"$f_{xx}$", c(h) * (1 - s(h) * s(3 * h))), (r"$f_{yy}$", c(h) * (1 + s(h) * s(3 * h))), (r"$f_{xy}$", c(h) * s(h) * c(3 * h))],
           "polar": [(r"$f_{rr}$", 1.25 * c(h) - .25 * c(3 * h)), (r"$f_{\theta\theta}$", .75 * c(h) + .25 * c(3 * h)), (r"$f_{r\theta}$", .25 * s(h) + .25 * s(3 * h))]},
     "II": {"cart": [(r"$f_{xx}$", -s(h) * (2 + c(h) * c(3 * h))), (r"$f_{yy}$", s(h) * c(h) * c(3 * h)), (r"$f_{xy}$", c(h) * (1 - s(h) * s(3 * h)))],
            "polar": [(r"$f_{rr}$", -1.25 * s(h) + .75 * s(3 * h)), (r"$f_{\theta\theta}$", -.75 * s(h) - .75 * s(3 * h)), (r"$f_{r\theta}$", .25 * c(h) + .75 * c(3 * h))]},
     "III": {"cart": [(r"$f_{xz}$", -s(h)), (r"$f_{yz}$", c(h))], "polar": [(r"$f_{rz}$", s(h)), (r"$f_{\theta z}$", c(h))]}}
fig, axs = plt.subplots(3, 2, figsize=(9, 7.5), sharex=True, constrained_layout=True)
for i, mo in enumerate(["I", "II", "III"]):
    for j, sy in enumerate(["cart", "polar"]):
        ax = axs[i, j]
        for (lab, v), col in zip(F[mo][sy], [BLUE, ORANGE, RED]): ax.plot(np.degrees(t), v, color=col, lw=1.8, label=lab)
        ax.axhline(0, color="k", lw=.6); ax.set_ylim(-2.6, 2.1); ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_title(f"Mode {mo}, " + ("Cartesian" if sy == "cart" else "polar") + " components", fontsize=10); ax.legend(fontsize=8, ncol=3, loc="lower center")
for ax in axs[-1]: ax.set_xlabel(r"$\theta$ (degrees)")
save(fig, "angular_functions")

with open(os.path.join(OUT, "numbers.txt"), "w") as fh:
    for k, v in res.items(): fh.write(f"{k}: {v:.3f}\n")
    for k, v in acs.items(): fh.write(f"ac {k}: {v}\n")
    fh.write(f"Kref {Kref:.3f}\n")
print(open(os.path.join(OUT, "numbers.txt")).read())

# ---------- small-scale yielding: Irwin stress profile ----------
ry = 1.0; r = np.linspace(0.02, 6, 600)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(r, np.sqrt(ry / r), "--", color=GREY, lw=1.5, label=r"elastic $K/\sqrt{2\pi r}$")
plastic = np.where(r <= 2 * ry, 1.0, np.sqrt(ry / np.maximum(r - ry, 1e-9)))
ax.plot(r, plastic, color=ORANGE, lw=2.2, label="Irwin: yielded zone + shifted elastic field")
ax.axhline(1, color="k", lw=.6); ax.axvline(ry, ls=":", color=BLUE); ax.axvline(2 * ry, ls=":", color=RED)
ax.fill_between(r, plastic, np.minimum(np.sqrt(ry / r), 3), where=r <= ry, color=BLUE, alpha=.15, label="load shed by yielding")
ax.fill_between(r, 1, np.sqrt(ry / r), where=(r >= ry) & (r <= 2 * ry), color=ORANGE, alpha=.12)
ax.text(ry + .05, 2.6, r"$r_y$", color=BLUE); ax.text(2 * ry + .05, 2.6, r"$r_p=2r_y$", color=RED)
ax.set_xlim(0, 6); ax.set_ylim(0, 3); ax.set_xlabel(r"$r/r_y$"); ax.set_ylabel(r"$\sigma_{yy}/\sigma_{ys}$ on $\theta=0$")
ax.legend(fontsize=8, loc="upper right")
save(fig, "irwin_profile")
