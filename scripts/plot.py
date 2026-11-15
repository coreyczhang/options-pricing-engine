"""Render the figures in docs/img/ from the CSVs written by ope_export.

    ./build/ope_export build/plots
    python3 scripts/plot.py build/plots docs/img

Needs numpy and matplotlib. All numbers come from the C++ engine; this script
only draws them.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Chart chrome (light surface) and palette. Expiries are ordered, so they get
# steps of one sequential blue ramp (light = short-dated, dark = long-dated);
# the single-series validation panels use categorical slot 1 (blue) for the
# measurement and slot 2 (orange) for the reference it is compared against.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
RAMP_4 = ["#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"]
RAMP_5 = ["#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "axes.titlesize": 11,
    "axes.titleweight": "semibold",
    "axes.labelsize": 9,
    "axes.grid": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "legend.labelcolor": INK_2,
    "lines.linewidth": 1.5,
    "font.size": 9,
})


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def expiry_label(T):
    weeks = T * 52
    if abs(weeks - round(weeks)) < 1e-3 and weeks < 4:
        return f"{round(weeks)}w"
    months = T * 12
    if abs(months - round(months)) < 1e-3 and months < 12:
        return f"{round(months)}m"
    return f"{T:g}y"


def suptitle(fig, title, subtitle):
    fig.suptitle(title, x=0.01, ha="left", y=0.995, fontsize=13, fontweight="semibold", color=INK)
    fig.text(0.01, 0.945, subtitle, ha="left", va="top", fontsize=9, color=INK_2,
             linespacing=1.4)


# ---- Analytic Greeks across spot ------------------------------------------------

def greeks_vs_spot(data_dir, out_dir):
    rows = read(data_dir / "greeks_vs_spot.csv")
    by_T = defaultdict(list)
    for r in rows:
        by_T[float(r["expiry"])].append(r)
    expiries = sorted(by_T)

    # (column, title, scale to a trader's unit)
    panels = [("price", "Price", 1.0),
              ("delta", "Delta", 1.0),
              ("gamma", "Gamma", 1.0),
              ("vega", "Vega (per vol point)", 0.01),
              ("theta", "Theta (per calendar day)", 1 / 365),
              ("rho", "Rho (per 1% rate)", 0.01)]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), constrained_layout=True)
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.915))
    for ax, (col, title, scale) in zip(axes.flat, panels):
        for T, color in zip(expiries, RAMP_4):
            s = np.array([float(r["spot"]) for r in by_T[T]])
            y = np.array([float(r[col]) for r in by_T[T]]) * scale
            ax.plot(s, y, color=color, label=expiry_label(T))
        ax.axvline(100.0, color=AXIS, linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
        ax.set_title(title, loc="left")
        ax.set_xlabel("spot")
    # Direct labels where the curves separate: gamma peaks, tallest = shortest.
    ax = axes.flat[2]
    for T, color in zip(expiries, RAMP_4):
        s = np.array([float(r["spot"]) for r in by_T[T]])
        g = np.array([float(r["gamma"]) for r in by_T[T]])
        i = int(np.argmax(g))
        ax.annotate(expiry_label(T), (s[i], g[i]), xytext=(6, 0), textcoords="offset points",
                    fontsize=8, color=INK_2, va="center")
    axes.flat[0].legend(title="time to expiry", title_fontsize=8, loc="upper left")
    suptitle(fig, "Black-Scholes Greeks of a K = 100 call across spot",
             "r = 5%, σ = 20%, q = 0. Dashed line: the strike. Gamma and vega concentrate at "
             "the money and sharpen as expiry nears.")
    fig.savefig(out_dir / "greeks_vs_spot.png", dpi=150)
    plt.close(fig)


# ---- Monte Carlo Greeks against the closed form ----------------------------------

def mc_vs_analytic(data_dir, out_dir):
    rows = read(data_dir / "mc_greeks.csv")
    by_g = defaultdict(list)
    for r in rows:
        by_g[r["greek"]].append(r)
    paths = int(rows[0]["paths"])

    names = [("price", "Price"), ("delta", "Delta"), ("gamma", "Gamma"),
             ("vega", "Vega"), ("theta", "Theta"), ("rho", "Rho")]
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), constrained_layout=True)
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.89))
    worst_z = 0.0
    for ax, (key, title) in zip(axes.flat, names):
        s = np.array([float(r["spot"]) for r in by_g[key]])
        mc = np.array([float(r["mc"]) for r in by_g[key]])
        se = np.array([float(r["mc_se"]) for r in by_g[key]])
        exact = np.array([float(r["analytic"]) for r in by_g[key]])
        z = (mc - exact) / se
        worst_z = max(worst_z, float(np.max(np.abs(z))))
        ax.axhline(0.0, color=ORANGE, linewidth=1.5, zorder=1)
        ax.fill_between(s, -2, 2, color=GRID, alpha=0.6, linewidth=0, zorder=0)
        ax.plot(s, z, marker="o", markersize=5, linestyle="none", color=BLUE, zorder=3,
                markeredgecolor=SURFACE, markeredgewidth=1.0)
        ax.set_ylim(-4, 4)
        ax.set_title(title, loc="left")
        ax.set_xlabel("spot")
        if key in ("price", "vega"):
            ax.set_ylabel("(MC − analytic) / std error")
    suptitle(fig, "Monte Carlo against the closed form: errors in units of the reported "
                  "standard error",
             f"K = 100 call, T = 0.25, {paths:,} paths per point, an independent seed per spot; "
             f"bump-and-reprice Greeks with common random numbers.\nShaded band: ±2 s.e., where "
             f"~95% of points should fall if the error bars are right. Worst |z| = {worst_z:.2f}.")
    fig.savefig(out_dir / "mc_vs_analytic.png", dpi=150)
    plt.close(fig)


# ---- Implied-vol surface ------------------------------------------------------------

def vol_surface(data_dir, out_dir):
    nodes = read(data_dir / "iv_nodes.csv")
    grid = read(data_dir / "iv_surface.csv")
    quoted = sorted({float(r["expiry"]) for r in nodes})

    # 3D view of the regular (K, T) grid.
    reg = [r for r in grid if all(abs(float(r["expiry"]) - T) > 1e-9 for T in quoted)
           or abs(round(float(r["expiry"]) / 0.05) * 0.05 - float(r["expiry"])) < 1e-9]
    Ts = sorted({float(r["expiry"]) for r in reg
                 if abs(round(float(r["expiry"]) / 0.05) * 0.05 - float(r["expiry"])) < 1e-9})
    Ks = sorted({float(r["strike"]) for r in reg if float(r["strike"]).is_integer()})
    lookup = {(round(float(r["expiry"]), 4), float(r["strike"])): float(r["vol"]) for r in grid}
    K, T = np.meshgrid(Ks, Ts)
    V = np.vectorize(lambda k, t: lookup[(round(t, 4), k)])(K, T) * 100

    fig = plt.figure(figsize=(9, 6.5))
    ax = fig.add_subplot(projection="3d", facecolor=SURFACE)
    ax.plot_surface(K, T, V, cmap="Blues", vmin=V.min() - 8, vmax=V.max(),
                    linewidth=0.2, edgecolor=SURFACE, antialiased=True, alpha=0.95)
    nk = np.array([float(r["strike"]) for r in nodes])
    nt = np.array([float(r["expiry"]) for r in nodes])
    nv = np.array([float(r["vol"]) for r in nodes]) * 100
    ax.scatter(nk, nt, nv, color=ORANGE, s=10, depthshade=False, label="quoted (inverted) points")
    ax.set_xlabel("strike")
    ax.set_ylabel("expiry (years)")
    ax.set_zlabel("implied vol (%)")
    ax.view_init(elev=22, azim=-55)
    ax.set_box_aspect(None, zoom=0.88)  # keep the z label inside the figure
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(SURFACE)
        axis.pane.set_edgecolor(GRID)
    ax.legend(loc="upper right", bbox_to_anchor=(0.98, 0.88))
    fig.suptitle("Implied-volatility surface built from 105 synthetic quotes",
                 x=0.03, ha="left", y=0.995, fontsize=13, fontweight="semibold", color=INK)
    fig.text(0.03, 0.945, "Quotes priced from an arbitrage-free SSVI smile (S = 100, r = 3%, "
             "q = 1%), inverted, then interpolated\nlinearly in total variance; constant vol "
             "outside the quoted expiries and strikes.", ha="left", va="top", fontsize=9,
             color=INK_2, linespacing=1.4)
    fig.savefig(out_dir / "iv_surface.png", dpi=150, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)

    # Smiles per quoted expiry, and the interpolation error against the truth.
    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.85))
    for T, color in zip(quoted, RAMP_5):
        fine = sorted((float(r["strike"]), float(r["vol"]), float(r["true_vol"]))
                      for r in grid if abs(float(r["expiry"]) - T) < 1e-9)
        k = np.array([x[0] for x in fine])
        v = np.array([x[1] for x in fine]) * 100
        tv = np.array([x[2] for x in fine]) * 100
        pts = [(float(r["strike"]), float(r["vol"]) * 100) for r in nodes
               if abs(float(r["expiry"]) - T) < 1e-9]
        left.plot(k, v, color=color, label=expiry_label(T))
        left.plot([p[0] for p in pts], [p[1] for p in pts], "o", color=color, markersize=4,
                  markeredgecolor=SURFACE, markeredgewidth=0.8)
        left.annotate(expiry_label(T), (k[-1], v[-1]), xytext=(4, 0),
                      textcoords="offset points", fontsize=8, color=INK_2, va="center")
        right.plot(k, (v - tv) * 100, color=color, label=expiry_label(T))
    left.set_title("Smile at each quoted expiry", loc="left")
    left.set_xlabel("strike")
    left.set_ylabel("implied vol (%)")
    left.set_xlim(58, 168)
    right.axhline(0.0, color=AXIS, linewidth=0.8)
    right.set_title("Interpolation error: surface − true SSVI vol", loc="left")
    right.set_xlabel("strike")
    right.set_ylabel("error (basis points of vol)")
    right.legend(title="expiry", title_fontsize=8, loc="upper right")
    suptitle(fig, "What linear interpolation in total variance costs",
             "Dots: inverted quotes ($5 strikes). Error is zero at every quote and peaks "
             "between them, worst near the money at 1 month, where the smile curves most.\nIt is never negative: "
             "total variance is convex in log-moneyness, so a straight chord always lies "
             "above it.")
    fig.savefig(out_dir / "iv_smiles.png", dpi=150)
    plt.close(fig)


# ---- SVI fit vs linear interpolation on noisy quotes -------------------------------

def svi_fit(data_dir, out_dir):
    quotes = read(data_dir / "svi_quotes.csv")
    curves = read(data_dir / "svi_curves.csv")
    expiries = sorted({float(r["expiry"]) for r in quotes})

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.15, 1]})
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.88))
    for row, T in enumerate(expiries):
        q = [r for r in quotes if abs(float(r["expiry"]) - T) < 1e-9]
        c = [r for r in curves if abs(float(r["expiry"]) - T) < 1e-9]
        K = np.array([float(r["strike"]) for r in c])
        svi = np.array([float(r["svi_vol"]) for r in c]) * 100
        lin = np.array([float(r["linear_vol"]) for r in c]) * 100
        true = np.array([float(r["true_vol"]) for r in c]) * 100
        qk = [float(r["strike"]) for r in q]
        lo, hi = min(qk), max(qk)

        left, right = axes[row]
        for ax in (left, right):
            ax.axvspan(K[0], lo, color=GRID, alpha=0.5, linewidth=0, zorder=0)
            ax.axvspan(hi, K[-1], color=GRID, alpha=0.5, linewidth=0, zorder=0)
            ax.set_xlim(K[0], K[-1])
        left.plot(qk, [float(r["quote_vol"]) * 100 for r in q], "o", color=MUTED, markersize=3.5,
                  label="noisy quotes")
        left.plot(K, lin, color=ORANGE, linewidth=1.0, label="linear interpolation")
        left.plot(K, svi, color=BLUE, label="SVI fit")
        left.plot(K, true, color=INK, linewidth=1.0, linestyle=(0, (4, 3)), label="true smile")
        left.set_title(f"{expiry_label(T)} smile", loc="left")
        left.set_ylabel("implied vol (%)")

        right.axhline(0.0, color=AXIS, linewidth=0.8)
        right.plot(K, (lin - true) * 100, color=ORANGE, linewidth=1.0, label="linear")
        right.plot(K, (svi - true) * 100, color=BLUE, label="SVI")
        right.set_title(f"{expiry_label(T)}: error vs true smile", loc="left")
        right.set_ylabel("basis points of vol")
        right.set_ylim(-80, 80)
        if row == 1:
            left.set_xlabel("strike")
            right.set_xlabel("strike")
        if row == 0:
            left.legend(loc="upper right")
            right.legend(loc="upper right")
    suptitle(fig, "Fitting SVI instead of interpolating noisy quotes",
             "SSVI smile + 20 bp of Gaussian noise on every quote, $2.50 strikes. Shaded: "
             "outside the quoted strikes, where both curves extrapolate\n(linear: flat vol; "
             "SVI: its fitted wings). Inside, SVI's error is about half the noise; linear "
             "interpolation passes the noise straight through.")
    fig.savefig(out_dir / "svi_fit.png", dpi=150)
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: plot.py <csv_dir> <out_dir>")
    data_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    greeks_vs_spot(data_dir, out_dir)
    mc_vs_analytic(data_dir, out_dir)
    vol_surface(data_dir, out_dir)
    svi_fit(data_dir, out_dir)
    print(f"wrote 5 figures to {out_dir}")


if __name__ == "__main__":
    main()
