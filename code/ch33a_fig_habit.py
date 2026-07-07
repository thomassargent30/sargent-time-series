"""
Chapter 33a, habit-persistence figure (reconstruction of Hansen-Sargent 1993,
Figs. 4 and 8).

TRUE model: time-invariant with seasonal habit persistence (Section 5.2 / 5.5),
   gamma=0.1, phi1=0.005, lambda=0.8, deltah=0.9, deltak=0.95, beta=1/1.05,
   b_t: (1-0.2L), sigma_w1=0.25;   d_t: (1-0.7L), sigma_w2=1.
The seasonal habit induces sharp spectral peaks at the quarterly seasonal
frequencies (omega = pi/2, pi) in consumption and investment.

Two approximating models are overlaid:
  * correctly specified (5.2): identical parameters -> spectrum lies atop the
    truth (this is why every data treatment recovers the truth in Table of 5.2);
  * habit omitted (5.5): lambda=0, deltah=0, and the reported estimates
    phi1=0.0043, a1(endowment)=0.4019 -> misses the seasonal peaks.

Output: ../figures/ch33a_habit.png
"""
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ch33a_seasonality_lq import build_model, solve, spectrum, SEASONAL_OMEGA

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

BETA, DELTAK = 1 / 1.05, 0.95
omega = np.linspace(0.02, np.pi, 1400)

# --- true model: seasonal habit persistence ---
true = solve(build_model(0.1, 0.005, DELTAK, BETA, [0.2], 0.25, [0.7], 1.0,
                         lam=0.8, deltah=0.9, p=4))
# --- correctly specified approximating model (coincides with truth) ---
correct = true
# --- misspecified: habit omitted, endowment AR(1) at reported estimate ---
nohabit = solve(build_model(0.1, 0.0043, DELTAK, BETA, [0.2], 0.25, [0.4019], 1.0))

plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, key, title in [(axes[0], "c", "Consumption"),
                       (axes[1], "i", "Investment")]:
    St = np.real(spectrum(true['Acl'], true['C'], true[key], omega)[:, 0, 0])
    Sn = np.real(spectrum(nohabit['Acl'], nohabit['C'], nohabit[key], omega)[:, 0, 0])
    for w0 in SEASONAL_OMEGA:
        ax.axvline(w0, color="0.75", ls=":", lw=0.9, zorder=0)
    ax.semilogy(omega, St, color="C0", lw=2.0, label="True (habit persistence)")
    ax.semilogy(omega, St, color="k", lw=1.0, ls=(0, (1, 2)),
                label="Approx.: correctly specified", zorder=3)
    ax.semilogy(omega, Sn, color="C3", lw=1.8, ls="--",
                label="Approx.: habit omitted (5.5)")
    ax.set_title(title)
    ax.set_xlabel(r"Frequency $\omega$")
    ax.set_xlim(0, np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi])
    ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])

axes[0].set_ylabel("Spectral density (log scale)")
axes[0].legend(loc="upper right", fontsize=9, framealpha=0.9)
fig.suptitle("Seasonal habit persistence: correctly specified vs. habit omitted "
             "(Hansen–Sargent 1993, Figs. 4 & 8)", fontsize=12)
fig.tight_layout()
out = FIG_DIR / "ch33a_habit.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print("wrote", out)
