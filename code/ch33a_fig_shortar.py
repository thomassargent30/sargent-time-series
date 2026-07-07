"""
Chapter 33a, short-autoregression figure (reconstruction of Hansen-Sargent 1993,
Figs. 6 and 7, Section 5.4).

TRUE model: time-invariant with a *seasonal* endowment process,
   gamma=0.1, phi1=0.3, deltak=0.95, beta=1/1.05,  b_t: (1-0.2L), sigma=0.25,
   endowment  d'_t = 0.1 d'_{t-1} + 0.5 d'_{t-4} - 0.4 d'_{t-5} + w2_t  (sigma=1),
whose lag-4/5 terms put seasonal peaks (notably at omega=pi) into consumption,
investment, and the endowment itself.

The econometrician misspecifies by fitting a first-order autoregression for the
endowment.  Reported estimates:
   seasonally adjusted   : gamma=0.1000, phi1=0.3011, a1= 0.1631
   seasonally unadjusted  : gamma=0.1003, phi1=0.2964, a1=-0.2018
The unadjusted fit drives a1 negative in a vain attempt to match the peak at
omega=pi; neither AR(1) can reproduce the seasonal peak.

Output: ../figures/ch33a_shortar.png
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

true = solve(build_model(0.1, 0.3, DELTAK, BETA, [0.2], 0.25,
                         [0.1, 0.0, 0.0, 0.5, -0.4], 1.0))
adj = solve(build_model(0.1000, 0.3011, DELTAK, BETA, [0.2], 0.25, [0.1631], 1.0))
una = solve(build_model(0.1003, 0.2964, DELTAK, BETA, [0.2], 0.25, [-0.2018], 1.0))

plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, key, title in [(axes[0], "c", "Consumption"),
                       (axes[1], "d", "Endowment")]:
    St = np.real(spectrum(true['Acl'], true['C'], true[key], omega)[:, 0, 0])
    Sa = np.real(spectrum(adj['Acl'], adj['C'], adj[key], omega)[:, 0, 0])
    Su = np.real(spectrum(una['Acl'], una['C'], una[key], omega)[:, 0, 0])
    for w0 in SEASONAL_OMEGA:
        ax.axvline(w0, color="0.75", ls=":", lw=0.9, zorder=0)
    ax.semilogy(omega, St, color="C0", lw=2.0, label="True (seasonal endowment)")
    ax.semilogy(omega, Sa, color="C2", lw=1.7, ls="--",
                label="Approx.: AR(1), adjusted")
    ax.semilogy(omega, Su, color="C3", lw=1.7, ls="-.",
                label="Approx.: AR(1), unadjusted")
    ax.set_title(title)
    ax.set_xlabel(r"Frequency $\omega$")
    ax.set_xlim(0, np.pi)
    ax.set_xticks([0, np.pi / 2, np.pi])
    ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])

axes[0].set_ylabel("Spectral density (log scale)")
axes[0].legend(loc="lower center", fontsize=9, framealpha=0.9)
fig.suptitle("Too short an endowment autoregression: AR(1) misses the seasonal "
             "peak (Hansen–Sargent 1993, Figs. 6 & 7)", fontsize=12)
fig.tight_layout()
out = FIG_DIR / "ch33a_shortar.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print("wrote", out)
