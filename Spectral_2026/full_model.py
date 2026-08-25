"""Driver for partial_nitritation(): reproduces Fig. 1(a-d) and Fig. 2(a-b)
of Masic & Eberl (2014), Bull. Math. Biol. 76:27-58, for the D=5/day,
c=50-carrier scenario.

Requires the fixed partial_nitritation.py (nc=50, Table-1-consistent
Rxa/Rxn -- see accompanying summary for why).
"""

import matplotlib.pyplot as plt
import numpy as np

from partial_nitritation import partial_nitritation

RHO = 10000.0          # biomass density (g/m^3), must match partial_nitritation.py
AR, AC, NC = 0.17, 0.0068, 50
A = AR + AC * NC        # colonizable area (m^2) for c=50


# ---------------------------------------------------------------- Fig. 1a-c
def plot_fig1a(t, S):
    """Effluent NH4 and NO2 vs time. (NO3 is not tracked -- see note below.)"""
    plt.figure()
    plt.plot(t, S[0, :], "k-", linewidth=2, label=r"NH$_4$")
    plt.plot(t, S[1, :], "k:", linewidth=2, label=r"NO$_2$")
    plt.xlabel("time (d)")
    plt.ylabel(r"effluent reactor concentrations (g/m$^3$)")
    plt.title("Fig. 1a: effluent concentrations")
    plt.legend()


def plot_fig1b(t, S):
    """Suspended AOB / NOB / inert biomass vs time."""
    plt.figure()
    plt.plot(t, S[2, :], "k-", linewidth=2, label="AOB")
    plt.plot(t, S[3, :], "k:", linewidth=2, label="NOB")
    plt.plot(t, S[4, :], "k-.", linewidth=2, label="inert")
    plt.xlabel("time (d)")
    plt.ylabel("suspended biomass (g)")
    plt.title("Fig. 1b: suspended biomass")
    plt.legend()


def plot_fig1c(t, S):
    """Biofilm thickness vs time."""
    plt.figure()
    plt.plot(t, S[5, :] * 1e6, "k-", linewidth=2)
    plt.xlabel("time (d)")
    plt.ylabel(r"biofilm thickness ($\mu$m)")
    plt.title("Fig. 1c: biofilm thickness")


# --------------------------------------------------------------------- Fig 1d
def integrate_biofilm_biomass(x, S, fa, fn, fi, rho=RHO, area=A):
    """Integrate fA/fN/fI over depth at every saved time point to get
    total biofilm AOB/NOB/inert biomass (g) as a function of time.

    fa, fn, fi: arrays of shape (N+1, n_saved_times) as returned by
    partial_nitritation(). x: depth-fraction grid (also returned).
    """
    order = np.argsort(x)
    xs = x[order]
    L = S[5, :]
    n_t = fa.shape[1]
    mA = np.empty(n_t)
    mN = np.empty(n_t)
    mI = np.empty(n_t)
    for j in range(n_t):
        avgA = np.trapezoid(fa[order, j], xs)
        avgN = np.trapezoid(fn[order, j], xs)
        avgI = np.trapezoid(fi[order, j], xs)
        mA[j] = rho * area * L[j] * avgA
        mN[j] = rho * area * L[j] * avgN
        mI[j] = rho * area * L[j] * avgI
    return mA, mN, mI


def plot_fig1d(t, S, x, fa, fn, fi):
    """Biofilm AOB / NOB / inert biomass vs time."""
    mA, mN, mI = integrate_biofilm_biomass(x, S, fa, fn, fi)
    plt.figure()
    plt.plot(t, mA, "k-", linewidth=2, label="AOB")
    plt.plot(t, mN, "k:", linewidth=2, label="NOB")
    plt.plot(t, mI, "k-.", linewidth=2, label="inert")
    plt.xlabel("time (d)")
    plt.ylabel("biofilm biomass (g)")
    plt.title("Fig. 1d: biofilm biomass")
    plt.legend()
    return mA, mN, mI


# --------------------------------------------------------------------- Fig 2a
def plot_fig2a(x, S, fa, fn, fi):
    """Cumulative biomass fractions vs biofilm depth, at steady state
    (last saved time point). x-axis in um, 0 = substratum, L = surface,
    matching the paper's convention."""
    L_um = S[5, -1] * 1e6
    z_um = x * L_um
    order = np.argsort(z_um)
    z = z_um[order]

    plt.figure()
    plt.stackplot(
        z, fn[order, -1], fi[order, -1], fa[order, -1],
        labels=["NOB", "inert", "AOB"],
        colors=["#1f3d7a", "white", "#b0b0b0"], edgecolor="k", linewidth=0.5,
    )
    plt.xlim(0, L_um)
    plt.ylim(0, 1)
    plt.xlabel(r"biofilm thickness ($\mu$m)")
    plt.ylabel("cumulative biomass fractions")
    plt.title("Fig. 2a: biomass fractions in the biofilm")
    plt.legend(loc="lower right")


# --------------------------------------------------------------------- Fig 2b
def plot_fig2b(x, S, CNH4, CNO2, CO2):
    """Dissolved substrate concentrations vs biofilm depth, at steady state.

    NO3 is not plotted: it's fully decoupled from every other state
    variable in this model (never feeds back into any growth or
    consumption term), so partial_nitritation() doesn't track it.
    """
    L_um = S[5, -1] * 1e6
    z_um = x * L_um
    order = np.argsort(z_um)
    z = z_um[order]

    plt.figure()
    plt.plot(z, CO2[order, -1], "-.", color="gray", linewidth=2, label=r"O$_2$")
    plt.plot(z, CNH4[order, -1], "-", color="k", linewidth=2, label=r"NH$_4$")
    plt.plot(z, CNO2[order, -1], ":", color="k", linewidth=2, label=r"NO$_2$")
    plt.xlim(0, L_um)
    plt.xlabel(r"biofilm thickness ($\mu$m)")
    plt.ylabel(r"dissolved components C(z) (g/m$^3$)")
    plt.title("Fig. 2b: dissolved components in the biofilm")
    plt.legend()


# ------------------------------------------------------------------------ main
def main():
    T = 90
    N = 32

    t, S, CNH4, CNO2, CO2, fa, fn, fi, x = partial_nitritation(T, N)

    plot_fig1a(t, S)
    plot_fig1b(t, S)
    plot_fig1c(t, S)
    plot_fig1d(t, S, x, fa, fn, fi)
    plot_fig2a(x, S, fa, fn, fi)
    plot_fig2b(x, S, CNH4, CNO2, CO2)

    plt.show()


if __name__ == "__main__":
    main()
