"""
Diagnostic 5: validating spatial profiles independently of bulk output.

Rationale
---------
This whole investigation's central lesson: bulk/integrated quantities
(effluent concentrations, total biomass, thickness) can match a target
closely while the *internal spatial structure* (Fig. 2a) is qualitatively
wrong. Aggregate agreement is not evidence the spatial field is right --
it's evidence the *integral* of the spatial field is right, which is a
much weaker statement. For any new biofilm/microbial-ecology system where
within-biofilm structure matters (species stratification, niche
separation by depth, etc.), the spatial profile needs its own, separate
validation, run every time, not just when a bulk number looks suspicious.

Two checks:

1. `grid_convergence_profile` -- the spatial-profile analogue of a
   standard grid-convergence study. It's cheap and should be routine: run
   at N and 2N, and check the *full profile* agrees, not just the
   integrated/bulk output. (This is exactly the check that, in this
   session, correctly ruled out grid resolution as the cause of the flat
   Fig. 2a profile -- N=16 vs N=32 agreed to 3 decimals everywhere, so the
   flatness was real physics/model behaviour, not an artifact. Knowing
   that took five minutes; not knowing it would have wasted much more
   time chasing a numerics explanation for a modelling issue.)

2. `boundary_condition_residual` -- directly checks that imposed boundary
   conditions (e.g. no-flux/Neumann at an impermeable wall) are actually
   satisfied by the computed profile to near machine precision, rather
   than assuming the solver enforced them correctly just because it
   converged.
"""

from __future__ import annotations

import numpy as np


def grid_convergence_profile(solve_fn, N_values, x_key="x", profile_keys=("fa", "fn", "fi"),
                              interp_to_common_grid=True, verbose=True):
    """Run `solve_fn(N)` at each N in N_values and compare the resulting
    spatial profiles pairwise (not just their integrals/bulk summaries).

    Parameters
    ----------
    solve_fn : callable
        `solve_fn(N)` -> dict containing at least `x_key` (the grid,
        normalized 0-1 or any fixed domain) and each of `profile_keys`
        (1D arrays of the same length as x, evaluated at steady state or
        whatever single instant you want to check).
    N_values : sequence of int
        E.g. (16, 32, 64). At least 2 values required.
    interp_to_common_grid : bool
        Since different N give different (Chebyshev) grid points, profiles
        must be interpolated onto a common grid before comparing pointwise.

    Returns
    -------
    dict
        For each consecutive (N_i, N_{i+1}) pair and each profile key, the
        max absolute difference on a common fine grid. Also prints a
        pass/fail-style summary.
    """
    results = {}
    fine_grid = np.linspace(0, 1, 501)
    runs = {N: solve_fn(N) for N in N_values}

    for i in range(len(N_values) - 1):
        N1, N2 = N_values[i], N_values[i + 1]
        r1, r2 = runs[N1], runs[N2]
        for key in profile_keys:
            x1, x2 = np.asarray(r1[x_key]), np.asarray(r2[x_key])
            y1, y2 = np.asarray(r1[key]), np.asarray(r2[key])
            order1, order2 = np.argsort(x1), np.argsort(x2)
            if interp_to_common_grid:
                y1i = np.interp(fine_grid, x1[order1], y1[order1])
                y2i = np.interp(fine_grid, x2[order2], y2[order2])
            else:
                y1i, y2i = y1[order1], y2[order2]
            max_diff = float(np.max(np.abs(y1i - y2i)))
            results[(N1, N2, key)] = max_diff
            if verbose:
                print(f"  N={N1:4d} vs N={N2:4d}, profile '{key}': "
                      f"max|Δ| over depth = {max_diff:.3e}")

    if verbose:
        worst = max(results.values())
        print(f"  -> Largest cross-N profile difference: {worst:.3e}. "
              "If this is not small compared to the *feature size you care "
              "about* (e.g. the spread you're trying to resolve in a "
              "stratification profile), you are NOT yet grid-converged -- "
              "increase N further before trusting the profile shape.")
    return results


def boundary_condition_residual(profile, D_matrix, bc_type="neumann", bc_node=-1, tol=1e-8,
                                 verbose=True, label=""):
    """Check that an imposed boundary condition is actually satisfied by
    a computed profile, to near machine precision -- not just assumed
    because the solver "converged".

    Parameters
    ----------
    profile : array_like
        The computed field (e.g. a substrate concentration profile) at
        the grid points corresponding to `D_matrix`'s rows/columns.
    D_matrix : array_like, shape (n, n)
        Differentiation matrix for the same grid (e.g. the Chebyshev D
        matrix), used to evaluate the derivative at the boundary node.
    bc_type : {"neumann", "dirichlet"}
        "neumann": checks D_matrix[bc_node, :] @ profile ~ 0 (no-flux).
        "dirichlet": checks profile[bc_node] against an expected value
        (pass `tol` generously and check the printed value by eye, or
        extend this for a specific expected value if you have one).
    bc_node : int
        Index of the boundary node to check (e.g. -1 or 0 depending on
        your grid convention).
    tol : float
        Tolerance for the Neumann residual.

    Raises
    ------
    RuntimeError
        If a Neumann condition is violated beyond tolerance.
    """
    profile = np.asarray(profile)
    D_matrix = np.asarray(D_matrix)
    if bc_type == "neumann":
        residual = float(D_matrix[bc_node, :] @ profile)
        if verbose:
            print(f"  Neumann BC residual{' (' + label + ')' if label else ''} "
                  f"at node {bc_node}: {residual:.3e} (tol {tol:.0e})")
        if abs(residual) > tol:
            raise RuntimeError(
                f"Neumann (no-flux) boundary condition violated{' for ' + label if label else ''}: "
                f"D @ profile = {residual:.3e} at node {bc_node}, expected ~0 "
                f"(tol {tol:.0e}). The solver converged numerically but is not "
                "actually satisfying the imposed physical boundary condition -- "
                "check the BC row substitution in your discretization."
            )
    elif bc_type == "dirichlet":
        if verbose:
            print(f"  Dirichlet BC value{' (' + label + ')' if label else ''} "
                  f"at node {bc_node}: {profile[bc_node]:.6e} -- compare by eye "
                  "to the value you imposed.")
    else:
        raise ValueError(f"Unknown bc_type: {bc_type!r}")


# --------------------------------------------------------------------------
# Demo against real saved data from this session.
# --------------------------------------------------------------------------
if __name__ == "__main__":
    import numpy as np

    print("=== grid_convergence_profile, using this session's saved N=16 "
          "and N=32 runs (fA profile) ===")

    d16 = np.load("/tmp/run_test2.npz") if __import__("os").path.exists("/tmp/run_test2.npz") else None
    d32 = np.load("/tmp/run_hires_N32_T60.npz") if __import__("os").path.exists("/tmp/run_hires_N32_T60.npz") else None

    if d16 is not None and d32 is not None:
        def fake_solve(N):
            d = d16 if N == 16 else d32
            return {"x": d["x"], "fa": d["fafa"][:, -1], "fn": d["fnfn"][:, -1], "fi": d["fifi"][:, -1]}
        grid_convergence_profile(fake_solve, N_values=(16, 32), profile_keys=("fa", "fn", "fi"))
    else:
        print("  (skipped -- saved runs not found in this environment)")

    print()
    print("=== boundary_condition_residual, using this session's Chebyshev "
          "D matrix and a converged CO2 profile (Neumann at the substratum) ===")
    import sys
    sys.path.insert(0, "/tmp/work")
    from cheb import cheb
    N = 32
    D, x = cheb(N)
    if d32 is not None:
        CO2_profile = d32["CO2C"][:, -1]
        # In this model's convention, the substratum (no-flux) node is index N (last row/col).
        boundary_condition_residual(CO2_profile, D, bc_type="neumann", bc_node=N,
                                     tol=1e-6, label="O2 at substratum")
    else:
        print("  (skipped -- saved run not found)")
