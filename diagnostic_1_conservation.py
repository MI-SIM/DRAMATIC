"""
Diagnostic 1: conservation and mass-balance self-consistency checks.

Rationale
---------
Bulk output matching a target (a paper's figure, an experiment) is
necessary but not sufficient evidence the code is right -- a bug can
still cancel out or stay too small to move the headline numbers you
happen to be checking. Conservation laws the *equations themselves*
guarantee mathematically (fractions summing to 1, total element mass
balance) should hold to numerical precision, independent of whether you
have a reference figure to compare against. Applicable to essentially
any reaction-transport model with compositional fractions and a
conserved element (N, C, COD, ...).

Two checks are implemented:

1. `check_fraction_conservation` -- for models where state variables are
   volume/area fractions that must sum to 1 pointwise (e.g. biofilm
   species fractions f_A + f_N + f_I = 1 in Masic & Eberl). Catches
   integration-scheme bugs (e.g. a back-substituted species computed
   inconsistently with the others) that a bulk check would never see,
   because the error can be exactly zero in the *sum* while very wrong
   in the *individual* profiles.

2. `check_element_mass_balance` -- for the whole system: over any time
   window, (element mass in) - (element mass out) should equal
   (element mass accumulated), to within integration error. This is a
   template: you supply small callables that compute each piece from
   your model's own state, because those are necessarily model-specific,
   but the checking logic itself is generic.

Usage pattern: call these as *assertions inside your timestepping loop*
(every K steps, not necessarily every step, for cost), not as a one-off
post-hoc audit. A check that only runs once, after the fact, only tells
you the final state was fine; a check embedded in the loop tells you the
first moment things went wrong, which is what you actually need to debug.
"""

from __future__ import annotations

import numpy as np


def check_fraction_conservation(fractions, tol=1e-8, t=None, label=""):
    """Assert that fraction-type state variables sum to 1 pointwise.

    Parameters
    ----------
    fractions : array_like, shape (n_points, n_species)
        E.g. np.column_stack([fA, fN, fI]) at one instant in time.
    tol : float
        Absolute tolerance on |sum - 1|.
    t : float, optional
        Current simulation time, only used for the error message.
    label : str
        Optional tag for the error message (e.g. "biofilm fractions").

    Raises
    ------
    RuntimeError
        If any point violates the tolerance, reporting the worst offender
        and its location -- not just "it failed somewhere".
    """
    fractions = np.asarray(fractions)
    totals = fractions.sum(axis=1)
    resid = np.abs(totals - 1.0)
    worst = np.argmax(resid)
    if resid[worst] > tol:
        when = f" at t={t}" if t is not None else ""
        raise RuntimeError(
            f"Fraction conservation violated{when}{': ' + label if label else ''}. "
            f"Worst residual |sum-1| = {resid[worst]:.3e} at grid index {worst} "
            f"(values there: {fractions[worst]}). Tolerance was {tol:.1e}."
        )


def check_element_mass_balance(
    element_in_rate,
    element_out_rate,
    element_stock_before,
    element_stock_after,
    dt,
    tol_rel=1e-4,
    t=None,
):
    """Assert (in - out)*dt equals the change in total element inventory.

    You provide the four physical quantities directly (in your model's
    own units) -- this function only does the bookkeeping and raises a
    clear, actionable error if they don't add up.

    Parameters
    ----------
    element_in_rate, element_out_rate : float
        Total mass rate of the tracked element (e.g. total N) entering
        / leaving the *whole system* (influent + any gas stripping etc,
        and effluent) at the *start* of the step, in mass/time.
    element_stock_before, element_stock_after : float
        Total mass of the element summed over *every* pool that holds it
        (bulk dissolved species, suspended biomass N content, biofilm N
        content, ...) before and after the step.
    dt : float
        Step size (same time units as the rates).
    tol_rel : float
        Relative tolerance on the imbalance, relative to the larger of
        the two stock values (avoids false alarms when both are tiny).

    Raises
    ------
    RuntimeError
        If the imbalance exceeds tolerance, printing the actual and
        expected change so you can see the sign and rough magnitude of
        whatever's wrong.
    """
    expected_change = (element_in_rate - element_out_rate) * dt
    actual_change = element_stock_after - element_stock_before
    scale = max(abs(element_stock_before), abs(element_stock_after), 1e-30)
    imbalance = abs(actual_change - expected_change)
    if imbalance / scale > tol_rel:
        when = f" at t={t}" if t is not None else ""
        raise RuntimeError(
            f"Element mass balance violated{when}. "
            f"Expected change (in-out)*dt = {expected_change:.6e}, "
            f"actual change in stock = {actual_change:.6e}, "
            f"relative imbalance = {imbalance/scale:.3e} (tol {tol_rel:.1e})."
        )


# --------------------------------------------------------------------------
# Worked example against the Masic & Eberl biofilm fractions, using data
# already on disk from this session, to show these actually catch things.
# --------------------------------------------------------------------------
if __name__ == "__main__":
    import glob

    npz_paths = ["/tmp/run_hires_N32_T60.npz"] if __import__("os").path.exists(
        "/tmp/run_hires_N32_T60.npz"
    ) else sorted(glob.glob("/tmp/work/run_asymmetric.npz"))
    if not npz_paths:
        print("No saved run found to demo against; run partial_nitritation() "
              "first and save fafa/fnfn/fifi to try this for real.")
    else:
        d = np.load(npz_paths[0])
        fafa, fnfn, fifi = d["fafa"], d["fnfn"], d["fifi"]
        tt = d["tt"]
        n_checked = 0
        for j in range(0, fafa.shape[1], 50):  # spot-check every 50th saved step
            frac = np.column_stack([fafa[:, j], fnfn[:, j], fifi[:, j]])
            check_fraction_conservation(frac, tol=1e-6, t=tt[j], label="fA+fN+fI")
            n_checked += 1
        print(f"check_fraction_conservation passed at {n_checked} time points "
              f"from {npz_paths[0]}.")
