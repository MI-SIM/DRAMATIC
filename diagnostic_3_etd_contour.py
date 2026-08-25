"""
Diagnostic 3: validating the exponential-time-differencing (ETD) contour.

Rationale
---------
The timestepping scheme evaluates coefficient functions like

    phi1(z) = (-4 - z + exp(z)*(4 - 3z + z^2)) / z^2

which have a removable singularity at z=0 and are evaluated by numerical
contour integration (Cauchy's integral formula on a circle of radius R
around the origin, M quadrature points) rather than direct floating-point
evaluation, following Kassam & Trefethen (2005). This is *only* accurate
if two things hold:

1. The contour radius R and point count M are large/fine enough to
   resolve the phi-function accurately near z=0 -- this is a pure
   numerical-quadrature question, independent of the model.
2. **Critically for reuse on a new model**: in the original code, the
   "z" being integrated over is a free quadrature variable, used to
   evaluate phi-functions at z=0 *only* because the actual linear
   operator for the biomass-fraction PDE turned out to be zero for this
   particular model (verified separately -- see the summary). If a new
   model's analogous PDE has a genuinely non-trivial linear part (e.g. a
   real diffusion operator, not zero), the *correct* ETD recipe requires
   evaluating phi(dt * L) for the actual operator L, which means the
   contour must enclose the eigenvalues of dt*L -- not just the origin.
   Reusing the "collapsed to scalar" version from this code on a new
   model without checking this would silently give wrong answers with
   no error message.

This module provides two checks:

- `contour_quadrature_converged`: confirms the phi-function evaluation
  itself is converged in M (contour resolution) -- a pure numerics check,
  always worth running once for a new dt/parameter regime.
- `linear_operator_is_trivial`: an explicit test of assumption (2) above
  -- checks whether the "linear part" of your PDE's right-hand side is
  actually (numerically) zero, which is the precondition for reusing the
  scalar-coefficient shortcut instead of the full matrix-exponential
  version. This is the check that would have caught a silent misuse of
  this code pattern on a new problem.
"""

from __future__ import annotations

import numpy as np


def _phi_coefficients(z, M=64, radius=15.0):
    """Evaluate FS1, FS2, FS3, QS (Kassam-Trefethen ETDRK4-style scalar
    coefficients) at a single scalar z via M-point contour quadrature on
    a circle of the given radius centred at z, as in the original code.
    """
    r = radius * np.exp(1j * np.pi * ((np.arange(1, M + 1) - 0.5) / M))
    FS1 = FS2 = FS3 = QS = 0.0
    for zz in z + r:
        zIA = 1.0 / zz
        QS = QS + zIA * (np.exp(zz / 2.0) - 1.0)
        FS1 = FS1 + zIA * (-4 - zz + np.exp(zz) * (4 - 3 * zz + zz**2)) / zz**2
        FS2 = FS2 + zIA * (2 + zz + np.exp(zz) * (zz - 2)) / zz**2
        FS3 = FS3 + zIA * (-4 - 3 * zz - zz**2 + np.exp(zz) * (4 - zz)) / zz**2
    return (np.real(c / M) for c in (FS1, FS2, FS3, QS))


def contour_quadrature_converged(z=0.0, radius=15.0, M_values=(16, 32, 64, 128, 256),
                                  rtol=1e-8, verbose=True):
    """Check the phi-function contour quadrature has converged in M.

    Evaluates the ETD coefficients at increasing M and reports the
    relative change between successive M -- if it hasn't flattened out
    by your chosen M (default 64, matching the original code), that M is
    not resolving the contour integral accurately for your dt/parameters
    and should be increased.

    Returns True if the *last two* M values agree to rtol.
    """
    prev = None
    converged = False
    for M in M_values:
        coeffs = np.array(list(_phi_coefficients(z, M=M, radius=radius)))
        if prev is not None:
            rel_change = np.max(np.abs(coeffs - prev) / (np.abs(prev) + 1e-300))
            if verbose:
                print(f"  M={M:4d}:  FS1={coeffs[0]: .8e}  rel.change vs prev M = {rel_change:.3e}")
            converged = rel_change < rtol
        else:
            if verbose:
                print(f"  M={M:4d}:  FS1={coeffs[0]: .8e}")
        prev = coeffs
    if verbose:
        print(f"  -> {'CONVERGED' if converged else 'NOT YET CONVERGED'} to rtol={rtol:.0e}")
    return converged


def linear_operator_is_trivial(rhs_function, n, at=None, eps=1e-6, tol=1e-6,
                                verbose=True):
    """Test whether a PDE right-hand side's Jacobian is a scalar multiple
    of the identity -- the actual precondition for the collapsed
    scalar-coefficient ETD shortcut used in the original code (there,
    `zI - A` collapsed to `zI` because the biomass-fraction PDE's linear
    part turned out to be exactly zero, i.e. a scalar multiple -- zero --
    of the identity).

    Method: compute the numerical Jacobian J = d(rhs)/du at a reference
    state by finite differences, then check whether J is (numerically)
    diagonal *and* has a constant diagonal. A real spatial operator
    (diffusion, advection via a differentiation matrix) shows up as
    off-diagonal coupling between grid points -- exactly what this test
    flags. A purely pointwise/local reaction term (no coupling between
    grid points) gives a diagonal J; if that diagonal is also constant
    across the domain, J actually is a scalar multiple of the identity
    and the shortcut is valid. A diagonal-but-non-constant J (e.g. a
    spatially-varying local reaction rate) is *not* safe for the scalar
    shortcut either, even though it's not a differential operator.

    Parameters
    ----------
    rhs_function : callable
        Called as `rhs_function(u)` with u of shape (n,), returning the
        same shape.
    n : int
        State dimension (e.g. number of grid points).
    at : array_like, optional
        Reference state to linearize around (default: zeros).
    eps : float
        Finite-difference step for the Jacobian.
    tol : float
        Relative tolerance (to ||J||) for both the off-diagonal norm and
        the diagonal spread.

    Returns
    -------
    bool
        True if J is a scalar multiple of the identity to within `tol`
        (scalar ETD shortcut is valid); False otherwise (use the full
        matrix-exponential ETD treatment).
    """
    if at is None:
        at = np.zeros(n)
    f0 = rhs_function(at)
    J = np.zeros((n, n))
    for j in range(n):
        pert = np.zeros(n)
        pert[j] = eps
        J[:, j] = (rhs_function(at + pert) - f0) / eps

    diag = np.diag(J)
    off_diag_norm = np.linalg.norm(J - np.diag(diag))
    diag_spread = diag.max() - diag.min()
    scale = np.linalg.norm(J) + 1e-300

    off_diag_rel = off_diag_norm / scale
    diag_spread_rel = diag_spread / scale
    is_trivial = (off_diag_rel < tol) and (diag_spread_rel < tol)

    if verbose:
        print(f"  ||J||={np.linalg.norm(J):.4e}   "
              f"off-diagonal norm (relative)={off_diag_rel:.3e}   "
              f"diagonal spread (relative)={diag_spread_rel:.3e}   (tol={tol:.0e})")
        print(f"  -> J is {'a scalar multiple of I (scalar ETD shortcut OK)' if is_trivial else 'NOT a scalar multiple of I -- use full matrix ETD, not the scalar shortcut'}")
    return is_trivial


# --------------------------------------------------------------------------
# Demo
# --------------------------------------------------------------------------
if __name__ == "__main__":
    print("Contour quadrature convergence at z=0, radius=15 (matches the "
          "original PartialNitritation.m/partial_nitritation.py settings):")
    contour_quadrature_converged(z=0.0, radius=15.0)

    print()
    N = 33
    print("linear_operator_is_trivial demo, case A: purely local/pointwise "
          "reactive RHS (like this model's biomass-fraction PDE, after the "
          "advection/reaction terms cancel to leave no coupling between "
          "grid points):")
    def reactive_only_rhs(u):
        return u**2 - 0.1 * u          # nonlinear reaction only, no coupling
    linear_operator_is_trivial(reactive_only_rhs, n=N)

    print()
    print("linear_operator_is_trivial demo, case B: RHS with a real "
          "diffusion operator added (what a NEW model might actually have "
          "if the analogous PDE keeps a genuine spatial coupling term):")
    D2 = -2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)  # toy Laplacian
    def diffusive_rhs(u):
        return u**2 - 0.1 * u + 50.0 * (D2 @ u)   # + real diffusion term
    linear_operator_is_trivial(diffusive_rhs, n=N)
