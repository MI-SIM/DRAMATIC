"""
Diagnostic 2: fixed-point iteration monitoring for the inner BVP solve.

Rationale
---------
The biofilm substrate profile is found each timestep by a fixed-point
(Anderson-accelerated) iteration: `biofilmbvp()` in the original code
loops `while change > tol`, with no cap on iteration count and no record
of the convergence history. Two failure modes this misses completely:

1. Slow convergence that never technically "fails" -- it just eats an
   enormous amount of wall-clock time on a harder problem (more species,
   stiffer kinetics, different diffusion/reaction ratio) without any
   signal that something is off, until the whole run is impractically slow.
2. Silent non-convergence -- e.g. an off-by-one in the Anderson history
   window causes oscillation around a fixed cycle rather than convergence;
   `change` stops decreasing but never technically exceeds a runaway
   threshold either. With no cap, this can hang forever, or on some
   platforms silently return whatever the loop last computed if a NaN
   check downstream mistakes lack of change for convergence.

This wraps any "generic fixed point solve for a nonlinear system" pattern
with: (a) a hard iteration cap that raises with real diagnostic content
instead of hanging or looping forever, (b) a convergence-history record you
can inspect or plot, (c) a slow-convergence *warning* (not just outright
failure) so you notice degrading performance before it becomes a hang.

This is written as a decorator/wrapper around your existing iteration
body, so you can drop it onto a new model's BVP solver without rewriting
the solver itself -- you only need to expose one function that does "one
fixed-point update step" and returns the new iterate plus the residual.
"""

from __future__ import annotations

import warnings
import numpy as np


class FixedPointDivergedError(RuntimeError):
    """Raised when a fixed-point iteration fails to converge in time."""


class FixedPointMonitor:
    """Wraps a fixed-point iteration with iteration cap + convergence log.

    Example
    -------
    >>> monitor = FixedPointMonitor(max_iter=200, tol=1e-6, warn_iter=50)
    >>> C = C_initial_guess
    >>> while True:
    ...     C_new, residual = one_fixed_point_step(C)   # your existing update
    ...     if monitor.step(residual):                  # True -> converged
    ...         break
    ...     C = C_new
    >>> monitor.history          # list of residuals, for plotting/inspection
    """

    def __init__(self, max_iter=500, tol=1e-6, warn_iter=None, label=""):
        self.max_iter = max_iter
        self.tol = tol
        self.warn_iter = warn_iter or max(20, max_iter // 10)
        self.label = label
        self.history = []
        self._warned = False

    def step(self, residual):
        """Record one iteration's residual. Returns True if converged.

        Raises FixedPointDivergedError if max_iter is exceeded, or if the
        residual has clearly stopped improving (stalled) well above tol.
        """
        self.history.append(float(residual))
        n = len(self.history)

        if not self._warned and n == self.warn_iter and residual > self.tol:
            warnings.warn(
                f"Fixed-point iteration{' (' + self.label + ')' if self.label else ''} "
                f"has not converged after {n} iterations (residual={residual:.3e}, "
                f"tol={self.tol:.1e}). Continuing, but this is unusually slow -- "
                f"check for a stiffer-than-expected problem or a solver bug.",
                RuntimeWarning,
            )
            self._warned = True

        if residual <= self.tol:
            return True

        if n >= self.max_iter:
            raise FixedPointDivergedError(
                f"Fixed-point iteration{' (' + self.label + ')' if self.label else ''} "
                f"did not converge in {self.max_iter} iterations "
                f"(final residual={residual:.3e}, tol={self.tol:.1e}). "
                f"Residual history (last 10): {self.history[-10:]}. "
                "This usually means either the problem is stiffer than the "
                "solver/tolerance combination can handle, or the Anderson "
                "history window (mmax) needs to be larger."
            )

        # Stall detection: last 10 iterations made < 1% relative progress.
        if n >= 20:
            recent = self.history[-10:]
            if recent[0] > 0 and abs(recent[0] - recent[-1]) / recent[0] < 0.01:
                raise FixedPointDivergedError(
                    f"Fixed-point iteration{' (' + self.label + ')' if self.label else ''} "
                    f"appears stalled: residual only moved from {recent[0]:.3e} to "
                    f"{recent[-1]:.3e} over the last 10 iterations, well above "
                    f"tol={self.tol:.1e}. Likely a cycling/oscillation failure "
                    "in the acceleration scheme rather than slow-but-real progress."
                )
        return False


# --------------------------------------------------------------------------
# Demo: wrap a toy fixed-point iteration to show normal convergence,
# a stall, and a genuine non-convergence, each correctly distinguished.
# --------------------------------------------------------------------------
if __name__ == "__main__":
    # Case 1: healthy convergence
    monitor = FixedPointMonitor(max_iter=100, tol=1e-8, label="healthy")
    r = 1.0
    for _ in range(100):
        r *= 0.5
        if monitor.step(r):
            break
    print(f"[healthy] converged in {len(monitor.history)} iterations")

    # Case 2: stalled/cycling (residual barely moves)
    monitor = FixedPointMonitor(max_iter=100, tol=1e-8, warn_iter=15, label="stalled")
    try:
        r = 1.0
        for i in range(100):
            r = 0.5 + 0.001 * np.sin(i)  # never actually shrinks
            monitor.step(r)
    except FixedPointDivergedError as e:
        print(f"[stalled] correctly raised: {e}")

    # Case 3: genuinely too slow, hits the cap
    monitor = FixedPointMonitor(max_iter=30, tol=1e-8, warn_iter=10, label="too-slow")
    try:
        r = 1.0
        for _ in range(100):
            r *= 0.9  # converges, but not within 30 iterations at this rate
            monitor.step(r)
    except FixedPointDivergedError as e:
        print(f"[too-slow] correctly raised: {e}")
