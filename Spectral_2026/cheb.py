"""Chebyshev differentiation matrix (Trefethen)."""

import numpy as np


def cheb(N):
    """Return D, x for the Chebyshev collocation grid of order N."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * ((-1.0) ** np.arange(N + 1))
    X = np.tile(x[:, None], (1, N + 1))
    dX = X - X.T
    D = (c[:, None] * (1.0 / c)[None, :]) / (dX + np.eye(N + 1))
    D = D - np.diag(np.sum(D, axis=1))
    return D, x
