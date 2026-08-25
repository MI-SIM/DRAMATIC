"""Spectral Chebyshev integration matrices."""

import numpy as np


def intmat(N):
    """Return (Sl, Sr) of order N+1.

    Sl integrates from x to 1; Sr integrates from -1 to x.
    N must be even.
    """
    j = np.arange(N // 2 + 1)[:, None]
    Ut = np.zeros((N // 2 + 1, N + 1))
    Ut[:, 0] = np.cos(j[:, 0] * np.pi / N)
    Ut[:, 1] = 0.5 * np.cos(j[:, 0] * np.pi / N) ** 2
    cols_hi = np.arange(3, N + 2)
    cols_lo = np.arange(1, N)
    Ut[:, 2:] = (
        np.cos(j * cols_hi * np.pi / N) / cols_hi
        - np.cos(j * cols_lo * np.pi / N) / cols_lo
    ) / 2.0

    stacked = np.hstack([Ut, Ut[:, N - 1:0:-1]])
    B = np.real(np.fft.fft(stacked.T, axis=0)).T
    At = np.hstack(
        [B[: N // 2 + 1, 0:1] / 2.0, B[: N // 2 + 1, 1:N], B[: N // 2 + 1, N : N + 1] / 2.0]
    ) / N
    A = np.vstack([At, -np.rot90(At[: N // 2, :], 2)])
    Sl = A - np.ones((N + 1, 1)) * A[N, :]
    Sr = np.rot90(Sl, 2)
    return Sl, Sr
