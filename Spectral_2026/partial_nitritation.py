"""Anammox and partial nitritation in a biofilm (Masic and Eberl, 2015).

Python translation of PartialNitritation.m. Nested closures keep the same
parameter-passing pattern as the MATLAB function.
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import norm, solve

from cheb import cheb
from intmat import intmat


def _rcond(A):
    c = np.linalg.cond(A)
    if not np.isfinite(c) or c == 0.0:
        return 0.0
    return 1.0 / c


def _qrdelete(Q, R):
    """Givens QR update after dropping the first column of R."""
    m = R.shape[0]
    n = Q.shape[0]
    Q = Q.copy()
    R = R.copy()
    for kk in range(m - 1):
        temp = np.sqrt(R[kk, kk + 1] ** 2 + R[kk + 1, kk + 1] ** 2)
        c = R[kk, kk + 1] / temp
        s = R[kk + 1, kk + 1] / temp
        R[kk, kk + 1] = temp
        R[kk + 1, kk + 1] = 0.0
        if kk < m - 2:
            for jj in range(kk + 2, m):
                temp = c * R[kk, jj] + s * R[kk + 1, jj]
                R[kk + 1, jj] = -s * R[kk, jj] + c * R[kk + 1, jj]
                R[kk, jj] = temp
        for ll in range(n):
            temp = c * Q[ll, kk] + s * Q[ll, kk + 1]
            Q[ll, kk + 1] = -s * Q[ll, kk] + c * Q[ll, kk + 1]
            Q[ll, kk] = temp
    return Q[:, : m - 1], R[: m - 1, 1:m]


def partial_nitritation(T, N):
    """Integrate the biofilm model to time T on N Chebyshev points.

    Returns tt, SS, CNH4C, CNO2C, CO2C, fafa, fnfn, fifi, x
    with x mapped from Chebyshev [-1, 1] onto biofilm depth [0, 1].
    """
    alpha = 1.0
    Ac = 0.0068
    Ar = 0.17
    d = 5.0
    E = 1000.0
    eta = 0.5
    fxi = 0.1
    nc = 50
    rho = 10000.0
    V = 0.00603186

    A = Ar + Ac * nc

    dNO2 = 7.93e-5
    dNH4 = 9.13e-5
    dO2 = 9.93e-5

    kaO2 = 0.5
    knO2 = 0.5
    kaNH4 = 0.1690
    knNO2 = 0.302
    muamax = 0.3082
    munmax = 0.4015
    ya = 0.15
    yn = 0.041
    ba = 0.04
    bn = 0.08
    ia = 0.07
    ii = 0.02

    SNH4in = 30.0
    SNO2in = 0.0
    SO2in = 5.0

    D, x = cheb(N)
    Dx = D.copy()
    Dx[N, :] = 0.0
    D2 = D @ D
    D2[N, :] = D[N, :]
    D2 = D2[1:, 1:]
    D2inv = np.linalg.inv(D2)
    dt = 0.001
    t = 0.0
    Int, _ = intmat(N)

    tmax = T
    points = 9000
    tplot = T / points
    plotgap = max(1, int(round(tplot / dt)))
    dt = tplot / plotgap
    nplots = int(round(tmax / dt))
    tol = 1e-6

    # L=0 contour integrals reduce to scalar ETD coefficients times identity.
    M = 64
    r = 15.0 * np.exp(1j * np.pi * ((np.arange(1, M + 1) - 0.5) / M))
    h = dt
    FS1 = 0.0
    FS2 = 0.0
    FS3 = 0.0
    QS = 0.0
    for z in r:
        zIA = 1.0 / z
        QS = QS + h * zIA * (np.exp(z / 2.0) - 1.0)
        FS1 = FS1 + h * zIA * (-4 - z + np.exp(z) * (4 - 3 * z + z**2)) / z**2
        FS2 = FS2 + h * zIA * (2 + z + np.exp(z) * (z - 2)) / z**2
        FS3 = FS3 + h * zIA * (-4 - 3 * z - z**2 + np.exp(z) * (4 - z)) / z**2
    FS1 = np.real(FS1 / M)
    FS2 = np.real(FS2 / M)
    FS3 = np.real(FS3 / M)
    QS = np.real(QS / M)

    def muaO(y):
        return y / (kaO2 + y)

    def munO(y):
        return y / (knO2 + y)

    def mua(y, z):
        return muamax * y / (kaNH4 + y) * muaO(z)

    def mun(y, z):
        return munmax * y / (knNO2 + y) * munO(z)

    def Rxa(y, z):
        return mua(y, z) - (1 + eta) * ba * muaO(z)

    def Rxn(y, z):
        return mun(y, z) - (1 + eta) * bn * munO(z)

    def mui(xa, xn, y):
        return (fxi + eta) * (ba * xa * muaO(y) + bn * xn * munO(y))

    def Rf(C, F):
        return (
            F[:, 0] * Rxa(C[:, 0], C[:, 2])
            + F[:, 1] * Rxn(C[:, 1], C[:, 2])
            + mui(F[:, 0], F[:, 1], C[:, 2])
        )

    def v(C, _S, F):
        return Int @ Rf(C, F)

    def Lp(C, S, F):
        vel = v(C, S, F)
        return S[5] / 2.0 * vel[0] + alpha / (A * rho) * (S[2] + S[3] + S[4]) - E * S[5] ** 2

    def RDE(C, S, F):
        dz = S[5] ** 2 / 4.0
        C_in = C[1:, :]
        F_in = F[1:, :]
        RNH4a = (1.0 / ya + ia) * mua(C_in[:, 0], C_in[:, 2])
        RNH4b = (ia - ii * fxi) * ba * muaO(C_in[:, 2])
        RNH4n = ia * mun(C_in[:, 1], C_in[:, 2])
        RNH4m = (ia - ii * fxi) * bn * munO(C_in[:, 2])
        RDENH4 = (dz / dNH4) * (
            (RNH4a - RNH4b) * F_in[:, 0] * rho + (RNH4n - RNH4m) * F_in[:, 1] * rho
        )
        RNO2a = 1.0 / ya * mua(C_in[:, 0], C_in[:, 2])
        RNO2n = 1.0 / yn * mun(C_in[:, 1], C_in[:, 2])
        RDENO2 = (dz / dNO2) * (RNO2n * F_in[:, 1] * rho - RNO2a * F_in[:, 0] * rho)
        RO2a = (3.43 - ya) / ya * mua(C_in[:, 0], C_in[:, 2])
        RO2b = (1.0 - fxi) * ba * muaO(C_in[:, 2])
        RO2n = (1.14 - yn) / yn * mun(C_in[:, 1], C_in[:, 2])
        RO2m = (1.0 - fxi) * bn * munO(C_in[:, 2])
        RDEO2 = (dz / dO2) * (
            (RO2a + RO2b) * F_in[:, 0] * rho + (RO2n + RO2m) * F_in[:, 1] * rho
        )
        return np.column_stack([RDENH4, RDENO2, RDEO2])

    def biofilmbvp(C, S, F):
        change = 1.0
        maa = 0
        ell = 0
        droptol = 1e-12
        mmax = 1
        GNH4 = np.zeros((N + 1, 0))
        GNO2 = np.zeros((N + 1, 0))
        GO2 = np.zeros((N + 1, 0))
        QNH4 = np.zeros((N + 1, 0))
        QNO2 = np.zeros((N + 1, 0))
        QO2 = np.zeros((N + 1, 0))
        RNH4 = np.zeros((0, 0))
        RNO2 = np.zeros((0, 0))
        RO2 = np.zeros((0, 0))
        FNH4old = FNO2old = FO2old = None
        GNH4old = GNO2old = GO2old = None
        deltaFNH4 = deltaFNO2 = deltaFO2 = None

        while change > tol:
            ell += 1
            Gtemp = RDE(C, S, F)
            # NOTE: a "corrected" Neumann-BC variant (zeroing Gtemp[-1,:] here)
            # was tested and is documented in SUMMARY.md Section 7 -- it improves
            # thickness timing around t=20d but WORSENS the Fig 2b substrate
            # profile match (O2 at the substratum: 0.06 here vs 0.23 with that
            # "fix"), so it is NOT applied by default. See SUMMARY.md before
            # re-enabling it.
            GNH4new = np.concatenate([[0.0], D2inv @ Gtemp[:, 0]]) + S[0]
            GNO2new = np.concatenate([[0.0], D2inv @ Gtemp[:, 1]]) + S[1]
            GO2new = np.concatenate([[0.0], D2inv @ Gtemp[:, 2]]) + SO2
            Gnew = np.column_stack([GNH4new, GNO2new, GO2new])
            Fnew = Gnew - C
            FNH4new = GNH4new - C[:, 0]
            FNO2new = GNO2new - C[:, 1]
            FO2new = GO2new - C[:, 2]
            if np.max(np.abs(Fnew)) < tol:
                return Gnew
            if ell > 1:
                deltaFNH4 = FNH4new - FNH4old
                deltaGNH4 = GNH4new - GNH4old
                deltaFNO2 = FNO2new - FNO2old
                deltaGNO2 = GNO2new - GNO2old
                deltaFO2 = FO2new - FO2old
                deltaGO2 = GO2new - GO2old
                if maa < mmax:
                    GNH4 = np.column_stack([GNH4, deltaGNH4]) if GNH4.size else deltaGNH4[:, None]
                    GNO2 = np.column_stack([GNO2, deltaGNO2]) if GNO2.size else deltaGNO2[:, None]
                    GO2 = np.column_stack([GO2, deltaGO2]) if GO2.size else deltaGO2[:, None]
                else:
                    GNH4 = np.column_stack([GNH4[:, 1:maa], deltaGNH4[:, None]]) if maa > 1 else deltaGNH4[:, None]
                    GNO2 = np.column_stack([GNO2[:, 1:maa], deltaGNO2[:, None]]) if maa > 1 else deltaGNO2[:, None]
                    GO2 = np.column_stack([GO2[:, 1:maa], deltaGO2[:, None]]) if maa > 1 else deltaGO2[:, None]
                maa += 1
            FNH4old, FNO2old, FO2old = FNH4new, FNO2new, FO2new
            GNH4old, GNO2old, GO2old = GNH4new, GNO2new, GO2new
            if maa == 0:
                Cnew = Gnew
            else:
                if maa == 1:
                    n1 = norm(deltaFNH4)
                    n2 = norm(deltaFNO2)
                    n3 = norm(deltaFO2)
                    QNH4 = (deltaFNH4 / n1)[:, None]
                    RNH4 = np.array([[n1]])
                    QNO2 = (deltaFNO2 / n2)[:, None]
                    RNO2 = np.array([[n2]])
                    QO2 = (deltaFO2 / n3)[:, None]
                    RO2 = np.array([[n3]])
                else:
                    if maa > mmax:
                        QNH4, RNH4 = _qrdelete(QNH4, RNH4)
                        QNO2, RNO2 = _qrdelete(QNO2, RNO2)
                        QO2, RO2 = _qrdelete(QO2, RO2)
                        maa -= 1
                    for kk in range(maa - 1):
                        RNH4[kk, maa - 1] = QNH4[:, kk] @ deltaFNH4
                        RNO2[kk, maa - 1] = QNO2[:, kk] @ deltaFNO2
                        RO2[kk, maa - 1] = QO2[:, kk] @ deltaFO2
                        deltaFNH4 = deltaFNH4 - RNH4[kk, maa - 1] * QNH4[:, kk]
                        deltaFNO2 = deltaFNO2 - RNO2[kk, maa - 1] * QNO2[:, kk]
                        deltaFO2 = deltaFO2 - RO2[kk, maa - 1] * QO2[:, kk]
                    QNH4 = np.column_stack([QNH4, deltaFNH4 / norm(deltaFNH4)])
                    # Grow R if needed
                    RNH4_new = np.zeros((maa, maa))
                    RNH4_new[: RNH4.shape[0], : RNH4.shape[1]] = RNH4
                    RNH4_new[maa - 1, maa - 1] = norm(deltaFNH4)
                    RNH4 = RNH4_new
                    QNO2 = np.column_stack([QNO2, deltaFNO2 / norm(deltaFNO2)])
                    RNO2_new = np.zeros((maa, maa))
                    RNO2_new[: RNO2.shape[0], : RNO2.shape[1]] = RNO2
                    RNO2_new[maa - 1, maa - 1] = norm(deltaFNO2)
                    RNO2 = RNO2_new
                    QO2 = np.column_stack([QO2, deltaFO2 / norm(deltaFO2)])
                    RO2_new = np.zeros((maa, maa))
                    RO2_new[: RO2.shape[0], : RO2.shape[1]] = RO2
                    RO2_new[maa - 1, maa - 1] = norm(deltaFO2)
                    RO2 = RO2_new
                Condt = np.array([_rcond(RNH4), _rcond(RNO2), _rcond(RO2)])
                while maa > 1 and np.any(Condt < droptol):
                    QNH4, RNH4 = _qrdelete(QNH4, RNH4)
                    QNO2, RNO2 = _qrdelete(QNO2, RNO2)
                    QO2, RO2 = _qrdelete(QO2, RO2)
                    GNH4 = GNH4[:, 1:maa]
                    GNO2 = GNO2[:, 1:maa]
                    GO2 = GO2[:, 1:maa]
                    maa -= 1
                    Condt = np.array([_rcond(RNH4), _rcond(RNO2), _rcond(RO2)])
                gammaNH4 = solve(RNH4, QNH4.T @ FNH4new)
                gammaNO2 = solve(RNO2, QNO2.T @ FNO2new)
                gammaO2 = solve(RO2, QO2.T @ FO2new)
                Cnew = np.column_stack(
                    [
                        GNH4new - GNH4 @ gammaNH4,
                        GNO2new - GNO2 @ gammaNO2,
                        GO2new - GO2 @ gammaO2,
                    ]
                )
            change = norm(C - Cnew)
            C = Cnew
        return C

    def transport(C, S, F):
        vel = v(C, S, F)
        vp = Rf(C, F)
        dz = ((x + 1) * (Lp(C, S, F) / S[5]))[:, None] * Dx
        fat = F[:, 0] * Rxa(C[:, 0], C[:, 2]) - F[:, 0] * vp - dz @ F[:, 0] - vel * (Dx @ F[:, 0])
        fnt = F[:, 1] * Rxn(C[:, 1], C[:, 2]) - F[:, 1] * vp - dz @ F[:, 1] - vel * (Dx @ F[:, 1])
        fit = 0.0 - fat - fnt
        return np.column_stack([fat, fnt, fit])

    def planktonic(C, S, F):
        JNH4 = dNH4 * 2.0 / S[5] * (D @ C[:, 0])
        JNO2 = dNO2 * 2.0 / S[5] * (D @ C[:, 1])
        RNH4a = (1.0 / ya + ia) * mua(S[0], SO2)
        RNH4b = (ia - ii * fxi) * muaO(SO2)
        RNH4n = ia * mun(S[1], SO2)
        RNH4m = (ia - ii * fxi) * munO(SO2)
        return np.array(
            [
                d * (SNH4in - S[0])
                - 1.0 / V * ((RNH4a - RNH4b) * S[2] + (RNH4n - RNH4m) * S[3])
                - 1.0 / V * A * JNH4[0],
                d * (SNO2in - S[1])
                - 1.0 / V * 1.0 / yn * mun(S[1], SO2) * S[3]
                + 1.0 / (V * ya) * mua(S[0], SO2) * S[2]
                - 1.0 / V * A * JNO2[0],
                S[2] * (Rxa(S[0], SO2) - d - alpha) + A * rho * F[0, 0] * E * S[5] ** 2,
                S[3] * (Rxn(S[1], SO2) - d - alpha) + A * rho * F[0, 1] * E * S[5] ** 2,
                mui(S[2], S[3], SO2) - S[4] * (d + alpha) + A * rho * F[0, 2] * E * S[5] ** 2,
                Lp(C, S, F),
            ]
        )

    SO2 = SO2in
    SNO2 = SNO2in
    SNH4 = SNH4in
    Xa = 2e-5
    Xn = 2e-5
    Xi = 0.0
    L = 1e-5
    S = np.array([SNH4, SNO2, Xa, Xn, Xi, L], dtype=float)

    CO2 = SO2 * np.ones_like(x)
    CNO2 = SNO2 * np.ones_like(x)
    CNH4 = SNH4 * np.ones_like(x)
    C = np.column_stack([CNH4, CNO2, CO2])
    fa = 0.5 * np.ones_like(x)
    fn = 0.5 * np.ones_like(x)
    fi = np.zeros_like(x)
    F = np.column_stack([fa, fn, fi])
    C = biofilmbvp(C, S, F)

    nstore = nplots // plotgap + 1
    SS = np.zeros((6, nstore))
    SS[:, 0] = S
    tt = np.zeros(nstore)
    fafa = np.zeros((N + 1, nstore))
    fafa[:, 0] = fa
    fnfn = np.zeros((N + 1, nstore))
    fnfn[:, 0] = fn
    fifi = np.zeros((N + 1, nstore))
    fifi[:, 0] = fi
    CNH4C = np.zeros((N + 1, nstore))
    CNH4C[:, 0] = C[:, 0]
    CNO2C = np.zeros((N + 1, nstore))
    CNO2C[:, 0] = C[:, 1]
    CO2C = np.zeros((N + 1, nstore))
    CO2C[:, 0] = C[:, 2]

    def timestep(C, S, F):
        NSu = planktonic(C, S, F)
        Nu = transport(C, S, F)
        au = F + QS * Nu
        as_ = S + QS * NSu
        NSa = planktonic(C, as_, au)
        Na = transport(C, as_, au)
        bu = F + QS * Na
        bs = S + QS * NSa
        NSb = planktonic(C, bs, bu)
        Nb = transport(C, bs, bu)
        cu = au + QS * (2 * Nb - Nu)
        cs = as_ + QS * (2 * NSb - NSu)
        NSc = planktonic(C, cs, cu)
        Nc = transport(C, cs, cu)
        S = S + FS1 * NSu + 2 * FS2 * (NSa + NSb) + FS3 * NSc
        F = F + FS1 * Nu + 2 * FS2 * (Na + Nb) + FS3 * Nc
        return S, F

    for i in range(1, nplots + 1):
        t = t + dt
        S, F = timestep(C, S, F)
        C = biofilmbvp(C, S, F)
        if np.isnan(C[-1, 0]):
            return tt, SS, CNH4C, CNO2C, CO2C, fafa, fnfn, fifi, (1 + x) / 2.0
        if i % plotgap == 0:
            j = i // plotgap
            tt[j] = t
            SS[:, j] = S
            fafa[:, j] = F[:, 0]
            fnfn[:, j] = F[:, 1]
            fifi[:, j] = F[:, 2]
            CNH4C[:, j] = C[:, 0]
            CNO2C[:, j] = C[:, 1]
            CO2C[:, j] = C[:, 2]
            if j % max(1, points // 200) == 0:
                print("Percent done: {:7.3f}".format(100.0 * i / nplots))

    print("Percent done: 100")
    x = (1 + x) / 2.0
    return tt, SS, CNH4C, CNO2C, CO2C, fafa, fnfn, fifi, x
