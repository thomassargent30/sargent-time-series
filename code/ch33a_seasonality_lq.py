"""
Shared module for Chapter 33a — Hansen & Sargent (1993), "Seasonality and
approximation errors in rational expectations models," J. Econometrics 55, 21-55.

It solves the linear-quadratic planning model of Section 3 and returns the
equilibrium spectral densities of the observables (consumption, investment,
endowment), used to reconstruct the paper's Figures 4-8.

Model (mean-zero deviations; means are irrelevant for spectral densities):

    felicity   -1/2 [ g_t^2 + phi1^2 i_t^2 ],   g_t = c_t - h_t - b_t
    c_t = gamma k_{t-1} + d_t - i_t
    k_t = deltak k_{t-1} + i_t
    h_t = lambda (1-deltah) H_t,   H_t = c_{t-p} + deltah H_{t-p}    (seasonal habit)
    A_b(L) b_t = w1_t,   A_d(L) d_t = w2_t                            (AR shocks)

The planner chooses i_t.  The problem is a discounted stochastic optimal linear
regulator, solved here by Riccati iteration (numpy only).  The solver is
validated in `_selftest()` below by (i) the analytic Euler-equation orthogonality
condition for the no-habit case (residual ~1e-12) and (ii) equality of the
spectral-density area with the stationary variance.

Run `python ch33a_seasonality_lq.py` to execute the self-tests.
"""
import numpy as np


# --------------------------------------------------------------------------- #
#  Linear-quadratic machinery                                                  #
# --------------------------------------------------------------------------- #
def _ar_companion(alpha):
    """Companion form for x_t = sum_k alpha_k x_{t-k} + w_t."""
    n = len(alpha)
    M = np.zeros((n, n))
    M[0, :] = alpha
    if n > 1:
        M[1:, :-1] = np.eye(n - 1)
    e = np.zeros((n, 1))
    e[0, 0] = 1.0
    return M, e


def riccati(A, B, R, Q, W, beta, iters=50000, tol=1e-14):
    """Maximize -E sum beta^t (x'Rx + 2 x'Wu + u'Qu) s.t. x_{t+1}=Ax_t+Bu_t.

    Returns (P, F) with optimal policy u_t = -F x_t.
    """
    P = R.copy()
    for _ in range(iters):
        BPB = Q + beta * B.T @ P @ B
        F = np.linalg.solve(BPB, beta * B.T @ P @ A + W.T)
        Pn = R + beta * A.T @ P @ A - (beta * A.T @ P @ B + W) @ F
        Pn = 0.5 * (Pn + Pn.T)
        if np.max(np.abs(Pn - P)) < tol:
            P = Pn
            break
        P = Pn
    BPB = Q + beta * B.T @ P @ B
    F = np.linalg.solve(BPB, beta * B.T @ P @ A + W.T)
    return P, F


def build_model(gamma, phi1, deltak, beta, alpha_b, sig_b, alpha_d, sig_d,
                lam=0.0, deltah=0.9, p=4):
    """Assemble the state-space (A, B, C) and cost (R, Q, W).

    State x_t = [ k_{t-1}, b-block(t), d-block(t), (c-lags, H-lags if habit) ];
    control u_t = i_t; shocks enter x_{t+1} through C.
    """
    Mb, eb = _ar_companion(alpha_b)
    Md, ed = _ar_companion(alpha_d)
    nb, nd = len(alpha_b), len(alpha_d)
    habit = lam != 0.0
    nc = p if habit else 0
    nH = p if habit else 0
    n = 1 + nb + nd + nc + nH
    ik, ib, idd, ic, iH = 0, 1, 1 + nb, 1 + nb + nd, 1 + nb + nd + (p if habit else 0)

    A = np.zeros((n, n))
    B = np.zeros((n, 1))
    C = np.zeros((n, 2))
    A[ik, ik] = deltak
    B[ik, 0] = 1.0                                   # k_t = deltak k_{t-1} + i_t
    A[ib:ib + nb, ib:ib + nb] = Mb
    C[ib:ib + nb, 0:1] = eb * sig_b
    A[idd:idd + nd, idd:idd + nd] = Md
    C[idd:idd + nd, 1:2] = ed * sig_d

    lh = lam * (1.0 - deltah)
    eg = np.zeros((n, 1))
    eg[ik, 0] = gamma
    eg[idd, 0] = 1.0
    eg[ib, 0] = -1.0
    if habit:
        eg[ic + (p - 1), 0] += -lh                   # c_{t-p}
        eg[iH + (p - 1), 0] += -lh * deltah          # H_{t-p}
        # c-lag register: new c_t = gamma k + d0 - i_t
        A[ic, ik] = gamma
        A[ic, idd] = 1.0
        B[ic, 0] = -1.0
        for j in range(1, p):
            A[ic + j, ic + j - 1] = 1.0
        # H-lag register: new H_t = c_{t-p} + deltah H_{t-p}
        A[iH, ic + (p - 1)] = 1.0
        A[iH, iH + (p - 1)] = deltah
        for j in range(1, p):
            A[iH + j, iH + j - 1] = 1.0

    R = 0.5 * (eg @ eg.T)
    W = -0.5 * eg
    Q = np.array([[0.5 * (1.0 + phi1 ** 2)]])
    return dict(n=n, ik=ik, idd=idd, eg=eg, A=A, B=B, C=C, R=R, Q=Q, W=W,
                gamma=gamma, phi1=phi1, deltak=deltak, beta=beta)


def solve(info):
    """Solve the regulator and return closed loop plus observation rows."""
    A, B, C, R, Q, W, beta = (info[k] for k in "ABCRQW") if False else (
        info['A'], info['B'], info['C'], info['R'], info['Q'], info['W'], info['beta'])
    P, F = riccati(A, B, R, Q, W, beta)
    Acl = A - B @ F                                   # x_{t+1} = Acl x_t + C w_{t+1}
    n, ik, idd = info['n'], info['ik'], info['idd']
    u = -F                                            # i_t = u x_t
    c = np.zeros((1, n)); c[0, ik] = info['gamma']; c[0, idd] = 1.0
    c = c - u                                         # c_t = gamma k + d - i
    d = np.zeros((1, n)); d[0, idd] = 1.0
    return dict(P=P, F=F, Acl=Acl, C=C, c=c, i=u.copy(), d=d,
                g=info['eg'].T - u, info=info)


def spectrum(Acl, C, rows, omega):
    """Spectral density matrix of `rows @ x_t` at frequencies `omega`.

    S(w) = H(w) H(w)^*,  H(w) = rows (I - Acl e^{-i w})^{-1} C.
    Returns array of shape (len(omega), r, r); diagonals are real densities.
    """
    n = Acl.shape[0]
    I = np.eye(n)
    r = rows.shape[0]
    out = np.zeros((len(omega), r, r), dtype=complex)
    for m, w in enumerate(omega):
        Tr = np.linalg.solve(I - Acl * np.exp(-1j * w), C)
        H = rows @ Tr
        out[m] = H @ H.conj().T
    return out


def lyap(Acl, Q):
    """Exact stationary covariance V = Acl V Acl' + Q."""
    n = Acl.shape[0]
    return np.linalg.solve(np.eye(n * n) - np.kron(Acl, Acl),
                           Q.reshape(-1)).reshape(n, n)


# quarterly seasonal frequencies (period p = 4): omega = pi/2 (period 4), pi (period 2)
SEASONAL_OMEGA = [np.pi / 2, np.pi]


# --------------------------------------------------------------------------- #
#  Self-tests                                                                  #
# --------------------------------------------------------------------------- #
def _selftest():
    # (1) no-habit: analytic Euler-equation orthogonality residual should be ~0
    info = build_model(0.1, 0.3, 0.95, 1 / 1.05, [0.2], 0.25, [0.4], 1.0)
    s = solve(info)
    beta, gamma, deltak, phi1 = 1 / 1.05, 0.1, 0.95, 0.3
    V = lyap(s['Acl'], s['C'] @ s['C'].T)
    Etg1 = s['g'] @ s['Acl']
    Eti1 = s['i'] @ s['Acl']
    eta = (s['g'] - phi1 ** 2 * s['i']
           - beta * (gamma + deltak) * Etg1 + beta * deltak * phi1 ** 2 * Eti1)
    euler = np.max(np.abs(eta.flatten() * np.sqrt(np.diag(V))))
    # (2) spectrum area equals variance
    w = np.linspace(-np.pi, np.pi, 20000, endpoint=False)
    area = np.real(spectrum(s['Acl'], s['C'], s['c'], w)[:, 0, 0]).mean()
    var = float(s['c'] @ V @ s['c'].T)
    print("self-test  Euler residual   = %.2e  (expect ~1e-12)" % euler)
    print("self-test  spectrum area(c) = %.6f   Var(c) = %.6f" % (area, var))
    assert euler < 1e-8 and abs(area / var - 1) < 1e-4
    print("OK")


if __name__ == "__main__":
    _selftest()
