"""Classical fourth-order Runge–Kutta step."""


def rk4(f, x, h):
    k1 = f(x)
    k2 = f(x + 0.5 * h * k1)
    k3 = f(x + 0.5 * h * k2)
    k4 = f(x + h * k3)
    return x + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
