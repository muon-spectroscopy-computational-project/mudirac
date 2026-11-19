"""Reference values for the uniform-sphere Uehling potential (test_potential.cpp)

V(r) = K/r * int_0^1 sqrt(1 - u^2) (1 + u^2/2) y(u; r, R) du,   K = -2 alpha^2 rho/3,   rho = 3 Z/(4 pi R^3)

with the kernel y of UehlingSpherePotential in its cancellation-free form (a = R u alpha/2, b = (u alpha)^2/4):
    r <= R:  y = r u alpha - (exp(-2c(R-r)/u) - exp(-2c(R+r)/u)) (a + b)
    r >  R:  y = exp(-2c(r-R)/u) (a - b) + exp(-2c(r+R)/u) (a + b)
which is algebraically identical to the original form. The integral is done with mpmath at 40 digits, with
u = sin(theta) to remove the square-root endpoint, and with the same double-precision constants and radii
as the C++ test. Run with: python3 uehling_sphere_ref.py
"""
from mpmath import mp, mpf, quad, exp, sin, cos, pi

mp.dps = 40

ALPHA = mpf(7.2973525664e-3)
C = mpf(137.035999139)
ANG = 1.0 / 5.2917721067e-1
FM = ANG * 1e-5


def uehling_sphere(Z, R, r):
    Z, R, r = mpf(Z), mpf(R), mpf(r)
    rho = 3 * Z / (4 * pi * R**3)
    K = -2 * ALPHA**2 * rho / 3

    def f(th):
        u = sin(th)
        if u == 0:
            return mpf(0)
        a = R * u * ALPHA / 2
        b = (u * ALPHA) ** 2 / 4
        em = exp(-2 * C * abs(R - r) / u)
        ep = exp(-2 * C * (R + r) / u)
        if r <= R:
            y = r * u * ALPHA - (em - ep) * (a + b)
        else:
            y = em * (a - b) + ep * (a + b)
        # du = cos(theta) dtheta and sqrt(1 - u^2) = cos(theta)
        return cos(th) ** 2 * (1 + u * u / 2) * y

    return K / r * quad(f, [0, pi / 16, pi / 8, pi / 4, pi / 2])


def uehling_sphere_centre(Z, R):
    """V at r = 0: K times the integral of the r -> 0 limit of y/r,
    u alpha - (4c/u) exp(-2cR/u) (a + b)."""
    Z, R = mpf(Z), mpf(R)
    rho = 3 * Z / (4 * pi * R**3)
    K = -2 * ALPHA**2 * rho / 3

    def f(th):
        u = sin(th)
        if u == 0:
            return mpf(0)
        a = R * u * ALPHA / 2
        b = (u * ALPHA) ** 2 / 4
        g0 = u * ALPHA - 4 * C / u * exp(-2 * C * R / u) * (a + b)
        return cos(th) ** 2 * (1 + u * u / 2) * g0

    return K * quad(f, [0, pi / 16, pi / 8, pi / 4, pi / 2])


if __name__ == "__main__":
    for Z, Rfm in ((82, 7.1), (14, 4.0)):
        R = Rfm * FM
        print(f"Z = {Z}, R = {Rfm} fm, r = 0: V = {mp.nstr(uehling_sphere_centre(Z, R), 17)}")
        for x in (0.01, 0.5, 0.99, 1.0, 1.5, 4.0):
            r = x * R
            print(f"Z = {Z}, R = {Rfm} fm, r/R = {x}: V = {mp.nstr(uehling_sphere(Z, R, r), 17)}")
