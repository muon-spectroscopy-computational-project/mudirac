"""Reference values for the Fermi 2-term Coulomb potential (test_potential.cpp)

For the charge density rho(r) ~ f(r) = 1/(1 + exp(4 ln3 (r - c)/t)), normalised to Z,

V(r) = -(Z/N) [ (1/r) int_0^r 4 pi s^2 f(s) ds + int_r^inf 4 pi s f(s) ds ],   N = int_0^inf 4 pi s^2 f(s) ds,

integrated with mpmath at 40 digits, for the parameters of the C++ test (Pb-208, R = 7.102018661385413 fm,
t = 2.3 fm, c from MuDirac's rule c = sqrt(R^2 - 7/3 (pi t/(4 ln3))^2)), with the same double-precision
constants. Run with: python3 fermi2_coulomb_ref.py
"""
import math

from mpmath import mp, mpf, quad, exp, log, pi, inf

mp.dps = 40

ANG = 1.0 / 5.2917721067e-1
FM = ANG * 1e-5


def fermi2_coulomb(Z, c, t, r):
    Z, c, t, r = mpf(Z), mpf(c), mpf(t), mpf(r)
    a = 4 * log(3) / t

    def f(s):
        return 1 / (1 + exp(a * (s - c)))

    pts = sorted(set([mpf(0), c - 10 / a, c, c + 10 / a, c + 40 / a]))
    pts = [p for p in pts if p >= 0]
    N = quad(lambda s: 4 * pi * s**2 * f(s), pts + [inf])
    inner = [p for p in pts if p < r] + [r]
    outer = [r] + [p for p in pts if p > r] + [inf]
    q_in = quad(lambda s: 4 * pi * s**2 * f(s), inner)
    q_out = quad(lambda s: 4 * pi * s * f(s), outer)
    return -(Z / N) * (q_in / r + q_out)


if __name__ == "__main__":
    R = 7.102018661385413 * FM
    t = 2.3 * FM
    # c exactly as CoulombFermi2Potential computes it in double precision
    c = math.sqrt(R * R - 7.0 / 3.0 * math.pow(math.pi * t / (4 * math.log(3.0)), 2))
    print(f"c = {c!r}")
    for x in (0.1, 0.5, 1.0, 1.5, 3.0):
        r = x * c
        print(f"r/c = {x}: V = {mp.nstr(fermi2_coulomb(82, c, t, r), 17)}")
