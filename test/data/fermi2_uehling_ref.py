"""Reference values for the Uehling potential of the Fermi 2-term charge density (test_potential.cpp)

Computed by a route independent of UehlingFermi2Potential: the direct convolution of the density with the
point-charge Uehling potential U(d) = -(2 alpha/(3 pi d)) int_0^1 sqrt(1 - u^2)(1 + u^2/2)/u exp(-2 c d/u) du.
The angular integral and the integral over the distance d are done analytically, which gives

V(r) = -(2 alpha/(3 c r)) int_0^inf s rho(s) ds int_0^1 sqrt(1 - u^2)(1 + u^2/2) [exp(-2c|r-s|/u) - exp(-2c(r+s)/u)] du

with rho(s) = Z f(s)/N, f(s) = 1/(1 + exp(4 ln3 (s - c_F)/t)), N = int_0^inf 4 pi s^2 f(s) ds. The double integral
is done with mpmath (u = sin(theta); s split at r and around c_F), with the same double-precision constants and
2pF parameters as the C++ test. Run with: python3 fermi2_uehling_ref.py [dps]
"""
import math
import sys

from mpmath import mp, mpf, quad, exp, log, pi, sin, cos, inf

mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20

ALPHA = mpf(7.2973525664e-3)
C = mpf(137.035999139)
ANG = 1.0 / 5.2917721067e-1
FM = ANG * 1e-5


def uehling_fermi2(Z, cF, t, r):
    Z, cF, t, r = mpf(Z), mpf(cF), mpf(t), mpf(r)
    a = 4 * log(3) / t

    def f(s):
        return 1 / (1 + exp(a * (s - cF)))

    pts = [mpf(0)] + [p for p in (cF - 12 / a, cF - 3 / a, cF, cF + 3 / a, cF + 12 / a, cF + 45 / a) if p > 0]
    N = quad(lambda s: 4 * pi * s**2 * f(s), pts + [inf])

    def inner_u(s):
        d1 = abs(r - s)
        d2 = r + s

        def g(th):
            u = sin(th)
            if u == 0:
                return mpf(0)
            return cos(th) ** 2 * (1 + u * u / 2) * (exp(-2 * C * d1 / u) - exp(-2 * C * d2 / u))

        return quad(g, [0, pi / 32, pi / 8, pi / 2])

    spts = sorted(set([p for p in pts] + [r]))
    integral = quad(lambda s: s * f(s) * inner_u(s), spts + [cF + 80 / a])
    return -(2 * ALPHA / (3 * C * r)) * (Z / N) * integral


def uehling_uniform(Z, cF, r):
    """The t = 0 limit: the same convolution for a uniform sphere of radius cF."""
    Z, cF, r = mpf(Z), mpf(cF), mpf(r)
    N = 4 * pi * cF**3 / 3

    def inner_u(s):
        d1 = abs(r - s)
        d2 = r + s

        def g(th):
            u = sin(th)
            if u == 0:
                return mpf(0)
            return cos(th) ** 2 * (1 + u * u / 2) * (exp(-2 * C * d1 / u) - exp(-2 * C * d2 / u))

        return quad(g, [0, pi / 32, pi / 8, pi / 2])

    spts = sorted(set([mpf(0), cF] + ([r] if r < cF else [])))
    integral = quad(lambda s: s * inner_u(s), spts)
    return -(2 * ALPHA / (3 * C * r)) * (Z / N) * integral


if __name__ == "__main__":
    cU = 6.6 * FM
    for x in (0.01, 0.5, 1.5, 4.0):
        print(f"t = 0, Z = 82, c = 6.6 fm, r/c = {x}: V = {mp.nstr(uehling_uniform(82, cU, x * cU), 17)}", flush=True)
    t = 2.3 * FM
    for Z, Rfm in ((82, 7.102018661385413), (14, 4.031001066732679)):
        R = Rfm * FM
        cF = math.sqrt(R * R - 7.0 / 3.0 * math.pow(math.pi * t / (4 * math.log(3.0)), 2))
        print(f"Z = {Z}: c = {cF!r}", flush=True)
        for x in (0.1, 0.5, 1.0, 1.5, 3.0, 10.0):
            r = x * cF
            print(f"Z = {Z}, r/c = {x}: V = {mp.nstr(uehling_fermi2(Z, cF, t, r), 17)}", flush=True)
