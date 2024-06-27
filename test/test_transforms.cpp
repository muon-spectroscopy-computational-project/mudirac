#include <vector>
#include <iostream>
#include <cmath>
#include "../lib/transforms.hpp"

#include "../vendor/catch/catch.hpp"
#include "macros.hpp"

#define ERRTOL_HIGH 1e-5
#define ERRTOL_LOW 1e-3

TEST_CASE("Discrete cosine transform", "[dctIV]")
{
    // Simple case
    int N = 100;
    int k = 3;
    vector<double> f(N), cf(N);

    for (int i = 0; i < N; ++i)
    {
        f[i] = cos((i + 0.5) * (k + 0.5) * M_PI / N);
    }

    cf = dctIV(f);

    REQUIRE(cf[k] == Approx(N / 2.0));

    // More complex

    for (int i = 0; i < N; ++i)
    {
        f[i] = exp(-pow(4.0 * i / N, 2.0));
    }

    cf = dctIV(f);
    cf = invDctIV(cf);

    CHECK_APPROX(f, cf, ERRTOL_HIGH);
}