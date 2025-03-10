#include "../lib/experiment.hpp"
#include "../lib/utils.hpp"
#include "../vendor/catch/catch.hpp"
using namespace std;

TEST_CASE("test fermi parameterisation c", "[fermiParameters]")
{
    double output_c, output_t;
    tie(output_c, output_t) = fermiParameters(5.88, 0);
    REQUIRE( output_c == Approx(7.59).epsilon(0.01) );
    REQUIRE(output_t == Approx(0.0).epsilon(0.01));
}

TEST_CASE("test fermi parameterisation 0 ", "[fermiParameters]")
{
    double output_c, output_t;
    tie(output_c, output_t) = fermiParameters(0, 0);
    REQUIRE( output_c == Approx(0).epsilon(0.01) );
    REQUIRE(output_t == Approx(0.0).epsilon(0.01));
}

TEST_CASE("test fermi parameterisation t", "[fermiParameters]")
{
    double output_c, output_t;
    tie(output_c, output_t) = fermiParameters(5.88, M_PI/2);
    REQUIRE( output_c == Approx(0.0).margin(0.001));
    REQUIRE(output_t == Approx(6.95).epsilon(0.01));
}

