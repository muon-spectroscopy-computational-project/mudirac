#include "../lib/experiment.hpp"
#include "../lib/utils.hpp"
#include "../vendor/catch/catch.hpp"
using namespace std;


TEST_CASE("test fermi parameterisation c", "[fermiParameters]")
{
    double output_c, output_t, rms_radius;
    tie(output_c, output_t) = fermiParameters(5.88, 0.0);
    rms_radius = rmsRadius(output_c, output_t);
    REQUIRE( output_c == Approx(7.59).epsilon(0.01) );
    REQUIRE(output_t == Approx(0.0).epsilon(0.01));
    REQUIRE ( rms_radius == Approx(5.88).epsilon(0.01));
}

TEST_CASE("test fermi parameterisation 0 ", "[fermiParameters]")
{
    double output_c, output_t, rms_radius;
    tie(output_c, output_t) = fermiParameters(0, 0);
    rms_radius = rmsRadius(output_c, output_t);
    REQUIRE( output_c == Approx(0).epsilon(0.01));
    REQUIRE(output_t == Approx(0.0).epsilon(0.01));
    REQUIRE(rms_radius == Approx(0.0).margin(0.00001));
}

TEST_CASE("test fermi parameterisation t", "[fermiParameters]")
{
    double output_c, output_t, rms_radius;
    tie(output_c, output_t) = fermiParameters(5.88, M_PI/2);
    rms_radius = rmsRadius(output_c, output_t);
    REQUIRE( output_c == Approx(0.0).margin(0.00001));
    REQUIRE(output_t == Approx(6.95).epsilon(0.01));
    REQUIRE(rms_radius == Approx(5.88).margin(0.00001));
}

TEST_CASE("test fermi parameterisation normal use", "[fermiParameters]")
{
    double output_c, output_t, rms_radius;
    tie(output_c, output_t) = fermiParameters(5.88, M_PI/12);
    rms_radius = rmsRadius(output_c, output_t);
    REQUIRE( output_c == Approx(7.285559152).margin(0.00001));
    REQUIRE(output_t == Approx(1.952159691).margin(0.00001));
    REQUIRE(rms_radius == Approx(5.88).margin(0.00001));
}
