// Main file for catch.hpp
#define CATCH_CONFIG_MAIN
#include "../vendor/catch/catch.hpp"

#define CHECK_APPROX(v1, v2)                               \
    {                                                      \
        bool ans = false;                                  \
        if ((v1).size() == (v2).size())                    \
        {                                                  \
            for (int i = 0; i < (v1).size(); ++i)          \
            {                                              \
                ans = ans && ((v1)[i] == Approx((v2)[i])); \
            }                                              \
        }                                                  \
        REQUIRE(ans);                                      \
    }
