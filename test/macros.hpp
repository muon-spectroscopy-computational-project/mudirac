#define CHECK_APPROX(v1, v2, eps)                                         \
    {                                                                     \
        bool ans = ((v1).size() == (v2).size());                          \
        if (ans)                                                          \
        {                                                                 \
            for (int i = 0; i < (v1).size(); ++i)                         \
            {                                                             \
                ans = ans && ((v1)[i] == Approx((v2)[i]).epsilon((eps))); \
            }                                                             \
        }                                                                 \
        REQUIRE(ans);                                                     \
    }
