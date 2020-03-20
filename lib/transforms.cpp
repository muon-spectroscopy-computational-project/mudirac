/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 * 
 * transforms.hpp
 * 
 * Functions to perform various kinds of function transforms
 * 
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "transforms.hpp"

/**
 * @brief Perform a Discrete Cosine Transform of type IV
 * @note  Perform a Discrete Cosine Transform of type IV on
 * the function passed. The transform is defined as:
 * 
 * X_k = sum_{n=0}^{N-1} x_n cos[pi/N*(n+1/2)*(k+1/2)]
 * 
 * @param f:                Function to transform
 * @return vector<double>   Transformed function
 */
vector<double> dctIV(vector<double> f)
{
    int N = f.size();
    double w = M_PI / N;
    vector<double> cf(N, 0.0);

    for (int k = 0; k < N; ++k)
    {
        for (int n = 0; n < N; ++n)
        {
            cf[k] += cos(w * (n + 0.5) * (k + 0.5)) * f[n];
        }
    }

    return cf;
}

/**
 * @brief Perform an inverse Discrete Cosine Transform of type IV
 * @note  Perform an inverse Discrete Cosine Transform of type IV on
 * the function passed. This is basically the same as a DCT-IV with
 * a prefactor of 2/N.
 * 
 * @param f:                Function to transform
 * @return vector<double>   Transformed function
 */
vector<double> invDctIV(vector<double> f)
{
    int N = f.size();
    vector<double> cf = dctIV(f);

    for (int i = 0; i < N; ++i)
    {
        cf[i] = cf[i] * 2.0 / N;
    }

    return cf;
}