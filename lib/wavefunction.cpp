/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * wavefunction.cpp
 *
 * Functions useful for operations on wavefunctions
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "wavefunction.hpp"

/**
  * @brief  Compute the bra-ket of two radial wavefunctions on a logarithmic grid
  * @note   Compute the bra-ket integral of two radial wavefunctions defined on a 
  * logarithmic grid:
  * 
  *     <psiL|psiR> = \int_0^\infty psiL(r)*psiR(r)*dr
  * 
  * @param  psiL:   Bra wavefunction
  * @param  psiR:   Ket wavefunction
  * @param  r:      Radial grid (space units)
  * @param  dx:     Logarithmic grid step
  * @retval 
 */
double braKetLog(vector<double> psiL, vector<double> psiR, vector<double> r, double dx)
{
    int N = psiL.size();
    vector<double> ker(N);

    if (psiR.size() != N)
    {
        throw invalid_argument("Bra and ket dimensions do not match");
    }

    for (int i = 0; i < N; ++i)
    {
        ker[i] = psiL[i] * psiR[i] * r[i];
    }

    return trapzInt(dx, ker);
}

/**
  * @brief  Compute the bra-ket of two radial wavefunctions on a logarithmic grid
  * @note   Compute the bra-ket integral of two radial wavefunctions defined on a 
  * logarithmic grid:
  * 
  *     <psiL|psiR> = \int_0^\infty psiL(r)*psiR(r)*dr
  * 
  * This version is designed for Dirac spinors, accepting two-component wavefunctions.
  * 
  * @param  psiL:   Bra wavefunction
  * @param  psiR:   Ket wavefunction
  * @param  r:      Radial grid (space units)
  * @param  dx:     Logarithmic grid step
  * @retval 
 */
double braKetLog(vector<vector<double>> psiL, vector<vector<double>> psiR, vector<double> r, double dx)
{
    int N = psiL[0].size();
    vector<double> ker(N);

    if (psiL.size() != 2 || psiR.size() != 2 || psiR[0].size() != N)
    {
        throw invalid_argument("Bra and ket dimensions do not match");
    }

    for (int i = 0; i < N; ++i)
    {
        ker[i] = (psiL[0][i] * psiR[0][i] + psiL[1][i] * psiR[1][i]) * r[i];
    }

    return trapzInt(dx, ker);
}

/**
  * @brief  Compute the matrix element of two radial wavefunctions on a logarithmic grid
  * @note   Compute the matrix element of two radial wavefunctions defined on a 
  * logarithmic grid with a given operator:
  * 
  *     <psiL|op|psiR> = \int_0^\infty psiL(r)*psiR(r)*op(r)*dr
  * 
  * This version assumes that the operator is diagonal in the space basis, like V or r.
  * 
  * @param  psiL:   Bra wavefunction
  * @param  op:     Operator
  * @param  psiR:   Ket wavefunction
  * @param  r:      Radial grid (space units)
  * @param  dx:     Logarithmic grid step
  * @retval 
 */
double braOpKetLog(vector<double> psiL, vector<double> op, vector<double> psiR, vector<double> r, double dx)
{
    int N = psiL.size();
    vector<double> ker(N);

    if (psiR.size() != N || op.size() != N)
    {
        throw invalid_argument("Bra, operator and ket dimensions do not match");
    }

    for (int i = 0; i < N; ++i)
    {
        ker[i] = psiL[i] * psiR[i] * op[i] * r[i];
    }

    return trapzInt(dx, ker);
}

/**
  * @brief  Compute the matrix element of two radial wavefunctions on a logarithmic grid
  * @note   Compute the matrix element of two radial wavefunctions defined on a 
  * logarithmic grid with a given operator:
  * 
  *     <psiL|op|psiR> = \int_0^\infty psiL(r)*psiR(r)*op(r)*dr
  * 
  * This version accepts a function as operator. The function must take the ket and operate on
  * it, returning the result, for example op(psi) = d psi / dx.
  * 
  * @param  psiL:   Bra wavefunction
  * @param  op:     Operator
  * @param  psiR:   Ket wavefunction
  * @param  r:      Radial grid (space units)
  * @param  dx:     Logarithmic grid step
  * @retval 
 */
double braOpKetLog(vector<double> psiL, vector<double> (*op)(vector<double> psi), vector<double> psiR, vector<double> r, double dx)
{
    vector<double> opval = op(psiR);
    return braKetLog(psiL, opval, r, dx);
}

/**
  * @brief  Compute the matrix element of two radial wavefunctions on a logarithmic grid
  * @note   Compute the matrix element of two radial wavefunctions defined on a 
  * logarithmic grid with a given operator:
  * 
  *     <psiL|op|psiR> = \int_0^\infty psiL(r)*psiR(r)*op(r)*dr
  * 
  * This version assumes that the operator is diagonal in the space basis, like V or r.
  * This version is designed for Dirac spinors, accepting two-component wavefunctions.
  * 
  * @param  psiL:   Bra wavefunction
  * @param  op:     Operator
  * @param  psiR:   Ket wavefunction
  * @param  r:      Radial grid (space units)
  * @param  dx:     Logarithmic grid step
  * @retval 
 */
double braOpKetLog(vector<vector<double>> psiL, vector<double> op, vector<vector<double>> psiR, vector<double> r, double dx)
{
    int N = psiL[0].size();
    vector<double> ker(N);

    if (psiL.size() != 2 || psiR.size() != 2 || psiR[0].size() != N || op.size() != N)
    {
        throw invalid_argument("Bra, operator and ket dimensions do not match");
    }

    for (int i = 0; i < N; ++i)
    {
        ker[i] = (psiL[0][i] * psiR[0][i] + psiL[1][i] * psiR[1][i]) * op[i] * r[i];
    }

    return trapzInt(dx, ker);
}

/**
  * @brief  Compute the matrix element of two radial wavefunctions on a logarithmic grid
  * @note   Compute the matrix element of two radial wavefunctions defined on a 
  * logarithmic grid with a given operator:
  * 
  *     <psiL|op|psiR> = \int_0^\infty psiL(r)*psiR(r)*op(r)*dr
  * 
  * This version accepts a function as operator. The function must take the ket and operate on
  * it, returning the result, for example op(psi) = d psi / dx.
  * This version is designed for Dirac spinors, accepting two-component wavefunctions.
  * 
  * @param  psiL:   Bra wavefunction
  * @param  op:     Operator
  * @param  psiR:   Ket wavefunction
  * @param  r:      Radial grid (space units)
  * @param  dx:     Logarithmic grid step
  * @retval 
 */
double braOpKetLog(vector<vector<double>> psiL, vector<vector<double>> (*op)(vector<vector<double>>), vector<vector<double>> psiR, vector<double> r, double dx)
{
    vector<vector<double>> opval = op(psiR);
    return braKetLog(psiL, opval, r, dx);
}