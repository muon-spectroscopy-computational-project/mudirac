/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * integrate.hpp
 * 
 * Routines for integrating differential equations with the shooting method
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "integrate.hpp"

/**
 * @brief  Integrate two coupled ODEs
 * @note   Integrate a system of coupled equations of the form:
 * 
 *      Q' = AA*Q+AB*P
 *      P' = BA*Q+BB*P
 * 
 * with a shooting method, up to a given index, either forward or backwards.
 * 
 * @param  &Q: Vector for Q. Will return the integrated values, must contain already the first two as boundary conditions.
 * @param  &P: Vector for P. Will return the integrated values, must contain already the first two as boundary conditions.
 * @param  AA: Vector for AA (see definition above). Same size as Q and P.
 * @param  AB: Vector for AB.
 * @param  BA: Vector for BA.
 * @param  BB: Vector for BB.
 * @param  h:  Integration step (default = 1).
 * @param  stop_i: Index to stop integration. Run the whole range if -1 (default = -1).
 * @param  dir: Integration direction, either forward 'f' or backwards 'b' (default = 'f').
 * @retval None
 */
void shoot_QP(vector<double> &Q, vector<double> &P, vector<double> AA, vector<double> AB, vector<double> BA, vector<double> BB,
              double h, int stop_i, char dir)
{
    int N = Q.size();
    int step = (dir == 'f') ? 1 : -1;
    int from_i = (step == 1) ? 2 : N - 3;
    double QA, QB, QC, PA, PB, PC;

    // First, check size
    if (P.size() != N || AA.size() != N || AB.size() != N || BA.size() != N || BB.size() != N)
    {
        throw "Invalid size for one or more arrays passed to shoot_QP";
    }

    if (stop_i == -1)
    {
        stop_i = (step == 1) ? N - 1 : 0;
    }

    for (int i = from_i; step * (i - stop_i) <= 0; i += step)
    {
        QC = 1 / h * (2 * Q[i - step] - 0.5 * Q[i - 2 * step]) * step;
        PC = 1 / h * (2 * P[i - step] - 0.5 * P[i - 2 * step]) * step;
        QA = (step * 1.5 / h - AA[i]);
        QB = -AB[i];
        PA = (step * 1.5 / h - BB[i]);
        PB = -BA[i];
        Q[i] = (PC * QB - PA * QC) / (PB * QB - PA * QA);
        P[i] = (PB * QC - PC * QA) / (PB * QB - PA * QA);
    }

    return;
}