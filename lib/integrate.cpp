/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 * 
 * integrate.cpp
 * 
 * Routines for integrating differential equations with the shooting method
 * 
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include "integrate.hpp"

// Exceptions
TurningPointError::TurningPointError(TPEType t)
{

    type = t;

    switch (type)
    {
    case RMIN_BIG:
        msg = "Range for integration does not include turning point (r_min > r_turn)";
        break;
    case RMAX_SMALL:
        msg = "Range for integration does not include turning point (r_max < r_turn)";
        break;
    default:
        msg = "Unknown turning point error";
        break;
    }
}

/**
 * @brief  Integrate a function with the trapezoidal rule
 * @note   Perform a trapezoidal rule integration of function y over range x:
 * 
 * I = sum_i (y[i]+y[i-1])/2*(x[i]-x[i-1])
 * 
 * @param  x: x axis values
 * @param  y: function values
 * @retval Integral
 */
double trapzInt(vector<double> x, vector<double> y)
{
    int N = y.size();
    double ans = 0.0;

    if (x.size() != N)
    {
        throw invalid_argument("Invalid size for arrays passed to trapzInt");
    }

    for (int i = 1; i < N; ++i)
    {
        ans += (y[i] + y[i - 1]) / 2.0 * (x[i] - x[i - 1]);
    }

    return ans;
}

/**
 * @brief  Integrate a function with the trapezoidal rule
 * @note   Perform a trapezoidal rule integration of function y with step dx:
 * 
 * I = sum_i (y[i]+y[i-1])/2*dx
 * 
 * @param  dx:  x step
 * @param  y:   function values
 * @retval Integral
 */
double trapzInt(double dx, vector<double> y)
{
    int N = y.size();
    double ans = 0.0;

    for (int i = 1; i < N; ++i)
    {
        ans += (y[i] + y[i - 1]) / 2.0 * dx;
    }

    return ans;
}

/**
 * @brief A single Runge-Kutta step
 * @note  Perform a single Runge-Kutta integration step for methods like shootRungeKutta, for
 * a differential equation of the form:
 * 
 *      Q' = A*Q+B
 * 
 * @param   Q0:   Value of the function at step i
 * @param   A0:   Value of A at step i
 * @param   A1:   Value of A at step i+1
 * @param   B0:   Value of B at step i
 * @param   B1:   Value of B at step i+1
 * @param   h:    Integration step
 * @param   step: Direction of integration (1 or -1)
 * @retval  Value of the function at i+1
 */ 
double stepRungeKutta(double Q0, double A0, double A1, double B0, double B1, double h, int step) 
{
    double Amid, Bmid;
    double k1, k2, k3, k4;

    Amid = (A0+A1)/2.0;
    Bmid = (B0+B1)/2.0;

    k1 = (A0 * Q0 + B0) * h * step;
    k2 = (Amid * (Q0 + k1 / 2) + Bmid) * h * step;
    k3 = (Amid * (Q0 + k2 / 2) + Bmid) * h * step;
    k4 = (A1 * (Q0 + k3) + B1) * h * step;

    return Q0 + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
}

/**
 * @brief  Integrate a single ODE
 * @note   Integrate one differential equation of the form:
 * 
 *      Q' = A*Q+B
 * 
 * with a fourth order Runge-Kutta method, up to a given index, either forward or backwards.
 * 
 * @param  &Q: Vector for Q. Will return the integrated values, must contain already the first two as boundary conditions.
 * @param  A: Vector for A (see definition above). Same size as Q.
 * @param  B: Vector for B.
 * @param  h: Integration step (default = 1).
 * @param  stop_i: Index to stop integration. Run the whole range if -1 (default = -1).
 * @param  dir: Integration direction, either forward 'f' or backwards 'b' (default = 'f').
 * @retval None
 */
void shootRungeKutta(vector<double> &Q, vector<double> A, vector<double> B, double h, int stop_i, char dir)
{
    int N = Q.size();
    int step = (dir == 'f') ? 1 : -1;
    int from_i = (step == 1) ? 1 : N - 2;

    // First, check size
    if (A.size() != N || B.size() != N)
    {
        throw invalid_argument("Invalid size for one or more arrays passed to shootQ");
    }

    if (stop_i == -1)
    {
        stop_i = (step == 1) ? N - 1 : 0;
    }

    for (int i = from_i; step * (i - stop_i) <= 0; i += step)
    {       
        Q[i] = stepRungeKutta(Q[i-step], A[i-step], A[i], B[i-step], B[i], h, step);
    }

    return;
}

/**
 * @brief  Integrate two coupled ODEs
 * @note   Integrate a system of coupled differential equations of the form:
 * 
 *      Q' = AA*Q+AB*P
 *      P' = BA*Q+BB*P
 * 
 * with a second order shooting method, up to a given index, either forward or backwards.
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
void shootQP(vector<double> &Q, vector<double> &P, vector<double> AA, vector<double> AB, vector<double> BA, vector<double> BB,
             double h, int stop_i, char dir)
{
    int N = Q.size();
    int step = (dir == 'f') ? 1 : -1;
    int from_i = (step == 1) ? 1 : N - 2;
    double QA, QB, QC, PA, PB, PC;

    double Qp, Pp;
    double AAmid, ABmid, BAmid, BBmid;
    double k1A, k1B, k2A, k2B, k3A, k3B, k4A, k4B;

    // First, check size
    if (P.size() != N || AA.size() != N || AB.size() != N || BA.size() != N || BB.size() != N)
    {
        throw invalid_argument("Invalid size for one or more arrays passed to shootQP");
    }

    if (stop_i == -1)
    {
        stop_i = (step == 1) ? N - 1 : 0;
    }

    for (int i = from_i; step * (i - stop_i) <= 0; i += step)
    {
        AAmid = (AA[i] + AA[i - step]) / 2;
        ABmid = (AB[i] + AB[i - step]) / 2;
        BAmid = (BA[i] + BA[i - step]) / 2;
        BBmid = (BB[i] + BB[i - step]) / 2;
        Pp = P[i - step];
        Qp = Q[i - step];
        k1A = (AA[i - step] * Qp + AB[i - step] * Pp) * h * step;
        k1B = (BA[i - step] * Qp + BB[i - step] * Pp) * h * step;
        k2A = (AAmid * (Qp + k1A / 2.0) + ABmid * (Pp + k1B / 2.0)) * h * step;
        k2B = (BAmid * (Qp + k1A / 2.0) + BBmid * (Pp + k1B / 2.0)) * h * step;
        k3A = (AAmid * (Qp + k2A / 2.0) + ABmid * (Pp + k2B / 2.0)) * h * step;
        k3B = (BAmid * (Qp + k2A / 2.0) + BBmid * (Pp + k2B / 2.0)) * h * step;
        k4A = (AA[i] * (Qp + k3A) + AB[i] * (Pp + k3B)) * h * step;
        k4B = (BA[i] * (Qp + k3A) + BB[i] * (Pp + k3B)) * h * step;

        Q[i] = Qp + 1.0 / 6.0 * (k1A + 2 * k2A + 2 * k3A + k4A);
        P[i] = Pp + 1.0 / 6.0 * (k1B + 2 * k2B + 2 * k3B + k4B);
    }

    return;
}

/**
 * @brief  Integrate a single 2nd order ODE with Numerov's method
 * @note   Integrate one differential equation of the form:
 * 
 *      Q'' = A*Q+B
 * 
 * with Numerov's method, up to a given index, either forward or backwards.
 * 
 * @param  &Q: Vector for Q. Will return the integrated values, must contain already the first two as boundary conditions.
 * @param  A:  Vector for A (see definition above). Same size as Q.
 * @param  B:  Vector for B.
 * @param  h:  Integration step (default = 1).
 * @param  stop_i: Index to stop integration. Run the whole range if -1 (default = -1).
 * @param  dir: Integration direction, either forward 'f' or backwards 'b' (default = 'f').
 * @retval None
 */
void shootNumerov(vector<double> &Q, vector<double> A, vector<double> B, double h, int stop_i, char dir)
{

    int N = Q.size();
    int step = (dir == 'f') ? 1 : -1;
    int from_i = (step == 1) ? 2 : N - 3;
    double h2_12 = h * h / 12;
    double QA0, QA1, QA2, QB;

    // First, check size
    if (A.size() != N || B.size() != N)
    {
        throw invalid_argument("Invalid size for one or more arrays passed to shootNumerov");
    }

    if (stop_i == -1)
    {
        stop_i = (step == 1) ? N - 1 : 0;
    }

    for (int i = from_i; step * (i - stop_i) <= 0; i += step)
    {
        QA0 = (1 - h2_12 * A[i]);
        QA1 = (1 + 5 * h2_12 * A[i - step]);
        QA2 = (1 - h2_12 * A[i - 2 * step]);
        QB = h2_12 * (B[i] + 10 * B[i - step] + B[i - 2 * step]);
        Q[i] = (2 * Q[i - step] * QA1 - Q[i - 2 * step] * QA2 + QB) / QA0;
    }

    return;
}

/**
 * @brief  Integrate a Coulomb potential from a radial background charge density on a logarithmic grid
 * @note   Integrate a Coulomb potential from a radial background charge density on a logarithmic grid.
 * The following ODE is integrated:
 *  
 *      V'' + V' = rho(x)
 * 
 * with the derivatives being in x (r = r0*exp(x)) and the density being already integrated over the angular
 * coordinates (so for example a constant density would have rho ~ r^2).
 * @param  &V:   Vector for V. Will return the integrated potential. 
 * @param  rho:  Vector for the charge density.
 * @param  h:    Integration step (default = 1)
 * @retval None
 */
void shootPotentialLog(vector<double> &V, vector<double> rho, double h)
{
    int N = V.size();
    double h2 = h * h;
    double A0, A1, A2, A3;

    A0 = (2 / h2 + 11 / (6.0 * h));
    A1 = (-5 / h2 - 3 / h);
    A2 = (4 / h2 + 1.5 / h);
    A3 = (-1 / h2 + -1 / (3 * h));

    // Constant charge assumption
    V[0] = rho[0] / 6;
    V[1] = (rho[1] - V[0] * (A1 + A2 * exp(-2 * h) + A3 * exp(-4 * h))) / A0;
    V[2] = (rho[1] - V[1] * A1 - V[0] * (A2 + A3 * exp(-2 * h))) / A0;

    for (int i = 3; i < N; ++i)
    {
        V[i] = (rho[i] - V[i - 1] * A1 - V[i - 2] * A2 - V[i - 3] * A3) / A0;
    }
}

/**
 * @brief  Integrate the radial Dirac equation on a logarithmic grid
 * @note   Perform integration of the radial Dirac equation on a logarithmic grid, forward and backwards, up to the turning point.
 * The coupled equations have the form:
 * 
 *      Q' = k/r*Q+[mc-(E-V)/c]*P
 *      P' = -k/r*P+[mc+(E-V)/c]*Q
 * 
 * With k the quantum number: if j=|l+s|, then k = -(j+1/2)*sign(j-l), and E the expected energy (including the rest mass term).
 * The function will return a struct containing the index of the 'turning point', where the forward and backwards integration meet,
 * as well as the values of Q and P integrated forward (Qi, Pi) and backwards (Qe, Pe) at it.
 *   
 * @param  &Q: Vector for Q. Will return the integrated values, must contain already the first and last two as boundary conditions.
 * @param  &P: Vector for P. Will return the integrated values, must contain already the first and last two as boundary conditions.
 * @param  r:  Radial (logarithmic) grid
 * @param  V:  Potential
 * @param  E:  Energy (binding + mc^2)
 * @param  k:  Quantum number (default = -1)
 * @param  m:  Mass of the particle (default = 1)
 * @param  dx: Integration step (default = 1)
 * @retval turn_i: Turning point index
 */
TurningPoint shootDiracLog(vector<double> &Q, vector<double> &P, vector<double> r, vector<double> V,
                           double E, int k, double m, double dx)
{

    int N = Q.size(), turn_i;
    double B; // Binding energy
    TurningPoint out;
    vector<double> AA(N), AB(N), BA(N), BB(N); // Define the equation

    // Check size
    if (P.size() != N || r.size() != N || V.size() != N)
    {
        throw invalid_argument("Invalid size for one or more arrays passed to shootDiracLog");
    }

    B = E - m * pow(Physical::c, 2);

    // Find the turning point
    for (turn_i = 0; turn_i < V.size(); ++turn_i)
    {
        if (V[turn_i] > B)
            break;
    }
    if (turn_i >= V.size() - 1)
    {
        LOG(ERROR) << "Turning point not included in range: r_max too small\n";
        // Turning point not included in range
        throw TurningPointError(TurningPointError::TPEType::RMAX_SMALL);
    }
    else if (turn_i == 0)
    {
        LOG(ERROR) << "Turning point not included in range: r_min too big\n";
        // Turning point not included in range
        throw TurningPointError(TurningPointError::TPEType::RMIN_BIG);
    }

    // Now define the other arrays
    for (int i = 0; i < N; ++i)
    {
        AA[i] = k;
        AB[i] = -r[i] * (B - V[i]) * Physical::alpha;
        BB[i] = -k;
        BA[i] = r[i] * ((B - V[i]) * Physical::alpha + 2 * m * Physical::c);
    }

    // Integrate forward
    shootQP(Q, P, AA, AB, BA, BB, dx, turn_i + 1);
    out.Qi = Q[turn_i];
    out.Pi = P[turn_i];
    // Integrate backwards
    shootQP(Q, P, AA, AB, BA, BB, dx, turn_i, 'b');
    out.Qe = Q[turn_i];
    out.Pe = P[turn_i];

    out.i = turn_i;

    return out;
}

/**
 * @brief  Integrate d/dE (Q/P) for a Dirac wavefunction based on a Coulomb potential
 * @note   Integrate zeta = d/dE (Q/P) for a Dirac wavefunction based on a Coulomb potential. 
 * Zeta is used to find the optimal energy correction at each step when converging a solution. 
 * 
 * @param  &zeta:   Vector for zeta. Will return the integrated values, must contain already the first and last two as boundary conditions.
 * @param  y:       Vector containing values of Q/P.
 * @param  r:       Radial (logarithmic) grid
 * @param  V:       Potential
 * @param  turn_i:  Index of the turning point. Will mark the stop of the integration.
 * @param  E:       Energy (binding + mc^2)
 * @param  k:       Quantum number (default = -1)
 * @param  m:       Mass of the particle (default = 1)
 * @param  dx:      Integration step (default = 1)
 * @param  dir:     Integration direction, either forward 'f' or backwards 'b' (default = 'f').
 * @retval None
 */
void shootDiracErrorDELog(vector<double> &zeta, vector<double> y, vector<double> r, vector<double> V, int turn_i,
                          double E, int k, double m, double dx, char dir)
{
    int N = zeta.size();
    int step = (dir == 'f') ? 1 : -1;
    int from_i = (step == 1) ? 1 : N - 2;
    double mc = m * Physical::c;
    double g, A0, A1, B0, B1, y02, y12;

    // Check size
    if (y.size() != N || r.size() != N || V.size() != N)
    {
        throw invalid_argument("Invalid size for one or more arrays passed to shootDiracErrorDELog");
    }

    for (int i = from_i; step * (i - turn_i) <= 0; i += step)
    {
        if (abs(y[i]) < Physical::alpha || zeta[i-step] == 0) {
            g = (mc + (E - V[i]) * Physical::alpha);
            A0 = 2*(k-g*r[i-step]*y[i-step]);
            A1 = 2*(k-g*r[i]*y[i]);
            B0 = -r[i-step]*(1+pow(y[i-step], 2))*Physical::alpha;
            B1 = -r[i]*(1+pow(y[i], 2))*Physical::alpha;

            zeta[i] = stepRungeKutta(zeta[i-step], A0, A1, B0, B1, dx, step);
        }
        else {
            g = (mc - (E - V[i]) * Physical::alpha);
            y02 = pow(y[i-step], 2);
            y12 = pow(y[i], 2);
            A0 = -2*(k+g*r[i-step]/y[i-step]);
            A1 = -2*(k+g*r[i]/y[i]);
            B0 = r[i-step]*(1+1/y02)*Physical::alpha;
            B1 = r[i]*(1+1/y12)*Physical::alpha;

            zeta[i] = -y12*stepRungeKutta(-zeta[i-step]/y02, A0, A1, B0, B1, dx, step);
        }
    }
}