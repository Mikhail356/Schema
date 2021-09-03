#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

#define EPS  1e-16

typedef struct
{
  double T;
  double X;
  double C;     // if p (rho) = C*rho
  double gamma; // if p (rho) = rho^gamma
  double mu;
} data_t;

typedef struct
{
  int N;
  int M;
  double tau;
  double h;
} params_t;

using arr_t = std::vector<double>;
using namespace std;

double rho (double t, double x);
double rho_t (double t, double x);
double rho_x (double t, double x);

double g (double t, double x);
double g_t (double /*t*/, double /*x*/);
double g_x (double t, double x);

double u (double t, double x);
double u_t (double t, double x);
double u_x (double t, double x);
double u_xx (double t, double x);

double f0 (double t, double x, const data_t & d);
double f (double t, double x, const data_t & d);

double tilde_p (const data_t & data, double g);

void tridiagonal_solve (arr_t & a, arr_t & b, arr_t & c, arr_t & d, int n);

void scheme (const data_t & data, const params_t & params,
             arr_t & V_curr, arr_t & G_curr, arr_t & V_next, arr_t & G_next);


int
main(void)
{
  double mu_set[3] {0.1, 0.01, 0.001};
  double C_set [3] {1., 10., 100.};
  double N_set [4] {10, 100, 1000, 10000};
  double M_set [4] {100, 1000, 10000, 100000};
  double mu, C;
  int N, M;

  FILE * file_v = fopen ("V_res.tex", "w");
  FILE * file_g = fopen ("G_res.tex", "w");

  for (int i = 0; i < 3; i++)
    {
      C = C_set[i];
      for (int j = 0; j < 3; j++)
        {
          mu = mu_set[j];
          data_t data;
          data.T = 1.;
          data.X = 10.;
          data.C = C;
          data.gamma = 0.;
          data.mu = mu;

          fprintf (file_v, "\\begin{table}\n\\centering\n\\begin{tabular}"
                           "{|c|cccc|}\n\\hline\n");
          fprintf (file_g, "\\begin{table}\n\\centering\n\\begin{tabular}"
                           "{|c|cccc|}\n\\hline\n");
          fprintf (file_v, "{\\diagbox{\\boldmath$\\tau$}{\\boldmath$ h$}} & "
                           "\\boldmath $10^{-1}$ & \\boldmath $10^{-2}$ & "
                           "\\boldmath $10^{-3}$ & \\boldmath $10^{-4}$ "
                           "\\\\\n\\hline\n");
          fprintf (file_g, "{\\diagbox{\\boldmath$\\tau$}{\\boldmath$ h$}} "
                           "& \\boldmath $10^{-1}$ & \\boldmath $10^{-2}$ & "
                           "\\boldmath $10^{-3}$ & \\boldmath $10^{-4}$ \\\\\n"
                           "\\hline\n");

          printf ("mu = %e ; C = %e\n", mu, C);

          for (int k = 0; k < 4; k++)
            {
              N = N_set[k];
              printf ("N = %d | ", N);
              fprintf (file_v, "\\boldmath $10^{-%d}$", k+1);
              fprintf (file_g, "\\boldmath $10^{-%d}$", k+1);

              for (int l = 0; l < 4; l++)
                {
                  M = M_set[l];
                  params_t params;
                  params.N = N;
                  params.M = M;
                  params.tau = data.T / N;
                  params.h = data.X / M;

                  vector<double> V_curr (M + 1);
                  vector<double> G_curr (M + 1);
                  vector<double> V_next (M + 1);
                  vector<double> G_next (M + 1);

                  for (int m = 0; m <= M; m++)
                    {
                      double t = 0.;
                      double x = m * params.h;
                      V_curr[m] = u (t, x);
                      G_curr[m] = g (t, x);
                    }

                  scheme (data, params, V_curr, G_curr, V_next, G_next);

                  double t = N * params.tau;
                  double residual_v = 0.;
                  double residual_g = 0.;

                  for (int m = 0; m <= M; m++)
                    {
                      double x = m * params.h;

                      double val_v = fabs (V_curr[m] - u (t, x));
                      double val_g = fabs (G_curr[m] - g (t, x));

                      residual_v = std::max (residual_v, val_v);
                      residual_g = std::max (residual_g, val_g);
                    }

                  if (residual_v < 0.)
                    {
                      residual_v = NAN;
                    }

                  if (residual_g < 0.)
                    {
                      residual_g = NAN;
                    }

                  printf ("%e (%e) | ", residual_v, residual_g);
                  fprintf (file_v, " & \\texttt{%e}", residual_v);
                  fprintf (file_g, " & \\texttt{%e}", residual_g);

                  fflush (file_v);
                  fflush (file_g);
                  fflush (stdout);
                }

              printf("\n\n");
              fprintf (file_v, " \\\\\n");
              fprintf (file_g, " \\\\\n");

              fflush (file_v);
              fflush (file_g);
              fflush (stdout);
            }

          fprintf (file_v, "\\hline\n\\end{tabular}\n\\caption{Ошибка решения для $V$ при $\\mu = 10^{%d}$ и $C = 10^{%d}$}\n\\end{table}\n\n\n", (int) log10 (mu), (int) log10 (C));
          fprintf (file_g, "\\hline\n\\end{tabular}\n\\caption{Ошибка решения для $G$ при $\\mu = 10^{%d}$ и $C = 10^{%d}$}\n\\end{table}\n\n\n", (int) log10 (mu), (int) log10 (C));

          fflush (file_v);
          fflush (file_g);
          fflush (stdout);
        }
    }
  fclose (file_v);
  fclose (file_g);
}

double
rho (double t, double x)
{
  return exp (t) * (cos (M_PI * x / 10.) + 1.5);
}

double
rho_t (double t, double x)
{
  return rho (t, x);
}

double
rho_x (double t, double x)
{
  return - M_PI * exp (t) * sin (M_PI * x / 10.) / 10.;
}

double
g (double t, double x)
{
  return log (rho (t, x));
}

double
g_t (double /*t*/, double /*x*/)
{
  return 1.;
}

double
g_x (double t, double x)
{
  return rho_x (t, x) / rho (t, x);
}

double
u (double t, double x)
{
  return cos (2. * M_PI * t) * sin (M_PI * x * x / 100.);
}

double
u_t (double t, double x)
{
  return -2. * M_PI * sin (2. * M_PI * t) * sin (M_PI * x * x / 100.);
}

double
u_x (double t, double x)
{
  return M_PI * x * cos (2. * M_PI * t) * cos (M_PI * x * x / 100.) / 50.;
}

double
u_xx (double t, double x)
{
  double ang = M_PI * x * x / 100.;
  double tmp = cos (ang) - (ang * sin (ang) * 2);
  return M_PI * cos (2. * M_PI * t) * tmp / 50.;
}

double
f0 (double t, double x, const data_t & /*data*/)
{
  return g_t (t, x) + u (t, x) * g_x (t, x) + u_x (t, x);
}

double
f (double t, double x, const data_t & data)
{
  double p_tilde = tilde_p (data, g (t, x));
  double rhs = data.mu * exp (-g (t, x)) * u_xx (t, x);
  return u_t (t, x) + u (t, x) * u_x (t, x) + p_tilde * g_x (t, x) - rhs;
}

double
tilde_p (const data_t & data, double /*g*/)
{
//  double C = data.C;
//  return (C < 0) ? 1.4 * exp (0.4 * g) : C;
  return data.C;
}

void
tridiagonal_solve (arr_t & a, arr_t & b, arr_t & c, arr_t & d, int n)
{
  /*
  source: https://gist.github.com/nlw0/8edc1241bd05d5a9e5483bee763696a8
  // n is the number of unknowns
  |b0 c0 0 ||x0| |d0|
  |a1 b1 c1||x1|=|d1|
  |0  a2 b2||x2| |d2|
  1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->
      x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0
  2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
      from 1st it.: -| a1x0 + a1g0x1        = a1r0
                  -----------------------------
                        (b1 - a1g0)x1 + c1x2 = d1 - a1r0
      x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)
  3rd iteration:      | a2x1 + b2x2   = d2
      from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                     -----------------------
                     (b2 - a2g1)x2 = d2 - a2r2
      x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
  Finally we have a triangular matrix:
  |1  g0 0 ||x0| |r0|
  |0  1  g1||x1|=|r1|
  |0  0  1 ||x2| |r2|
  Condition: ||bi|| > ||ai|| + ||ci||
  in this version the c matrix reused instead of g
  and             the d matrix reused instead of r and x matrices to report results
  Written by Keivan Moradi, 2014
  */
  n--; // since we start from x0 (not x1)
  c[0] /= b[0];
  d[0] /= b[0];
  for (int i = 1; i < n; i++)
    {
      c[i] /= b[i] - a[i] * c[i - 1];
      d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }
  d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);
  for (int i = n; i-- > 0;)
    {
      d[i] -= c[i] * d[i + 1];
    }
}

void
scheme (const data_t & data, const params_t & params, arr_t & V_curr,
        arr_t & G_curr, arr_t & V_next, arr_t & G_next)
{
  int N = params.N;
  int M = params.M;
  double tau = params.tau;
  double h = params.h;

  arr_t a (M + 1);
  arr_t b (M + 1);
  arr_t c (M + 1);

  for (int n = 0; n < N; n++)
    {
      double t = n * tau;

      /* === INIT SYSTEM FOR G === */

      // first equation
      b[0] = 1.;
      c[0] = 0.;
      G_next[0] = G_curr[0] - (tau * (V_curr[1] - V_curr[0]) / h) + tau * f0 (t, 0, data);

      for (int m = 1; m < M; m++)
        {
          double x = m * h;
          a[m] = - (V_curr[m] + fabs (V_curr[m])) / (2. * h);
          b[m] = 1. / tau + fabs (V_curr[m]) / h;
          c[m] = (V_curr[m] - fabs (V_curr[m])) / (2. * h);
          G_next[m] = G_curr[m] / tau - (V_curr[m + 1] - V_curr[m - 1]) / (2. * h) + f0 (t, x, data);
        }

      // last equation
      a[M] = 0.;
      b[M] = 1.;
      G_next[M] = tau * f0 (t, M * h, data) + G_curr[M] - tau * (V_next[M] - V_next[M - 1]) / h;


      // для отладки (точное значение)
//      for (int m = 0; m <= M; m++)
//        {
//          double t = (n + 1) * tau;
//          double x = m * h;
//          G_next[m] = g (t, x);
//        }

      /* === SOLVE SYSTEM FOR G === */
      tridiagonal_solve (a, b, c, G_next, M + 1);

      double power = 0.;

      for (int m = 0; m <= M; m++)
        {
          power = std::min (power, G_next[m]);
        }

      double mu_tilde = data.mu * exp (-power);


      /* === INIT SYSTEM FOR V === */

      // first equation
      b[0] = 1.;
      c[0] = 0.;
      V_next[0] = 0.;

      for (int m = 1; m < M; m++)
        {
          double x = m * h;
          double p_tilde = tilde_p (data, G_next[m]);
          a[m] = - (V_curr[m] + fabs (V_curr[m])) / (2. * h) - mu_tilde / (h * h);
          b[m] = 1. / tau + fabs (V_curr[m]) / h + 2. * mu_tilde / (h * h);
          c[m] = (V_curr[m] - fabs (V_curr[m])) / (2. * h) - mu_tilde / (h * h);
          V_next[m] =   V_curr[m] / tau
                      - p_tilde * (G_next[m + 1] - G_next[m - 1]) / (2. * h)
                      - (mu_tilde - data.mu * exp (-G_next[m])) * (V_curr[m - 1] - 2. * V_curr[m] + V_curr[m + 1]) / (h * h)
                      + f (t, x, data);
        }

      // last equation
      a[M] = 0.;
      b[M] = 1.;
      V_next[M] = 0.;

      // для отладки (точное значение)
//      for (int m = 0; m <= M; m++)
//        {
//          double t = (n + 1) * tau;
//          double x = m * h;
//          V_next[m] = u (t, x);
//        }

      /* === SOLVE SYSTEM FOR V === */
      tridiagonal_solve (a, b, c, V_next, M + 1);

      std::swap (V_curr, V_next);
      std::swap (G_curr, G_next);

    }

}
