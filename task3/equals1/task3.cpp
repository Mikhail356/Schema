#include <stdio.h>
#include <math.h>

#include <vector>

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
  double eps;
  double tau;
  double h;
} params_t;

using arr_t = std::vector<double>;

// Подзадача 1. Тут везде подразумевается t = 0
double rho_1 (double x, int N);
double g_1 (double x, int N);
double u_1 (double x, int N);

// Подзадача 2. Тут везде подразумевается t = 0
double rho_2 (double x, int N);
double g_2 (double x, int N);
double u_2 (double x, int N);

double tilde_p (const data_t & data, double g);

double residual (const arr_t & a, const arr_t & b, int M);
double norm (const arr_t & v, int M);

void tridiagonal_solve (arr_t & a, arr_t & b, arr_t & c, arr_t & d, int n);

int scheme (const data_t & data, const params_t & params, arr_t & V_curr, arr_t & G_curr, arr_t & V_next, arr_t & G_next);

int
main (int argc, char * argv [])
{

  if (argc != 8)
    {
      fprintf (stderr, "Use: %s task_No Mu C M tau N eps\n", argv[0]);
      return -1;
    }

  int    task = atoi (argv[1]);
  double mu   = atof (argv[2]);
  double C    = atof (argv[3]);
  int    M    = atoi (argv[4]);
  double tau  = atof (argv[5]);
  int    N    = atoi (argv[6]);
  double eps  = atof (argv[7]);

  if (task != 1 && task != 2)
    {
      fprintf (stderr, "Unknown task (task must be equal 1 or 2)\n");
      return -1;
    }

  if (mu <= 0.)
    {
      fprintf (stderr, "Input param mu incorrect (mu must be great 0)\n");
      return -1;
    }

  if (M <= 0)
    {
      fprintf (stderr, "Input param M incorrect (M must be great 0)\n");
      return -1;
    }

  if (tau <= 0.)
    {
      fprintf (stderr, "Input param tau incorrect (tau must be great 0)\n");
      return -1;
    }

  if (eps <= 0.)
    {
      fprintf (stderr, "Input param eps incorrect (eps must be great 0)\n");
      return -1;
    }

  if (N <= 0 || N > M / 10)
    {
      fprintf (stderr, "Input param N incorrect (N must be in (0, M/10] )\n");
      return -1;
    }

  data_t data;
  // data.T = 0.;
  data.X = 1.;
  data.C = C;
  data.gamma = 0.;
  data.mu = mu;

  params_t params;
  params.M = M;
  params.eps = eps;
  params.tau = tau;
  params.h = data.X / M;

  arr_t V_curr (M + 1);
  arr_t G_curr (M + 1);
  arr_t V_next (M + 1);
  arr_t G_next (M + 1);

  for (int m = 0; m <= M; m++)
    {
      double x = m * params.h;
      V_curr[m] = (task == 1) ? u_1 (x, N) : u_2 (x, N);
      G_curr[m] = (task == 1) ? g_1 (x, N) : g_2 (x, N);
    }

  int st = scheme (data, params, V_curr, G_curr, V_next, G_next);

  printf ("%f\n", st * tau);

}

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

double
rho_1 (double x, int N)
{
  return 2. + sin (N * M_PI * x);
}

double
g_1 (double x, int N)
{
  return log (rho_1 (x, N));
}

double
u_1 (double x, int N)
{
  return 0.;
}

/* ================================================================================================================== */

double
rho_2 (double x, int N)
{
  return 1.;
}

double
g_2 (double x, int N)
{
  return 0.;
}

double
u_2 (double x, int N)
{
  return sin (N * M_PI * x);
}

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

double
tilde_p (const data_t & data, double g)
{
  double C = data.C;
  return (C < 0) ? 1.4 * exp (0.4 * g) : C;
}

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

double
residual (const arr_t & a, const arr_t & b, int M)
{
  double residual = 0.;
  for (int m = 0; m <= M; m++)
    {
      residual += fabs (a[m] - b[m]);
    }
  return residual;
}

double
norm (const arr_t & v, int M)
{
  double norm = 0.;
  for (int m = 0; m <= M; m++)
    {
      norm = std::max (norm, fabs (v[m]));
    }
  return norm;
}

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

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

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

int
scheme (const data_t & data, const params_t & params, arr_t & V_curr, arr_t & G_curr, arr_t & V_next, arr_t & G_next)
{
  int N = params.N;
  int M = params.M;
  double tau = params.tau;
  double h = params.h;

  arr_t a (M + 1);
  arr_t b (M + 1);
  arr_t c (M + 1);

  for (int n = 0; ; n++)
    {
      double t = n * tau;

      /* === INIT SYSTEM FOR G === */

      // first equation
      b[0] = 1.;
      c[0] = 0.;
      G_next[0] = G_curr[0] - (tau * (V_curr[1] - V_curr[0]) / h) /*+ tau * f0 (t, 0, data)*/;

      for (int m = 1; m < M; m++)
        {
          double x = m * h;
          a[m] = - (V_curr[m] + fabs (V_curr[m])) / (2. * h);
          b[m] = 1. / tau + fabs (V_curr[m]) / h;
          c[m] = (V_curr[m] - fabs (V_curr[m])) / (2. * h);
          G_next[m] = G_curr[m] / tau - (V_curr[m + 1] - V_curr[m - 1]) / (2. * h) /*+ f0 (t, x, data)*/;
        }

      // last equation
      a[M] = 0.;
      b[M] = 1.;
      G_next[M] = /*tau * f0 (t, M * h, data) +*/ G_curr[M] - tau * (V_next[M] - V_next[M - 1]) / h;


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
                      /*+ f (t, x, data)*/;
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

//      double res = std::max (residual (V_curr, V_next, M), residual (G_curr, G_next, M));
      double res = norm (V_next, M);

      std::swap (V_curr, V_next);
      std::swap (G_curr, G_next);

      if (n >= 999 && res < params.eps)
        {
          return n + 1;
        }

    }

}
