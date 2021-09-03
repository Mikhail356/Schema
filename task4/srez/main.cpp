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
  double u_left;
  double g_left;
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

double tilde_p (const data_t & data, double g);

double residual (const arr_t & a, const arr_t & b, int M);
double norm (const arr_t & v, int M);

void tridiagonal_solve (arr_t & a, arr_t & b, arr_t & c, arr_t & d, int n);

void scheme (const data_t & data, const params_t & params, arr_t & V_curr, arr_t & G_curr, arr_t & V_next, arr_t & G_next);

int
main (int argc, char * argv [])
{

  if (argc != 8)
    {
      fprintf (stderr, "Use: %s Mu C M tau u_left g_left eps\n", argv[0]);
      return -1;
    }

  double mu     = atof (argv[1]);
  double C      = atof (argv[2]);
  int    M      = atoi (argv[3]);
  double tau    = atof (argv[4]);
  double u_left = atof (argv[5]);
  double g_left = atof (argv[6]);
  double eps    = atof (argv[7]);

//  printf ("%s Mu=%s C=%s M=%s tau=%s u_left=%s g_left=%s eps=%s\n",
//          argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);

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

  if (u_left <= 0)
    {
      fprintf (stderr, "Input param u_left incorrect (u_left must be great 0)\n");
      return -1;
    }

  if (g_left < 0)
    {
      fprintf (stderr, "Input param g_left incorrect (g_left must be great 0)\n");
      return -1;
    }

  if (eps <= 0.)
    {
      fprintf (stderr, "Input param eps incorrect (eps must be great 0)\n");
      return -1;
    }

  data_t data;
  // data.T = 0.;
  data.X = 10.;
  data.C = C;
  data.gamma = 0.;
  data.mu = mu;

  data.u_left = u_left;
  data.g_left = g_left;

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
      V_curr[m] = G_curr[m] = 0.;
//      if (m <= M/10)
//        {
//          V_curr[m] = u_left;
//          G_curr[m] = g_left;
//        }
    }

//  V_curr[0] = u_left;
//  G_curr[0] = g_left;

  scheme (data, params, V_curr, G_curr, V_next, G_next);

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

void
scheme (const data_t & data, const params_t & params, arr_t & V_curr, arr_t & G_curr, arr_t & V_next, arr_t & G_next)
{
  int N = params.N;
  int M = params.M;
  double tau = params.tau;
  double h = params.h;

  arr_t a (M + 1);
  arr_t b (M + 1);
  arr_t c (M + 1);


  FILE * f = fopen ("task2_log.txt", "w");

  if (f == nullptr)
    {
      return;
    }

  fprintf (f, "%e %e\n", tau, h);

  auto dump = [&] (const arr_t & a, double t)
    {
      for (int m = 0; m <= M; m++)
        {
          fprintf (f, "%+e ", a[m]);
        }
      fprintf (f, "\n");
    };

  dump (V_curr, 0.);
  dump (G_curr, 0.);

  for (int n = 0; n < 100000; n++)
    {
      double t = n * tau;

      /* === INIT SYSTEM FOR G === */

      // first equation
      a[0] = 0.;
      b[0] = 1.;
      c[0] = 1.;
      G_next[0] = data.g_left;

      for (int m = 1; m < M; m++)
        {
          double x = m * h;
          a[m] = - (V_curr[m] + fabs (V_curr[m])) / (2. * h);
          b[m] = 1. / tau + fabs (V_curr[m]) / h;
          c[m] = (V_curr[m] - fabs (V_curr[m])) / (2. * h);
          G_next[m] = G_curr[m] / tau - (V_curr[m + 1] - V_curr[m - 1]) / (2. * h) /*+ f0 (t, x, data)*/;
        }

      // last equation
//      a[M] = 0.;
//      b[M] = 1.;
//      G_next[M] =  G_curr[M] - tau * (V_curr[M] - V_curr[M - 1]) / h;
////      G_next[M] = G_curr[M] - tau * V_curr[M] * (G_curr[M] - G_curr[M - 1]) / h;
//      G_next[M] = (G_curr[M]/tau + (G_next[M-1]*V_curr[M])/h)/(1./tau + V_curr[M]/h);

      a[M] = - V_curr[M] / h;
      b[M] = 1. / tau + V_curr[M] / h;
      c[M] = 0;
      G_next[M] = G_curr[M] / tau;

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
      a[0] = 0.;
      b[0] = 1.;
      c[0] = 0.;
      V_next[0] = data.u_left;
//      b[0] = 1.;
//      c[0] = -1.;
//      V_next[0] = 0;

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
      a[M] = -1.;
      b[M] = 1.;
      c[M] = 0.;
      V_next[M] = /*V_next[M-1]*/0.;

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
//      double res = norm (V_next, M);

//      if (n > 90000)
        {
          dump (V_next, (n + 1) * tau);
          dump (G_next, (n + 1) * tau);
        }

      std::swap (V_curr, V_next);
      std::swap (G_curr, G_next);

    }

  fclose (f);

}
