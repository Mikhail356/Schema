import math
import os

eps = 1e-3
tau = 0.0001
M   = 1000

mu_set  = [0.1, 0.01]
C_set   = [1., 10.]

for No in [2]:
    for C in C_set:
        print (r'\begin{table}[H]')
        print (r'\centering')
        print (r'\begin{tabular}{|c|c|c|c|c|c|}')
        print (r'\hline')
        print (r'$\mu$ & $n = N_0/5$ & $n = 2N_0/5$ & $n = 3N_0/5$ & $n = 4N_0/5$ & $n = N_0$\\')
        print (r'\hline')

        for mu in mu_set:
            k = 0

            for t_it in [tau, tau / 2]:
                for m_it in [M, M * 2]:
                    s = f'./a.out {No} {mu} {C} {m_it} {t_it} {eps}'
                    y = [float (t) for t in os.popen (s).read ().split ()]

                    n = len (y) - 1
                    if k == 0:
                        print ('\\multirow{2}{*}{$10^{' + str (int (math.log10 (mu))) + '}$}', end = '')
                    for ind in [n/5, 2*n/5, 3*n/5, 4*n/5, n]:
                        print (' & \\texttt{' + '{:e}'.format (y[int (ind)]) + '}', end = '')
                    print (' \\\\')
                    k += 1

            print ('\\hline')

        print (r'\end{tabular}')
        print ('\\caption{Разность масс при $p(\\rho) = 10^{' + str (int (math.log10 (C))) + '}\\rho$}')
        print (r'\end{table}')
        print ()
        print ()

