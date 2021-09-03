import math
import os

eps = 1e-3
tau = 0.0001
M   = 1000

mu_set  = [0.1, 0.01]
C_set   = [1., 10.]

for No in [1]:
    for C in C_set:
        print (r'\begin{table}[H]')
        print (r'\centering')
        print (r'\begin{tabular}{|c|c|c|c|c|c|}')
        print (r'\hline')
        print (r'$\mu$ & $T = N_0\tau$ & $n = N_0/4$ & $n = N_0/2$ & $n = 3N_0/4$ & $n = N_0$ \\')
        print (r'\hline')

        for mu in mu_set:
            k = 0

            for t_it in [tau, tau / 2]:
                for m_it in [M, M * 2]:
                    s = f'./a.out {No} {mu} {C} {m_it} {t_it} {eps}'
                    y = [float (t) for t in os.popen (s).read ().split ()]
                    #x = [t_it * n for n in range (len (y))]
                    time = y[-1]

                    n = len (y) - 2
                    if k == 1:
                        print ('\\multirow{2}{*}{$10^{' + str (int (math.log10 (mu))) + '}$}', end = '')
                    print (' & \\texttt{' + str (time) + '}', end = '')
                    for ind in [n/4, n/2, 3*n/4, n]:
                        print (' & \\texttt{' + '{:e}'.format (y[int (ind)]) + '}', end = '')
                    print (' \\\\')
                    k += 1

            print ('\\hline')
        
        print (r'\end{tabular}')
        print ('\\caption{Точности решения для при $p(\\rho) = 10^{' + str (int (math.log10 (C))) + '}\\rho$}')
        print (r'\end{table}')
        print ()
        print ()

