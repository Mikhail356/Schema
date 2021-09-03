import math
import os

eps = 1e-3
tau = 0.0001
M   = 100

mu_set  = [0.1, 0.01]
C_set   = [1., 10.]

for No in [2]:
    for mu in mu_set:
        for C in C_set:
            print (f'({mu} {C}) ', end = '')
    print ()
    print (r'\begin{table}[H]')
    print (r'\centering')
    print (r'\begin{tabular}{|c|c|c|c|c|c|}')
    print (r'\hline')
    print (r'$C$ & $C = 10^0$ & $C = 10^1$ & $\gamma$ \\')
    print (r'\hline')
    for k in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        print (k, end = '')
        for mu in mu_set:
            for C in C_set:
                s = f'./a.out {No} {mu} {C} {M} {tau} {k} {eps}'
                y = [float (t) for t in os.popen (s).read ().split ()]
                print (' & \\texttt{' + str (y[0]) + '}', end = '')
        print (' \\\\')
        print ('\\hline')
    print ()
    print ()

