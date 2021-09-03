import math
import os

eps = 1e-3
tau = 0.0001
M   = 1000

mu_set  = [0.1, 0.01, 0.001]
C_set   = [1., 10., 100.]

with open ('task2_4.log', 'w') as f:
    for No in [1,2]:
        f.write (r'\begin{table}[H]' + '\n')
        f.write (r'\centering' + '\n')
        f.write (r'\begin{tabular}{|c|c|c|c|c|c|}' + '\n')
        f.write (r'\hline' + '\n')
        f.write (r'123123 & $10^0$ & $10^1$ & $10^2$ & $\gamma$ \\' + '\n')
        f.write (r'\hline' + '\n')

        for mu in mu_set:

            k = 0
            for t_it in [tau, tau / 2]:
                for m_it in [M, M * 2]:

                    f.write ('$10^{' + str (int (math.log10 (mu))) + '}$')

                    for C in C_set:
                        s = f'./a.out {No} {mu} {C} {m_it} {t_it} {eps}'
                        print (s)
                        y = float (os.popen (s).read ())

                        f.write (' & \\texttt{' + str (y) + '}')
                        f.flush ()
                    f.write (' \\\\' + '\n')
                    k += 1

            f.write ('\\hline' + '\n')

        f.write (r'\end{tabular}' + '\n')
        f.write (r'%\caption{Времена сходимости}' + '\n')
        f.write (r'\end{table}' + '\n\n\n')
        f.flush ()


