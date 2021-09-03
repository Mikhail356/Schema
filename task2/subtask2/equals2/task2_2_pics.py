import matplotlib.pyplot as plt
import math
import os

plt.rc ('text', usetex = True)
plt.rc ('text.latex', preamble = r'\usepackage[charter, expert, cal = cmcal]{mathdesign}')

eps = 1e-3
tau = 0.0001
M   = 1000

mu_set  = [0.1, 0.01]
C_set   = [1., 10.]

def to_num (s):
    num = float (s)
    if num > 1.:
        return 1.
    if num < -1.:
        return -1.
    return num

def get_aspect (x, y):
    x_min = x_max = x[0]
    y_min = y_max = y[0]

    for t in x:
        if t > x_max:
            x_max = t
        if t < x_min:
            x_min = t

    for t in y:
        if t > y_max:
            y_max = t
        if t < y_min:
            y_min = t

    return (1 / 3) * (x_max - x_min) / (y_max - y_min)

for No in [2]:
    k = 0
    for C in C_set:
        for mu in mu_set:

            k = k + 1

            fig = plt.figure ()
            ax = fig.add_subplot (111)

            x_data, y_data, lab_data = [], [], []

            for t_it in [1, 2]:
                for m_it in [1, 2]:
                    s = f'./a.out {No} {mu} {C} {m_it * M} {tau / t_it} {eps}'
                    Y = [to_num (t) for t in os.popen (s).read ().split ()]
                    X = [(tau / t_it) * n for n in range (len (Y))]

                    #n = int (len (Y) / 500)
                    #if n == 0:
                    #    n = 1
                    #x = [1000 * t for t in X[::n]] + [1000 * X[-1]]
                    #y = Y[::n] + [Y[-1]]
                    x = [1000 * t for t in X[::2]]
                    y = Y[::2]

                    p1 = '/2' if t_it == 2 else ''
                    p2 = '/2' if m_it == 2 else ''
                    lab = '$\\Omega_{\\tau' + p1 + ', \\, h' + p2 + '}$'
                    #ax.plot (x, y, label = lab)
                    x_data += [x]
                    y_data += [y]
                    lab_data += [lab]
            #ax.set_aspect (get_aspect (x_data[0] + x_data[1] + x_data[2] + x_data[3], y_data[0] + y_data[1] + y_data[2] + y_data[3]))
            ax.set_aspect ('auto')
            ax.plot (x_data[0], y_data[0], label = lab_data[0], linewidth = 1)
            ax.plot (x_data[1], y_data[1], label = lab_data[1], linewidth = 1)
            ax.plot (x_data[2], y_data[2], label = lab_data[2], linewidth = 1)
            ax.plot (x_data[3], y_data[3], label = lab_data[3], linewidth = 1)
            plt.legend ()
            plt.savefig (f'plot_2_{No}_{k}.pdf')
            mu_pow = int (math.log10 (mu))
            C_pow = int (math.log10 (C))
            print ('\\inserpicture{res/' + f'plot_2_{No}_{k}.pdf' + '}{Графики функций $\\Delta_m$ для $C = 10^{' + str (C_pow) + '}$ и $\mu = 10^{' + str (mu_pow) + '}$}')
            plt.cla ()





