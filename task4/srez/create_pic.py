import matplotlib.pyplot as plt
import numpy as np
import sys

plt.rc ('text', usetex = True)
plt.rc ('text.latex', preamble = r'\usepackage[charter, expert, cal = cmcal]{mathdesign}')

fig = plt.figure ()
ax = fig.add_subplot (111)

with open ('task2_log.txt', 'r') as f:
    t, h = [float (x) for x in next (f).split ()]
    Z = [[float (x) for x in line.split ()] for line in f]

Z_V = [Z[i] for i in range (0, len (Z), 2)]
Z_G = [Z[i] for i in range (1, len (Z), 2)]

x_val = [h * m for m in range (len (Z[0]))]
y_val = [t * n for n in range (len (Z_V))]

X, Y = np.meshgrid (x_val, y_val)

arr = [0, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 5000, 5100, 5200, 5300, 5400, 5500, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 30000, 50000, len (Z_V) - 1] + [t for t in range (5050, 5100, 10)]
arr.sort ()

for i in arr:
    if len (Z_V) > i:
        line1 = ax.plot (x_val, Z_V[i], label = '$V$', linewidth = 1)
        line2 = ax.plot (x_val, Z_G[i], label = '$G$', linewidth = 1)
        #ax.set_aspect (get_aspect (x_val, Z_V[i] + Z_G[i]))
        ax.set_aspect ('auto')
        plt.legend ()
        name = f'plot_{i}.pdf'
        plt.savefig (name)
        print ('\\inserpicture{res/' + name + '}{Срез для $t = ' + str (i * t) + '$}')
        plt.cla ()
print ()
print ()
