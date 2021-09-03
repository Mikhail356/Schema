import matplotlib.pyplot as plt
import numpy as np

plt.rc ('text', usetex = True)
plt.rc ('text.latex', preamble = r'\usepackage[charter, expert, cal = cmcal]{mathdesign}')

fig = plt.figure ()
ax = fig.add_subplot (111)
ax.set_aspect ('auto')

with open ('task2_log.txt', 'r') as f:
    t, h = [float (x) for x in next (f).split ()]
    Z = [[float (x) for x in line.split ()] for line in f]

Z_V = [Z[i] for i in range (0, len (Z), 2)]
Z_G = [Z[i] for i in range (1, len (Z), 2)]

x_val = [h * m for m in range (len (Z[0]))]
y_val = [t * n for n in range (len (Z_V))]

X, Y = np.meshgrid (x_val, y_val)

cp = plt.contourf (Y, X, Z_V, 1000, cmap = 'nipy_spectral')
plt.colorbar (cp)
plt.savefig ('plot_V.png', dpi = 300)

