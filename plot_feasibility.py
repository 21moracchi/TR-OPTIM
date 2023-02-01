
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

dump_dir = "/Users/leonard/opf/"

pg3_min = 0.1
pg3_max = 2.7
pg2_min = 0.1
pg2_max = 3.0

# Load data
conv = np.loadtxt(dump_dir + "feasibility.txt")
conv = conv[::-1, :]

# Resolution
resolution = np.size(conv, 1)
x = np.linspace(pg3_min, pg3_max, resolution)
y = np.linspace(pg2_min, pg2_max, resolution)

X, Y = np.meshgrid(x, y)

fig, ax = plt.subplots(figsize=(6, 5))

levels = [-0.1, 0.1]

c1 = np.copy(conv)
c1[c1 > 0] = 1.0
CS = ax.contourf(X, Y, c1, levels, corner_mask=False, colors=["b"], alpha=.4)

c2 = np.copy(conv)
c2[c2 != 2.0] = 0.0
CS2 = ax.contour(X, Y, c2, 0)
fmt = {}
fmt[CS2.levels[0]] = "Voltage"
ax.clabel(CS2, fontsize=9, inline=True, fmt=fmt)

c3 = np.copy(conv)
c3[c3 != 4.0] = 0.0
CS3 = ax.contour(X, Y, c3, 0)
fmt = {}
fmt[CS3.levels[0]] = "Power gen."
ax.clabel(CS3, fontsize=9, inline=True, fmt=fmt)

c4 = np.copy(conv)
c4[c4 != 8.0] = 0.0
CS4 = ax.contour(X, Y, c4, 0, color="yellow")
fmt = {}
fmt[CS4.levels[0]] = "Line flow"
ax.clabel(CS4, fontsize=9, inline=True, fmt=fmt)

ax.set_xlim(pg3_min, pg3_max)
ax.set_ylim(pg2_min, pg2_max)

ax.set_xlabel("Active power generation P3")
ax.set_ylabel("Active power generation P2")
legend_elements = [Patch(facecolor='blue', edgecolor='blue',  alpha=.4, label='Feasible space')]
ax.legend(handles=legend_elements)
ax.set_title("Feasibility space - case9")
plt.tight_layout()
plt.savefig("feasible_space.png")
