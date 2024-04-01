import raytrace as rt
import loadmap as me
import numpy as np
import matplotlib.pyplot as plt
from math import log10

xs = np.linspace(.01, 600, 1000)
ys = []

env = me.Environment()
env.load("single_tunnel.xml")

pos = [7, 2]

pwr = 3 # dBm

pwr = pow(10, pwr / 10) / 1000

zs = np.linspace(-2, 2, 100)

stns = []
for z in zs:
    row = []
    for x in xs:
        y_raw = rt.raytrace_loss((z, 0, 2), [z, x, 2], 2.4e9, env)
        y_raw *= pwr
        row.append(10 * log10(y_raw))
    stns.append(row)


line = np.polynomial.polynomial.Polynomial.fit(xs, ys, 1)
ly = []
for x in xs:
    ly.append(line(x))

plt.plot(xs, ys)
plt.plot(xs, ly)
plt.show()