import raytrace as rtloss
import loadmap as me
import numpy as np
from math import log10, sqrt, sin, cos, acos, pi
import cmath
from scipy import stats
from matplotlib import pyplot as plt

NUM_MODES = 25

def rt_tunnel(t:me.Tunnel, tx, rx, f):
    global NUM_MODES

    c = 299792458
    eps_zero = 8.854187817 * pow(10, -12)
    k = 2 * pi * f / c

    tx_pt = me.Point(tx[0], tx[1])
    rx_pt = me.Point(rx[0], rx[1])

    # d_tot, d_txc, d_rxc, ttx, trx = env.get_distance(tx_pt, rx_pt)
    d_tot = tx_pt.distance(rx_pt)

    if d_tot == float('inf') or d_tot == float('nan'):
        return 0.0

    a = t.width / 2
    b = t.height / 2

    eps_ceil = t.eps_ceil - (1j * t.sig_ceil / (2*pi*f*eps_zero))
    eps_wall = t.eps_wall - (1j * t.sig_wall / (2*pi*f*eps_zero))

    tx_rel = t.get_relative_coordinate(tx_pt)
    rx_rel = t.get_relative_coordinate(rx_pt)

    tx_x = tx_rel.x
    # tx_y = 0
    tx_z = (-t.height / 2) + tx[2]

    rx_x = rx_rel.x
    # rx_y = d_tot
    rx_z = (-t.height / 2) + rx[2]
    
    sum = 0
    for m in range(-NUM_MODES, NUM_MODES + 1):
        for n in range(-NUM_MODES, NUM_MODES + 1):
            x_m = (2 * m * a) + (pow(-1, m) * tx_x)
            z_n = (2 * n * b) + (pow(-1, n) * tx_z)

            r_mn = sqrt(pow(rx_x - x_m, 2) + pow(d_tot, 2) + pow(rx_z - z_n, 2))

            if r_mn == 0:
                return 1

            theta_perp = acos(abs(x_m - rx_x) / r_mn)
            theta_pll = acos(abs(z_n - rx_z) / r_mn)

            Delta_perp = cmath.sqrt(eps_wall - pow(sin(theta_perp), 2))
            Delta_pll = cmath.sqrt(eps_ceil - pow(sin(theta_pll), 2)) / eps_ceil

            rho_perp = (cos(theta_perp) - Delta_perp) / (cos(theta_perp) + Delta_perp)
            rho_pll = (cos(theta_pll) - Delta_pll) / (cos(theta_pll) + Delta_pll)

            sum += ((cmath.exp(-1j * k * r_mn) / r_mn) * pow(rho_perp,abs(m)) * pow(rho_pll,abs(n)))

    loss_pwr = pow((c/f) / (4 * pi), 2) * pow(abs(sum), 2)
    return loss_pwr

def get_rtslope(t:me.Tunnel, freq, height):
    # tx = (t.width/2, 0, height)
    tx = (0, 0, height)
    ys = np.linspace(1,600,1000)
    pows = []
    for y in ys:
        # rx = (t.width / 2, y, height)
        rx = (0, y, height)
        pows.append(10*log10(rt_tunnel(t, tx, rx, freq)))
    
    # plt.plot(ys,pows)
    slope, intercept, r_value, p_value, std_err = stats.linregress(ys,pows)
    # plt.plot(ys, ys * slope + intercept)
    # plt.title("Full Model vs Linear Approximation")
    # plt.ylabel("Loss (dB)")
    # plt.xlabel("Distance (meters)")
    # plt.legend(["Full Model", "Approximation"])
    # plt.savefig("full-vs-approx.png")
    # plt.show()
    return (slope, intercept)


# ys = np.linspace(1,600,1000)

# env = me.Environment()
# env.load("single_tunnel.xml")

# slopes = []

# f = 2.45e9
# h = 1
# values = []
# for y in ys:
#     loss = rtloss.raytrace_loss((0, y, .1), (0, 0, .1), f, env)
#     loss = 10*log10(loss)
#     values.append(loss)

# slope, intercept, r_value, p_value, std_err = stats.linregress(ys,values)

# rud = (slope * ys) + intercept
# plt.plot(ys, values, lw=0.3)
# plt.plot(ys, rud)
# plt.show()