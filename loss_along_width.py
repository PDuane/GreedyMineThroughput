from math import pi, sin, cos, acos, sqrt, log10
import cmath
import numpy as np
import matplotlib.pyplot as plt

CORNER_THRESHOLD = 30
CORNER_LOSS = 40

def simple_raytrace_loss(n1, n2, f, width, height, eps_ceil, eps_wall, sig_ceil, sig_wall):
    global CORNER_THRESHOLD
    global CORNER_LOSS
    # global NUM_MODES
    NUM_MODES = 25

    c = 299792458
    eps_zero = 8.854187817 * pow(10, -12)
    k = 2 * pi * f / c

    # if dist == float('inf') or dist == float('nan'):
    #     return 0.0

    a = width / 2
    b = height / 2

    eps_ceil = eps_ceil - (1j * sig_ceil / (2*pi*f*eps_zero))
    eps_wall = eps_wall - (1j * sig_wall / (2*pi*f*eps_zero))

    tx_z = n1[2] - (height / 2)
    rx_z = n2[2] - (height / 2)
    
    sum = 0
    for m in range(-NUM_MODES, NUM_MODES + 1):
        for n in range(-NUM_MODES, NUM_MODES + 1):
            x_m = (2 * m * a) + (pow(-1, m) * n1[0])
            z_n = (2 * n * b) + (pow(-1, n) * tx_z)

            # sqrt(pow(rx_x - x_m, 2) + pow(d_tot, 2) + pow(rx_z - z_n, 2))
            r_mn = sqrt(pow(n2[0] - x_m, 2) + pow(n2[1] - n1[1], 2) + pow(rx_z - z_n, 2))

            if r_mn == 0:
                return 1

            theta_perp = acos(abs(x_m) / r_mn)
            theta_pll = acos(abs(z_n - rx_z) / r_mn)

            Delta_perp = cmath.sqrt(eps_wall - pow(sin(theta_perp), 2))
            Delta_pll = cmath.sqrt(eps_ceil - pow(sin(theta_pll), 2)) / eps_ceil

            rho_perp = (cos(theta_perp) - Delta_perp) / (cos(theta_perp) + Delta_perp)
            rho_pll = (cos(theta_pll) - Delta_pll) / (cos(theta_pll) + Delta_pll)

            sum += ((cmath.exp(-1j * k * r_mn) / r_mn) * pow(rho_perp,abs(m)) * pow(rho_pll,abs(n)))

    loss_pwr = pow((c/f) / (4 * pi), 2) * pow(abs(sum), 2)

    return loss_pwr


t_width = 4
t_height = 2

eps = 8.9
sig = 0.15

# tx = (0, 0, t_height / 2)
xs = np.linspace(-2,2,500)
ds = np.linspace(20, 120, 100)
pows = []
count = 1
for x in xs:
    print("\r{}/500            ".format(count), end="")
    t_pow = 0
    for d in ds:
        tx = (x, 0, t_height / 2)
        rx = (x, d, t_height / 2)
        t_pow += 10*log10(simple_raytrace_loss(tx, rx, 2.4e9, t_width, t_height, eps, eps, sig, sig))
    t_pow /= len(ds)
    pows.append(t_pow)
    count += 1

plt.plot(xs, pows)
plt.xlabel("X coordinate")
plt.ylabel("Average power (dB)")
plt.savefig("pwr_along_width.png")
plt.show()