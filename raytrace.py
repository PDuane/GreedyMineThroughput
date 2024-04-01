import loadmap as me
from math import pi, sqrt, acos, sin, cos, exp
import cmath

CORNER_THRESHOLD = 30
CORNER_LOSS = 30
NUM_MODES = 25

# tx, rx in a tuple of (x, y, z) where z is the height
def raytrace_loss(tx:tuple, rx:tuple, f, env:me.Environment):
    global CORNER_THRESHOLD
    global CORNER_LOSS
    global NUM_MODES

    c = 299792458
    eps_zero = 8.854187817 * pow(10, -12)
    k = 2 * pi * f / c

    tx_pt = me.Point(tx[0], tx[1])
    rx_pt = me.Point(rx[0], rx[1])

    d_tot, d_txc, d_rxc, ttx, trx = env.get_distance(tx_pt, rx_pt)

    if d_tot == float('inf') or d_tot == float('nan'):
        return 0.0

    a = ttx.width / 2
    b = ttx.height / 2

    eps_ceil = ttx.eps_ceil - (1j * ttx.sig_ceil / (2*pi*f*eps_zero))
    eps_wall = ttx.eps_wall - (1j * ttx.sig_wall / (2*pi*f*eps_zero))

    tx_rel = ttx.get_relative_coordinate(tx_pt)
    rx_rel = trx.get_relative_coordinate(rx_pt)

    tx_x = tx_rel.x
    # tx_y = 0
    tx_z = (-ttx.height / 2) + tx[2]

    rx_x = rx_rel.x
    # rx_y = d_tot
    rx_z = (-trx.height / 2) + rx[2]

    if env.has_los(tx_pt, rx_pt):
        corner_loss = 1
    else:
        if d_txc < CORNER_THRESHOLD:
            ratio = d_txc / CORNER_THRESHOLD
            corner_loss = ratio * CORNER_LOSS
        else:
            corner_loss = CORNER_LOSS
        
        corner_loss = pow(10, -corner_loss / 10)
    
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

    loss_pwr = pow((c/f) / (4 * pi), 2) * pow(abs(sum) * corner_loss, 2)

    return loss_pwr

# dielectric parameters should be (epsilon, sigma)
# Since the environment is not provided, this assumes line of sight
def raytrace_loss_raw(tx:tuple, rx:tuple, f:float, width:float, height:float, de_ceil:tuple, de_wall:tuple):

    global CORNER_THRESHOLD
    global CORNER_LOSS
    global NUM_MODES

    c = 299792458
    eps_zero = 8.854187817 * pow(10, -12)
    k = 2 * pi * f / c

    # tx_pt = me.Point(tx[0], tx[1])
    # rx_pt = me.Point(rx[0], rx[1])

    # d_tot, d_txc, d_rxc, ttx, trx = env.get_distance(tx_pt, rx_pt)

    d_tot = abs(rx[1] - tx[1])

    if d_tot == float('inf') or d_tot == float('nan'):
        return 0.0

    a = width / 2
    b = height / 2

    eps_ceil = de_ceil[0] - (1j * de_ceil[1] / (2*pi*f*eps_zero))
    eps_wall = de_wall[0] - (1j * de_wall[1] / (2*pi*f*eps_zero))
    
    sum = 0
    for m in range(-NUM_MODES, NUM_MODES + 1):
        for n in range(-NUM_MODES, NUM_MODES + 1):
            x_m = (2 * m * a) + (pow(-1, m) * tx[0])
            z_n = (2 * n * b) + (pow(-1, n) * tx[2])

            r_mn = sqrt(pow(rx[0] - x_m, 2) + pow(d_tot, 2) + pow(rx[2] - z_n, 2))

            if r_mn == 0:
                return 1

            theta_perp = acos(abs(x_m - rx[0]) / r_mn)
            theta_pll = acos(abs(z_n - rx[2]) / r_mn)

            Delta_perp = cmath.sqrt(eps_wall - pow(sin(theta_perp), 2))
            Delta_pll = cmath.sqrt(eps_ceil - pow(sin(theta_pll), 2)) / eps_ceil

            rho_perp = (cos(theta_perp) - Delta_perp) / (cos(theta_perp) + Delta_perp)
            rho_pll = (cos(theta_pll) - Delta_pll) / (cos(theta_pll) + Delta_pll)

            sum += ((cmath.exp(-1j * k * r_mn) / r_mn) * pow(rho_perp,abs(m)) * pow(rho_pll,abs(n)))

    loss_pwr = pow((c/f) / (4 * pi), 2) * pow(abs(sum), 2)

    return loss_pwr