from math import sin, cos, acos, sqrt, pi
import cmath
import loadmap as me

def calc_loss(n1, n2, rf_prop:tuple, map, scale):
    p1 = me.Point(n1[0],n1[1])
    p2 = me.Point(n2[0],n2[1])

    if has_los(n1, n2, map, scale):
        dist = p1.distance(p2)
        return raytrace_loss(dist, a_height, f, width, height, eps_ceil, eps_wall, sig_ceil, sig_wall)
    else:
        dist = abs(p2.x - p1.x) + abs(p2.y - p1.y)
        loss = dist * rf_prop[0] + rf_prop[1]
        c1 = me.Point(p1.x, p2.y)
        c2 = me.Point(p2.x, p1.y)
        los_c1 = has_los((c1.x, c1.y), (p2.x, p2.y), map, scale) and has_los((c1.x, c1.y), (p1.x, p1.y), map, scale)
        los_c2 = has_los((c2.x, c2.y), (p2.x, p2.y), map, scale) and has_los((c2.x, c2.y), (p1.x, p1.y), map, scale)
        if not (los_c1 or los_c2):
            return float("-inf")
        elif (los_c1):
            c = c1
        elif (los_c2):
            c = c2
        else:
            if p1.distance(c1) > p1.distance(c2):
                c = c2
            else:
                c = c1
        dc = p1.distance(c)
        if dc > CORNER_THRESHOLD:
            return loss - CORNER_LOSS
        else:
            cl = (dc / CORNER_THRESHOLD) * CORNER_LOSS
            return loss - cl

def raytrace_loss(dist, a_height, f, width, height, eps_ceil, eps_wall, sig_ceil, sig_wall):
    global CORNER_THRESHOLD
    global CORNER_LOSS
    global NUM_MODES

    c = 299792458
    eps_zero = 8.854187817 * pow(10, -12)
    k = 2 * pi * f / c

    tx_pt = me.Point(0, 0)
    rx_pt = me.Point(0, dist)

    # d_tot, d_txc, d_rxc, ttx, trx = env.get_distance(tx_pt, rx_pt)

    if dist == float('inf') or dist == float('nan'):
        return 0.0

    a = width / 2
    b = height / 2

    eps_ceil = eps_ceil - (1j * sig_ceil / (2*pi*f*eps_zero))
    eps_wall = eps_wall - (1j * sig_wall / (2*pi*f*eps_zero))

    # tx_rel = ttx.get_relative_coordinate(tx_pt)
    # rx_rel = trx.get_relative_coordinate(rx_pt)

    # tx_x = tx_rel.x
    # # tx_y = 0
    tx_z = a_height - (height / 2)

    # rx_x = rx_rel.x
    # # rx_y = d_tot
    rx_z = a_height - (height / 2)

    # if env.has_los(tx_pt, rx_pt):
    #     corner_loss = 1
    # else:
    #     if d_txc < CORNER_THRESHOLD:
    #         ratio = d_txc / CORNER_THRESHOLD
    #         corner_loss = ratio * CORNER_LOSS
    #     else:
    #         corner_loss = CORNER_LOSS
        
    #     corner_loss = pow(10, -corner_loss / 10)
    
    sum = 0
    for m in range(-NUM_MODES, NUM_MODES + 1):
        for n in range(-NUM_MODES, NUM_MODES + 1):
            x_m = (2 * m * a) #+ (pow(-1, m) * tx_x)
            z_n = (2 * n * b) + (pow(-1, n) * tx_z)

            r_mn = sqrt(pow(dist, 2) + pow(rx_z - z_n, 2))

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