import numpy as np
from PIL import Image
from math import floor, sqrt, log10, log2, sin, cos, acos, pi
import loadmap as me
import rt_slope
import cmath

CORNER_THRESHOLD = 30
CORNER_LOSS = 40

def sign(val):
    if val == 0:
        return 1
    else:
        return val / abs(val)

def has_los(p1:tuple, p2:tuple, map, scale):
    px1 = (floor(p1[0] / scale), floor(p1[1] / scale))
    px2 = (floor(p2[0] / scale), floor(p2[1] / scale))

    rise = px2[1] - px1[1]
    run = px2[0] - px1[0]

    if rise == 0 and run == 0:
        return True
    elif (not (rise == 0 or run == 0)):
        slope = rise / run
        itcpt = px1[1] - slope * px1[0]

    if rise == 0:
        y_loc = px1[1]
        for i in range(0,floor(run), round(sign(run))):
            x_loc = px1[0] + i
            if map[x_loc, y_loc] == 0:
                return False
    elif run == 0:
        x_loc = px1[0]
        for i in range(0, floor(rise), round(sign(rise))):
            y_loc = px1[1] + i
            if map[x_loc, y_loc] == 0:
                return False
    elif abs(rise) > abs(run):
        for i in range(0,floor(rise), round(sign(rise))):
            y_loc = px1[1] + i
            x_loc = round((y_loc - itcpt) / slope)
            if map[x_loc, y_loc] == 0:
                return False
    else:
        for i in range(0, floor(run), round(sign(run))):
            x_loc = px1[0] + i
            y_loc = round(slope * x_loc + itcpt)
            if map[x_loc, y_loc] == 0:
                return False
    return True

def create_mine(size):
    img = np.zeros((size,size))

    for x in range(1,size-1):
        for y in range(1,size-1):
            if x % 2 == 1 or y % 2 == 1:
                img[x,y] = 1
    
    return img

def gen_rf_props(width, height, eps_wall, eps_ceil, sig_wall, sig_ceil):
    print("Estimating path loss constant:")
    slope_tun = me.Tunnel(me.Point(0, -2), me.Point(0, 602), width, height, eps_wall, sig_wall, eps_ceil, sig_ceil)

    slope, intercept = rt_slope.get_rtslope(slope_tun, 2.4e9, height/2)
    print("Path Loss Line: {}*d + {}".format(slope, intercept))
    return (slope, intercept)

def draw_binary_cmap(env, loc, c_range, scale):
    c_map = np.zeros(np.shape(env))

    for x in range(0,size):
        for y in range(0,size):
            if env[x,y] == 1:
                los = has_los((loc[0], loc[1]), ((x+0.5)*scale, (y+0.5)*scale), env, scale)
                if los and (floor(loc[0] / scale) == x or floor(loc[1] / scale) == y):
                    inrange = sqrt(pow(loc[0] - (x+0.5)*scale, 2) + pow(loc[1] - (y+0.5)*scale, 2)) < c_range[0]
                    if inrange:
                        c_map[x][y] = 1
                else:
                    # Determine if the point is within the diamond
                    dmnd = [
                        (loc[0], loc[1] + 2*c_range[1]),
                        (loc[0] + 2*c_range[1], loc[1]),
                        (loc[0], loc[1] - 2*c_range[1]),
                        (loc[0] - 2*c_range[1], loc[1])
                        ]
                    
                    x_upscale = (x + 0.5) * scale
                    y_upscale = (y + 0.5) * scale

                    area_tri = 0
                    for i in range(0,4):
                        j = (i + 1) % 4
                        side_a = sqrt((dmnd[i][0] - x_upscale)**2 + (dmnd[i][1] - y_upscale)**2)
                        side_b = sqrt((dmnd[j][0] - x_upscale)**2 + (dmnd[j][1] - y_upscale)**2)
                        side_c = sqrt((dmnd[i][0] - dmnd[j][0])**2 + (dmnd[i][1] - dmnd[j][1])**2)
                        # if side_a < 0 and abs(side_a) < 0.000001:
                        #     side_a = 0
                        # if side_b < 0 and abs(side_b) < 0.000001:
                        #     side_b = 0
                        # if side_c < 0 and abs(side_c) < 0.000001:
                        #     side_c = 0
                        sp = (side_a + side_b + side_c) / 2
                        area_tri += sqrt(sp * (sp - side_a) * (sp - side_b) * (sp - side_c))
                        
                    area_rect = 8*c_range[1]**2

                    if area_tri <= area_rect + 0.000001:
                        c_map[x][y] = 1
    return c_map

def euclidean(p1, p2):
    return sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

def manhattan(p1, p2):
    return abs(p2[0] - p1[0]) + abs(p2[1] - p1[1])

def calc_approx_loss(n1, n2, rf_prop:tuple, map, scale):
    # p1 = me.Point(n1[0],n1[1])
    # p2 = me.Point(n2[0],n2[1])

    if has_los(n1, n2, map, scale):
        # dist = p1.distance(p2)
        dist = euclidean(n1, n2)
        return dist * rf_prop[0] + rf_prop[1]
    else:
        dist = manhattan(n1, n2)
        # dist = abs(p2.x - p1.x) + abs(p2.y - p1.y)
        loss = dist * rf_prop[0] + rf_prop[1]
        c1 = (n1[0], n2[1])
        c2 = (n2[0], n1[1])
        los_c1 = has_los(c1, n2, map, scale) and has_los(c1, n1, map, scale)
        los_c2 = has_los(c2, n2, map, scale) and has_los(c2, n1, map, scale)
        if not (los_c1 or los_c2):
            return float("-inf")
        elif (los_c1 and los_c2):
            if euclidean(n2, c1) > euclidean(n2, c2):
                c = c2
            else:
                c = c1
        elif (los_c1):
            c = c1
        elif (los_c2):
            c = c2
        dc = euclidean(n2,c)
        if dc > CORNER_THRESHOLD:
            return loss - CORNER_LOSS
        else:
            cl = (dc / CORNER_THRESHOLD) * CORNER_LOSS
            return loss - cl

def draw_approx_cmap(mine, scale, loc, height, tx_pow, noise, onehop_thpt, rf_prop):
    c_map = np.zeros(np.shape(mine))
    for x in range(0,len(c_map)):
        for y in range(0,len(c_map[x])):
            if (mine[x][y] == 1):
                loss = calc_approx_loss(((x+0.5)*scale, (y+0.5)*scale, height), (loc[0], loc[1], height), rf_prop, mine, scale)
                if (loss == float("-inf")):
                    thpt = 0
                else:
                    rx_pow = tx_pow + loss
                    rx_snr = rx_pow - noise
                    thpt = bw * log2(1 + pow(10,(rx_snr / 10)))# / pow(2, n.hops)
                if thpt < onehop_thpt:
                    thpt = 0
                c_map[x][y] = thpt
    return c_map

def simple_raytrace_loss(dist, a_height, f, width, height, eps_ceil, eps_wall, sig_ceil, sig_wall):
    global CORNER_THRESHOLD
    global CORNER_LOSS
    # global NUM_MODES
    NUM_MODES = 25

    c = 299792458
    eps_zero = 8.854187817 * pow(10, -12)
    k = 2 * pi * f / c

    if dist == float('inf') or dist == float('nan'):
        return 0.0

    a = width / 2
    b = height / 2

    eps_ceil = eps_ceil - (1j * sig_ceil / (2*pi*f*eps_zero))
    eps_wall = eps_wall - (1j * sig_wall / (2*pi*f*eps_zero))

    tx_z = a_height - (height / 2)
    rx_z = a_height - (height / 2)
    
    sum = 0
    for m in range(-NUM_MODES, NUM_MODES + 1):
        for n in range(-NUM_MODES, NUM_MODES + 1):
            x_m = (2 * m * a) #+ (pow(-1, m) * tx_x)
            z_n = (2 * n * b) + (pow(-1, n) * tx_z)

            # sqrt(pow(rx_x - x_m, 2) + pow(d_tot, 2) + pow(rx_z - z_n, 2))
            r_mn = sqrt(pow(x_m, 2) + pow(dist, 2) + pow(rx_z - z_n, 2))

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

def calc_raytrace_loss(n1, n2, map, scale, a_height, f, t_width, t_height, eps_ceil, eps_wall, sig_ceil, sig_wall):
    if has_los(n1, n2, map, scale):
        dist = euclidean(n1, n2)
        return 10*log10(simple_raytrace_loss(dist, a_height, f, t_width, t_height, eps_ceil, eps_wall, sig_ceil, sig_wall))
    else:
        dist = manhattan(n1, n2)
        loss = 10*log10(simple_raytrace_loss(dist, a_height, f, t_width, t_height, eps_ceil, eps_wall, sig_ceil, sig_wall))
        c1 = (n1[0], n2[1])
        c2 = (n2[0], n1[1])
        los_c1 = has_los(c1, n2, map, scale) and has_los(c1, n1, map, scale)
        los_c2 = has_los(c2, n2, map, scale) and has_los(c2, n1, map, scale)
        if not (los_c1 or los_c2):
            return float("-inf")
        elif (los_c1 and los_c2):
            if euclidean(n2, c1) > euclidean(n2, c2):
                c = c2
            else:
                c = c1
        elif (los_c1):
            c = c1
        elif (los_c2):
            c = c2
        dc = euclidean(n2,c)
        if dc > CORNER_THRESHOLD:
            return loss - CORNER_LOSS
        else:
            cl = (dc / CORNER_THRESHOLD) * CORNER_LOSS
            return loss - cl

def draw_raytrace_cmap(mine, scale, loc, height, tx_pow, noise, onehop_thpt, width, tun_height, eps_ceil, eps_wall, sig_ceil, sig_wall):
    c_map = np.zeros(np.shape(mine))
    for x in range(0,len(c_map)):
        for y in range(0,len(c_map[x])):
            if (mine[x][y] == 1):
                loss = calc_raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (loc[0], loc[1], height), mine, scale, height, 2.4e9, width, tun_height, eps_ceil, eps_wall, sig_ceil, sig_wall)
                # loss = rt.raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), 2.4e9, xml_env)
                if (loss == float("-inf")):
                    thpt = 0
                else:
                    # loss = 10*log10(loss)
                    rx_pow = tx_pow + loss
                    rx_snr = rx_pow - noise
                    if rx_snr < 0:
                        thpt = 0
                    else:
                        thpt = rx_snr
                    # thpt = bw * log2(1 + pow(10,(rx_snr / 10)))# / pow(2, n.hops)
                # if thpt < onehop_thpt:
                if thpt < 45:
                    thpt = 0
                c_map[x][y] = thpt
    return c_map

def cmap_to_image(env, cmap, node_locs):
    cmap_img = cmap / np.max(cmap)
    cmap_inv = np.multiply(1 - cmap_img, env)
    image = np.dstack((cmap_inv, cmap_img, np.zeros(np.shape(cmap_img))))

    for n in node_locs:
        image[floor(n[0] / scale),floor(n[1] / scale),0] = 0
        image[floor(n[0] / scale),floor(n[1] / scale),1] = 0
        image[floor(n[0] / scale),floor(n[1] / scale),2] = 1

    image *= 255
    return Image.fromarray(np.uint8(image.swapaxes(0,1)))

if __name__ == "__main__":
    scale = 4
    noise = -120 # dBm
    tx_pow = 3 # dBm
    bw = 20 * 10**6 # Hz
    thpt_thresh = 300*10**6

    size = 43
    img = create_mine(size)

    # Image.fromarray(img.transpose() * 255).show()
    center = ((size / 2) * scale, (size / 2) * scale)

    # eps_ceil = "8.9" sig_ceil = "0.15" eps_wall = "8.9" sig_wall = "0.15"
    rf_props = gen_rf_props(4, 2, 8.9, 8.9, 0.15, 0.15)
    comm_thresh = noise + tx_pow + 10*log10(pow(2, thpt_thresh/bw) - 1)
    los_range = (comm_thresh - rf_props[1]) / rf_props[0]
    dmd_range = 17

    # binary_cmap = draw_binary_cmap(img, center, (los_range, dmd_range), scale)
    # cmap_to_image(img, binary_cmap, [center]).show()
    # cmap_to_image(img, binary_cmap, [center]).save("node_coverage_binary.png")

    # approx_cmap = draw_approx_cmap(img, scale, center, 1, tx_pow, noise, thpt_thresh, (rf_props[0] * 1.2, rf_props[1] + 3.3))
    # cmap_to_image(img, approx_cmap, [center]).show()
    # cmap_to_image(img, approx_cmap, [center]).save("node_coverage_approx.png")

    raytrace_cmap = draw_raytrace_cmap(img, scale, center, 1, tx_pow, noise, thpt_thresh, 4, 2, 8.9, 8.9, 0.15, 0.15)
    # cmap_to_image(img, raytrace_cmap, [center]).show()
    # cmap_to_image(img, raytrace_cmap, [center]).save("node_coverage_raytrace_snr.png")

    raytrace_solidified = np.copy(raytrace_cmap)

    for x in range(len(raytrace_solidified)):
        for y in range(len(raytrace_solidified[0])):
            if raytrace_solidified[x,y] != 0:
                raytrace_solidified[x, y] = 1
    
    cmap_to_image(img, raytrace_solidified, [center]).show()
    cmap_to_image(img, raytrace_solidified, [center]).save("node_coverage_rtsolid_snr.png")

