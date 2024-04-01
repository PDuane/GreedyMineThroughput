from PIL import Image
from math import pi, sqrt, log2, floor
import numpy as np
import gen_map
from matplotlib import pyplot as plt
import loadmap as me
import rt_slope
import raytrace as rt
import sympy.geometry as geo

CORNER_THRESHOLD = 30
CORNER_LOSS = 30

def has_los(c1, c2, map):
    r_c1x = round(c1[0])
    r_c1y = round(c1[1])
    r_c2x = round(c2[0])
    r_c2y = round(c2[1])
    if (r_c1x == r_c2x) and (r_c1y == r_c2y):
        if map[r_c1x, r_c1y] == 1 and map[r_c2x, r_c2y] == 1:
            return True
    slope = 0
    try:
        slope = (c2[1] - c1[1]) / (c2[0] - c1[0])
    except ZeroDivisionError:
        x = r_c1x
        if (c1[1] < c2[1]):
            start = r_c1y
            stop = r_c2y
        else:
            start = r_c2y
            stop = r_c1y
        for y in range(start, stop):
            if (x >= 0 and y >= 0 and map[x,y] == 0):
                return False
        
    intercept = c1[1] - slope * c1[0]
    if (abs(slope) < 1):
        # xs = np.asarray(range(round(c1[0]), round(c2[0])))
        # ys = ((xs) * slope) - c1[1]
        # rnd = [round(v) for v in ys]
        # plt.scatter(xs, rnd)
        # plt.show()
        if (c1[0] < c2[0]):
            start = r_c1x
            stop = r_c2x
        else:
            start = r_c2x
            stop = r_c1x
        for x in range(start, stop):
            y = round(((x) * slope) + intercept)
            if (x >= 0 and y >= 0 and map[x,y] == 0):
                return False
        return True
    else:
        # ys = np.asarray(range(round(c1[1]), round(c2[1])))
        # xs = (ys + c1[0]) / slope
        # rnd = [round(v) for v in xs]
        # plt.scatter(rnd, ys)
        # plt.show()
        if (c1[1] < c2[1]):
            start = r_c1y
            stop = r_c2y
        else:
            start = r_c2y
            stop = r_c1y
        
        for y in range(start, stop):
            x = round((y - intercept) / slope)
            if (x >= 0 and y >= 0 and map[x][y] == 0):
                return False
        return True

    # for x in range(round(c1[0]), round(c2[0])):
    #     for y in range (round(c1[1]), round(c2[1])):
    #         if (map[x][y]) <= 0.5:
    #             return False

    # return True

def dist(c1, c2):
    return sqrt((c2[0]-c1[0])**2 + (c2[1]-c1[1])**2)

def connected(n1, n2, comm_dist, map):
    # if dist(n1, n2) < comm_dist and has_los(n1, n2, map):
    #     pass
    n1_int = (round(n1[0]), round(n1[1]))
    n2_int = (round(n2[0]), round(n2[1]))
    is_free = map[n1_int[0],n1_int[1]] == 1 and map[n2_int[0],n2_int[1]] == 1
    return is_free and dist(n1, n2) < comm_dist and has_los(n1, n2, map)

def has_los(p1:me.Point, p2:me.Point, map, scale):
    px1 = (floor(p1.x / scale), floor(p1.y / scale))
    px2 = (floor(p2.x / scale), floor(p2.y / scale))
    rise = px2[1] - px1[1]
    run = px2[0] - px1[0]

    if rise == 0 and run == 0:
        return True
    elif (not (rise == 0 or run == 0)):
        slope = rise / run
        itcpt = px1[1] - slope * px1[0]

    if rise == 0:
        y_loc = px1[1]
        for i in range(0,floor(run), round(run/abs(run))):
            x_loc = px1[0] + i
            if map[x_loc, y_loc] == 0:
                return False
    elif run == 0:
        x_loc = px1[0]
        for i in range(0, floor(rise), round(rise / abs(rise))):
            y_loc = px1[1] + i
            if map[x_loc, y_loc] == 0:
                return False
    elif abs(rise) > abs(run):
        for i in range(0,floor(rise), round(rise/abs(rise))):
            y_loc = px1[1] + i
            x_loc = round((y_loc - itcpt) / slope)
            if map[x_loc, y_loc] == 0:
                return False
    else:
        for i in range(0, floor(run), round(run / abs(run))):
            x_loc = px1[0] + i
            y_loc = round(slope * x_loc + itcpt)
            if map[x_loc, y_loc] == 0:
                return False
    return True

# rf_prop should be a tuple of (slope, intercept) for the rf signal strength
#     in dB
# Returns result in dB
def calc_loss(n1, n2, rf_prop:tuple, env:me.Environment, map, scale):
    p1 = me.Point(n1[0],n1[1])
    p2 = me.Point(n2[0],n2[1])
    # dist, dc1, dc2, t1, t2 = env.get_distance(p1, p2)
    # loss = dist * rf_prop[0] + rf_prop[1]
    if has_los(p1, p2, map, scale):
        dist = p1.distance(p2)
        return dist * rf_prop[0] + rf_prop[1]
    else:
        dist = abs(p2.x - p1.x) + abs(p2.y - p1.y)
        loss = dist * rf_prop[0] + rf_prop[1]
        c1 = me.Point(p1.x, p2.y)
        c2 = me.Point(p2.x, p1.y)
        los_c1 = has_los(c1, p2, map, scale) and has_los(c1, p1, map, scale)
        los_c2 = has_los(c2, p2, map, scale) and has_los(c2, p1, map, scale)
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
        
# def calc_loss_img(n1, n2, rf_prop:tuple, mine, scale):
#     p1 = geo.Point(n1[0],n1[1])
#     p2 = geo.Point(n2[0],n2[1])

#     rise = n2[1] - n1[1]
#     run = n2[0] - n1[0]

#     dist, dc1, dc2, t1, t2 = env.get_distance(p1, p2)
#     loss = dist * rf_prop[0] + rf_prop[1]
#     if env.has_los(p1, p2):
#         return loss
#     else:
#         if dc1 > CORNER_THRESHOLD:
#             return loss - CORNER_LOSS
#         else:
#             cl = (dc1 / CORNER_THRESHOLD) * CORNER_LOSS
#             return loss - cl

def connected_corner_cvg(map, c_map, scale, loc, num_nodes, tx_pow, rf_prop, noise, snr):
    dropped = 0
    d_loc = [loc]

    c_map_best = np.zeros((len(c_map),len(c_map[0])))
    for x2 in range(len(c_map)):
        for y2 in range(len(c_map[0])):
            c_map_best[x2,y2] = c_map[x2,y2]

    while dropped < num_nodes:
        print("Node {}".format(dropped + 1))
        area = np.sum(c_map)
        c_map_best = np.zeros((len(c_map),len(c_map[0])))
        for x in range(0,len(map)):
            for y in range(0,len(map[x])):
                if map[x,y] == 1:
                    con = False
                    best_node = None
                    best_sig = float("-inf")
                    for n in d_loc:
                        loss = calc_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n[0], n[1], height), rf_prop, env, map, scale)
                        rx_pow = tx_pow + loss
                        rx_snr = rx_pow - noise
                        if (rx_snr > best_sig):
                            best_sig = rx_snr
                            best_node = n
                        # if rx_pow - noise > snr:
                        #     con = True
                        #     break
                    # if not con:
                    #     continue
                    c_map_new = np.zeros((len(c_map),len(c_map[0])))
                    for x2 in range(0,len(c_map_new)):
                        for y2 in range(0,len(c_map_new[x2])):
                            if (map[x2, y2] == 1):
                                loss = calc_loss(((x2+0.5)*scale, (y2+0.5)*scale, height), ((x+0.5)*scale, (y+0.5)*scale, height), rf_prop, env, map, scale)
                                rx_pow = tx_pow + loss
                                rx_snr = rx_pow - noise
                                thpt = 20*10**6 * log2(pow(10,(rx_snr / 10)) + 1)
                                if (map[x2,y2] == 1):
                                    c_map_new[x2,y2] = max(thpt, c_map[x2,y2])
                    
                    new_coverage = np.sum(c_map_new)
                    if new_coverage > area:
                        area = new_coverage
                        new_loc = ((x+0.5)*scale,(y+0.5)*scale)
                        for x2 in range(len(c_map_new)):
                            for y2 in range(len(c_map_new[0])):
                                c_map_best[x2,y2] = c_map_new[x2,y2]
        d_loc.append(new_loc)
        for x2 in range(0,len(c_map_best)):
            for y2 in range(0,len(c_map_best[x2])):
                    c_map[x2][y2] = c_map_best[x2,y2]
        dropped += 1
    return d_loc, c_map

def greedy_connected_cvg(map, c_map, loc, num_nodes, comm_dist):
    dropped = 0
    d_loc = [loc]

    c_map_best = np.zeros((len(c_map),len(c_map[0])))
    for x2 in range(len(c_map)):
        for y2 in range(len(c_map[0])):
            c_map_best[x2,y2] = c_map[x2,y2]

    while dropped < num_nodes:
        print("Node {}".format(dropped + 1))
        area = np.sum(c_map)

        for x in range(0,len(map)):
            for y in range(0,len(map[x])):
                if map[x,y] == 1:
                    con = False
                    for n in d_loc:
                        if connected((x, y), n, comm_dist, map):
                            con = True
                            break
                    if not con:
                        continue
                    c_map_new = np.zeros((len(c_map),len(c_map[0])))
                    # for x2 in range(len(c_map)):
                    #     for y2 in range(len(c_map[0])):
                    #         c_map_new[x2,y2] = c_map[x2,y2]
                    # c_map_new = c_map.copy()
                    for x2 in range(0,len(c_map_new)):
                        for y2 in range(0,len(c_map_new[x2])):
                            if (map[x2,y2] == 1):
                                if connected((x2, y2), (x,y), comm_dist, map) or c_map[x2,y2] == 1:
                                    c_map_new[x2][y2] = 1
                    
                    # Image.fromarray(np.uint8(c_map_new.transpose() * 255), 'L').show()
                    # diff_map = np.subtract(c_map_new, c_map)
                    # Image.fromarray(np.uint8(diff_map.transpose() * 255), 'L').show()
                    new_coverage = np.sum(c_map_new)
                    if new_coverage > area:
                        area = new_coverage
                        new_loc = (x,y)
        d_loc.append(new_loc)
        for x2 in range(0,len(c_map_new)):
            for y2 in range(0,len(c_map_new[x2])):
                if connected((x2, y2), new_loc, comm_dist, map):
                    c_map[x2][y2] = 1
        dropped += 1
        c_map = c_map_new
        # Image.fromarray(np.uint8(c_map.transpose() * 255), 'L').show()
    return d_loc, c_map


c1 = (7.5,7.5)
c2 = (5+2+4*24,5+2+3*24)
scale = 4

# psg = [(5,0,9,5)]
# mine = gen_map.gen_map((128, 128), 4, 20, 20, 1, (5, 5, 5, 5), psg, [])
env = me.Environment()
env.load("minexml_test.xml")
mine = env.draw_bitmap(scale)
cmap = np.zeros(np.shape(mine))

# Image.fromarray(np.uint8(cmap.transpose() * 255), 'L').show()

tun_ref:me.Tunnel = env.tunnels[0]
slope_tun = me.Tunnel(me.Point(0, -2), me.Point(0, 602), tun_ref.width, tun_ref.height, tun_ref.eps_wall, tun_ref.sig_wall, tun_ref.eps_ceil, tun_ref.sig_ceil)

slope, intercept = rt_slope.get_rtslope(slope_tun, 2.4e9, tun_ref.height/2)
l = me.Line(me.Point(0, intercept), me.Point(600, intercept + 600 * slope))
comm_dist = l.get_x(-80, True)

# rx_pow = tx_pow * loss

noise = -120 # dBm
snr = 38 # dB
tx_pow = 3 # dBm

# has_los(c1, c2, mine)

robot_loc = (22, 357)
height = 1
num_nodes = 4

for x in range(0,len(cmap)):
    for y in range(0,len(cmap[x])):
        if (mine[x][y] == 1):
            loss = calc_loss((robot_loc[0], robot_loc[1], height), ((x+0.5)*scale, (y+0.5)*scale, height), (slope, intercept), env, mine, scale)
            rx_pow = tx_pow + loss
            rx_snr = rx_pow - noise
            thpt = 20*10**6 * log2(1 + pow(10,(rx_snr / 10)))
            cmap[x][y] = thpt

Image.fromarray(np.uint8((cmap / np.max(cmap)).transpose() * 255), 'L').show()

node_locs, new_map = connected_corner_cvg(mine, cmap, scale, robot_loc, num_nodes, 3, (slope,intercept), -120, 40)

# mine = mine.transpose()
# new_map = new_map.transpose()

cmap = cmap / np.max(cmap)

image = np.repeat(mine[:,:,np.newaxis], 3, axis=2)  # green
image[:,:,2] = np.subtract(image[:,:,2], cmap)      # blue
image[:,:,0] = np.subtract(image[:,:,2], cmap)      # red
for x in range(len(cmap)):
    for y in range(len(cmap[0])):
        if mine[x,y] == 1:
            image[x,y,0] = (1-cmap[x,y])

# image[:,:,0] = np.subtract(image[:,:,0], cmap)

print("Node Locations: ")
for n in node_locs:
    print("\t({x},{y})".format(x=n[0], y=n[1]))
    image[floor(n[0] / scale),floor(n[1] / scale),0] = 0
    image[floor(n[0] / scale),floor(n[1] / scale),1] = 0
    image[floor(n[0] / scale),floor(n[1] / scale),2] = 1


image *= 255

img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
img_gen.show()
img_gen.save("C:\\Users\\patri\\OneDrive\\Documents\\NIOSH\\NodePlacement\\omnet_img.png")