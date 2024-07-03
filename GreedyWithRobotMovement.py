from PIL import Image
from math import pi, sin, cos, acos, sqrt, log2, log10, floor
import cmath
import numpy as np
import gen_map
from matplotlib import pyplot as plt
import loadmap as me
import rt_slope
import raytrace as rt
import sys
import time
import os.path

CORNER_THRESHOLD = 30
CORNER_LOSS = 30

class Node:
    location:tuple
    x:float
    y:float
    h:float
    connected = []
    hops:int
    rate:float

    def __init__(self, location, hops, rate) -> None:
        self.location = location
        self.x = location[0]
        self.y = location[1]
        self.h = location[2]
        self.hops = hops
        self.rate = rate

def euclidean(p1, p2):
    return sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

def manhattan(p1, p2):
    return abs(p2[0] - p1[0]) + abs(p2[1] - p1[1])

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

# rf_prop should be a tuple of (slope, intercept) for the rf signal strength
#     in dB
# Returns result in dB
# def calc_loss(n1, n2, rf_prop:tuple, map, scale, env:me.Environment):
def calc_loss(n1, n2, map, scale, a_height, f, t_width, t_height, eps_ceil, eps_wall, sig_ceil, sig_wall):

    # rt.raytrace_loss(n1, n2, 2.4e9. env)

    if has_los(n1, n2, map, scale):
        dist = euclidean(n1, n2)
        # return dist * rf_prop[0] + rf_prop[1]
        return 10*log10(simple_raytrace_loss(dist, a_height, f, t_width, t_height, eps_ceil, eps_wall, sig_ceil, sig_wall))
    else:
        dist = manhattan(n1, n2)
        # loss = dist * rf_prop[0] + rf_prop[1]
        loss = 10*log10(simple_raytrace_loss(dist, a_height, f, t_width, t_height, eps_ceil, eps_wall, sig_ceil, sig_wall))
        c1 = (n1[0], n2[1])
        c2 = (n2[0], n1[1])
        los_c1 = has_los(c1, n2, map, scale) and has_los(c1, n1, map, scale)
        los_c2 = has_los(c2, n2, map, scale) and has_los(c2, n1, map, scale)
        if not (los_c1 or los_c2):
            return float("-inf")
        elif (los_c1):
            c = c1
        elif (los_c2):
            c = c2
        else:
            if euclidean(n1, c1) > euclidean(n1, c2):
                c = c2
            else:
                c = c1
        dc = euclidean(n1,c)
        if dc > CORNER_THRESHOLD:
            return loss - CORNER_LOSS
        else:
            cl = (dc / CORNER_THRESHOLD) * CORNER_LOSS
            return loss - cl

def simple_raytrace_loss(dist, a_height, f, width, height, eps_ceil, eps_wall, sig_ceil, sig_wall):
    global CORNER_THRESHOLD
    global CORNER_LOSS
    # global NUM_MODES
    NUM_MODES = 25

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

def pathfind_bfs(env_map, start, dest):
    env = np.zeros(np.shape(env_map))
    for x in range(len(env)):
        for y in range(len(env[x])):
            env[x,y] = env_map[x,y]
    queue = [[start.copy(), []]]
    visited = []
    while (len(queue) > 0):
        item = queue.pop(0)
        loc = item[0]
        path = item[1].copy()
        path.append(loc)
        visited.append(loc)
        if loc[0] == dest[0] and loc[1] == dest[1]:
            return path, len(path)
        else:
            options = [
                [loc[0] + 1, loc[1]],
                [loc[0] - 1, loc[1]],
                [loc[0], loc[1] + 1],
                [loc[0], loc[1] - 1]
                ]
            for o in options:

                if o[0] >= 0 and o[1] >= 0:
                    if o[0] < len(env) and o[1] < len(env[0]):
                        if env[o[0], o[1]] == 1:
                            inQueue = False
                            for item in queue:
                                if item[0][0] == o[0] and item[0][1] == o[1]:
                                    inQueue = True
                                    break
                            if not inQueue:
                                queue.append([o, path.copy()])

            env[loc[0], loc[1]] = 0
    return [], 0


def pathfind(env, start, dest, parent):
    
    if start[0] >= len(env) or start[1] >= len(env[0]):
        return [], -1
    if env[start[0], start[1]] == 0:
        return [], -1
    elif start[0] == dest[0] and start[1] == dest[1]:
        return s, 1

    options = [
        [start[0] + 1, start[1]],
        [start[0] - 1, start[1]],
        [start[0], start[1] + 1],
        [start[0], start[1] - 1]
        ]
    
    searches = []
    for o in options:
        if o[0] >= 0 and o[1] >= 0:
            if o[0] < len(env) and o[1] < len(env[0]):
                if env[o[0], o[1]] == 1:
                    if not(o[0] == parent[0] and o[1] == parent[1]):
                        searches.append(o)
    
    new_env = np.zeros(np.shape(env))
    for i in range(len(env)):
        for j in range(len(env[i])):
            new_env[i,j] = env[i,j]
    new_env[start[0], start[1]] = 0

    best_path = []
    best_length = -1
    if len(searches) > 0:
        for s in searches:
            if s[0] == dest[0] and s[1] == dest[1]:
                return [s, start], 1
            path, pathlen = pathfind(new_env, s, dest, start)
            if pathlen > 0 and (best_length < 0 or pathlen < best_length):
                best_path = path
                best_length = pathlen

        if len(best_path) > 0:
            best_path.append(start)
        return best_path, best_length + 1
    else:
        return [], -1

def sign(val):
    if val == 0:
        return 1
    else:
        return val / abs(val)

def circle_points(loc, r):
    x = r
    y = 0

    points = [[int(loc[0]+x), int(loc[1])]]

    if r > 0:
        points.append([int(loc[0]-x), int(loc[1])])
        points.append([int(loc[0]), int(loc[1]+x)])
        points.append([int(loc[0]), int(loc[1]-x)])

    P = 1-r

    while x > y:
        y += 1

        if P <= 0:
            P = P + 2*y + 1
        else:
            x -= 1
            P = P + 2*y - 2*x + 1
        
        if (x < y):
            break

        points.append([int(loc[0]+x), int(loc[1]+y)])
        points.append([int(loc[0]-x), int(loc[1]+y)])
        points.append([int(loc[0]+x), int(loc[1]-y)])
        points.append([int(loc[0]-x), int(loc[1]-y)])

        if x != y:
            points.append([int(loc[0]+y), int(loc[1]+x)])
            points.append([int(loc[0]-y), int(loc[1]+x)])
            points.append([int(loc[0]+y), int(loc[1]-x)])
            points.append([int(loc[0]-y), int(loc[1]-x)])

        return points

def detect_obstacle(location, vis_range, og_env, obs_env):
    edge = circle_points(location, vis_range)
    e_len = len(edge)
    for i in range(len(edge)):
        if edge[e_len - i - 1][0] < 0 or edge[e_len - i - 1][1] < 0:
            del edge[e_len - i - 1]
    
    for p in edge:
        x = p[0]
        y = p[1]

        run = abs(x - location[0])
        rise = abs(y - location[1])
        if not (run == 0 or rise == 0):
            slope = rise / run
            itcp = location[1] - (rise / run) * location[0]

        if x == location[0]:
            dx = x
            for dy in range(location[1], y, int(sign(y - location[1]))):
                if (dy < len(obs_env[0])):
                    if obs_env[dx, dy] == 0:
                        og_env[dx, dy] = 0
                        break
        elif y == location[1]:
            dy = y
            for dx in range(location[0], x, int(sign(x - location[0]))):
                if (dx < len(obs_env)):
                    if obs_env[dx, dy] == 0:
                        og_env[dx, dy] = 0
                        break
        elif rise > run:
            for py in range(floor(rise)):
                dy = py + location[1]
                dx = round((dy - itcp) / slope)
                if (dx < len(obs_env) and dy < len(obs_env[0])):
                    if obs_env[dx, dy] == 0:
                        og_env[dx, dy] = 0
                        break
        elif run > rise:
            for px in range(floor(run)):
                dx = px + location[0]
                dy = round(dx * slope + itcp)
                if (dx < len(obs_env) and dy < len(obs_env[0])):
                    if obs_env[dx, dy] == 0:
                        og_env[dx, dy] = 0
                        break

    return og_env

def predict_node(map, c_map, scale, dropped, tx_pow, noise, onehop_thpt, xml_env:me.Environment):

    tun:me.Tunnel = xml_env.tunnels[0]

    c_map = np.zeros(np.shape(map))
    for x in range(0,len(c_map)):
        for y in range(0,len(c_map[x])):
            if (mine[x][y] == 1):
                best_thpt = 0
                for n in dropped:
                    loss = calc_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), mine, scale, height, 2.4e9, tun.width, tun.height, tun.eps_ceil, tun.eps_wall, tun.sig_ceil, tun.sig_wall)
                    # loss = rt.raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), 2.4e9, xml_env)
                    if (loss == float("-inf")):
                        thpt = 0
                    else:
                        # loss = 10*log10(loss)
                        rx_pow = tx_pow + loss
                        rx_snr = rx_pow - noise
                        thpt = bw * log2(1 + pow(10,(rx_snr / 10))) / pow(2, n.hops)
                    if thpt > best_thpt:
                        best_thpt = thpt
                c_map[x][y] = best_thpt

    # -------------------------------------------------------------
    # Save image of new coverage map
    # -------------------------------------------------------------
    # cmap_img = c_map / np.max(c_map)
    # cmap_inv = np.multiply(1 - cmap_img, mine)

    # image = np.dstack((cmap_inv, cmap_img, np.zeros(np.shape(cmap_img))))

    # for n in dropped_nodes:
    #     image[floor(n.x / scale),floor(n.y / scale),0] = 0
    #     image[floor(n.x / scale),floor(n.y / scale),1] = 0
    #     image[floor(n.x / scale),floor(n.y / scale),2] = 1


    # image *= 255

    # img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
    # img_gen.show()
    # img_gen.save("progress_cmap-pre.png")
    # -------------------------------------------------------------

    area = np.sum(c_map)
    c_map_best = np.zeros(np.shape(c_map))
    for x in range(0,len(map)):
        for y in range(0,len(map[x])):
            if map[x,y] == 1:
                best_sig = 0
                best_node = None
                for n in dropped:
                    loss = calc_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), mine, scale, height, 2.4e9, tun.width, tun.height, tun.eps_ceil, tun.eps_wall, tun.sig_ceil, tun.sig_wall)
                    # loss = rt.raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), 2.4e9, xml_env)
                    if (loss == float("-inf")):
                        thpt = 0
                    else:
                        # loss = 10*log10(loss)
                        rx_pow = tx_pow + loss
                        rx_snr = rx_pow - noise
                        thpt = bw * log2(1 + pow(10,(rx_snr / 10))) / pow(2, n.hops)

                    if (thpt > best_sig):
                        best_sig = thpt
                        best_node = n
                
                if best_sig >= onehop_thpt:
                    c_map_new = np.zeros(np.shape(c_map))
                    for x2 in range(0,len(c_map_new)):
                        for y2 in range(0,len(c_map_new[x2])):
                            if (map[x2, y2] == 1):
                                # loss = rt.raytrace_loss(((x2+0.5)*scale, (y2+0.5)*scale, height), ((x+0.5)*scale, (y+0.5)*scale, height), 2.4e9, xml_env)
                                loss = calc_loss(((x2+0.5)*scale, (y2+0.5)*scale, height), ((x+0.5)*scale, (y+0.5)*scale, height), mine, scale, height, 2.4e9, tun.width, tun.height, tun.eps_ceil, tun.eps_wall, tun.sig_ceil, tun.sig_wall)
                                if (loss == float("-inf")):
                                    thpt = 0
                                else:
                                    # loss = 10*log10(loss)
                                    rx_pow = tx_pow + loss
                                    rx_snr = rx_pow - noise
                                    thpt = bw * log2(1 + pow(10,(rx_snr / 10))) / pow(2, best_node.hops + 1)
                                if (map[x2,y2] == 1):
                                    c_map_new[x2,y2] = max(thpt, c_map[x2,y2])
                
                    new_coverage = np.sum(c_map_new)
                    if new_coverage > area:
                        area = new_coverage
                        new_loc = Node(((x+0.5)*scale,(y+0.5)*scale, height), best_node.hops + 1, best_sig)
                        for x2 in range(len(c_map_new)):
                            for y2 in range(len(c_map_new[0])):
                                c_map_best[x2,y2] = c_map_new[x2,y2]
    return new_loc, c_map_best

if __name__ == "__main__":
    c1 = (7.5,7.5)
    c2 = (5+2+4*24,5+2+3*24)

    cases = [
        ["ComplexRP_Case1", (19,205), 4,  4,  True],
        ["ComplexRP_Case3", (19,205), 4,  4,  True],
        ["ComplexRP_Case1", (19,205), 4, 10,  True],
        ["ComplexRP_Case3", (19,205), 4, 10,  True],
        ["SimRig_Case1",    ( 7,101), 4,  4,  True],
        ["SimRig_Case2",    ( 7,101), 4,  4,  True],
        ["SimRig_Case1",    ( 7,101), 4, 10,  True],
        ["SimRig_Case2",    ( 7,101), 4, 10,  True]
    ]

    avg_thpt_fname = "./results/raytracing_coverage_averages.txt"
    avg_time_fname = "./results/raytracing_time_averages.txt"
    has_of_been_opened = False

    nd_time = 0
    nd_attempts = 0

    for c in cases:
        print("Running simulation for {}, {} nodes".format(c[0], c[3]))

        env = me.Environment()

        # env.load("minexml_test.xml")  # For squre mine
        # env.load("SimRig.xml")        # For Simulation Rig mine
        # env.load("ComplexRP.xml")         # For Complex room and pillar
        # env.load("Mine1.xml")         # For Complex abnormal mine
        env.load("{}.xml".format(c[0]))

        # rx_loc = Node((22,357,1), 0, float("inf"))  # For squre mine
        # rx_loc = Node((255,53,1), 0, float("inf"))  # For Simulation Rig mine
        # rx_loc = Node((19,205,1), 0, float("inf"))  # For Complex room and pillar
        # rx_loc = Node((22,357,1), 0, float("inf"))  # For Complex abnormal mine
        rx_loc = Node((c[1][0],c[1][1],1), 0, float("inf"))

        # scale = 4                     # For square and simulation rig mines
        # scale = 6                     # For comple mines
        scale = c[2]

        #   Temporarily disable obstacle discovery and demonstrate how it
        #   performs in a fully known mine
        if c[4]:
            mine = env.draw_basic_bitmap(scale)
        else:
            mine = env.draw_obstacle_bitmap(scale)
        unmod = np.copy(mine)
        obs_mine = env.draw_obstacle_bitmap(scale)

        env.obstacles = []

        # Define environment and channel parameters
        noise = -120 # dBm
        tx_pow = 3 # dBm
        bw = 20 * 10**6 # Hz

        # rx_loc = Node((22,357,1), 0, float("inf"))
        robot_loc = [floor(rx_loc.x / scale), floor(rx_loc.y / scale)]
        height = 1
        num_nodes = c[3]

        dropped_nodes = [rx_loc]

        cmap = np.zeros(np.shape(mine))

        timer1 = time.time()

        for i in range(num_nodes):
            timer2 = time.time()
            
            while True:
                sum_before_path = np.sum(mine)
                findNewPath = True
                while findNewPath:
                    findNewPath = False
                    timer3 = time.time()
                    new_node, cmap = predict_node(mine, cmap, scale, dropped_nodes, 3, noise, 350*10**6, env)
                    attempt_time = time.time() - timer3
                    nd_time += attempt_time
                    nd_attempts += 1
                    print("\t\tTook {}s to predict node location".format(attempt_time))
                    path, length = pathfind_bfs(mine, robot_loc, [floor(new_node.x / scale), floor(new_node.y / scale)])

                    for l in range(len(path)):
                        robot_loc = path[l]
                        sum_before = np.sum(mine)
                        new_env = detect_obstacle(robot_loc, 50, np.copy(mine), obs_mine)
                        if not np.sum(new_env) == sum_before:
                            obs_locs = []
                            for x in range(len(new_env)):
                                for y in range(len(new_env[x])):
                                    if mine[x,y] - new_env[x,y] == 1:
                                        obs_locs.append((x,y))
                                    mine[x,y] = new_env[x,y]

                            for o in obs_locs:
                                tun_width = env.tunnels[0].width / 2
                                x_map = (o[0] + 0.5) * scale
                                y_map = (o[1] + 0.5) * scale
                                new_obs = me.Obstacle([
                                    me.Point(x_map - tun_width, y_map - tun_width),
                                    me.Point(x_map + tun_width, y_map + tun_width),
                                    me.Point(x_map - tun_width, y_map + tun_width),
                                    me.Point(x_map + tun_width, y_map - tun_width)
                                ])
                                env.add_obstacle(new_obs)
                            for p in path:
                                if mine[p[0], p[1]] == 0 or l == len(path) - 1:
                                    findNewPath = True
                                    break
                                    
                            if findNewPath:
                                break
                                # timer3 = time.time()
                                # new_node, cmap = predict_node(mine, cmap, scale, dropped_nodes, 3, noise, 350*10**6, env)
                                # print("\t\tTook {}s to predict node location".format(time.time() - timer3))
                                # path, length = pathfind_bfs(mine, robot_loc, [floor(new_node.x / scale), floor(new_node.y / scale)])

                                # l = -1
                if np.sum(mine) == sum_before_path:
                    break
                
            dropped_nodes.append(new_node)
            robot_loc = [floor(new_node.x / scale), floor(new_node.y / scale)]
            print("\tTook {}s to drop node".format(time.time() - timer2))

            # -------------------------------------------------------------
            # Save image of new coverage map
            # -------------------------------------------------------------
            # cmap_img = cmap / np.max(cmap)
            # cmap_inv = np.multiply(1 - cmap_img, mine)

            # image = np.dstack((cmap_inv, cmap_img, np.zeros(np.shape(cmap_img))))

            # for n in dropped_nodes:
            #     image[floor(n.x / scale),floor(n.y / scale),0] = 0
            #     image[floor(n.x / scale),floor(n.y / scale),1] = 0
            #     image[floor(n.x / scale),floor(n.y / scale),2] = 1


            # image *= 255

            # img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
            # img_gen.show()
            # img_gen.save("progress_cmap-{}.png".format(i))
            # -------------------------------------------------------------

        nd_time_total = time.time() - timer1
        print("Took {}s to drop all nodes".format(time.time() - timer1))

        cmap = cmap / np.max(cmap)
        cmap_inv = np.multiply(1 - cmap, mine)

        image = np.dstack((cmap_inv, cmap, np.zeros(np.shape(cmap))))

        print("Node Locations: ")
        for n in dropped_nodes:
            print("\t({x},{y})".format(x=n.x, y=n.y))
            image[floor(n.x / scale),floor(n.y / scale),0] = 0
            image[floor(n.x / scale),floor(n.y / scale),1] = 0
            image[floor(n.x / scale),floor(n.y / scale),2] = 1


        image *= 255

        img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
        # img_gen.show()
        print("Saving observed map for raytracing for {}, {} nodes.png".format(c[0], c[3]))
        if c[4]:
            img_gen.save("results/ratracing_observed_{}-{}nodes.png".format(c[0], c[3]))
        else:
            img_gen.save("results/ratracing_known_{}-{}nodes.png".format(c[0], c[3]))


        # !!! Not needed when obstacle is fully known. Uncomment before     !!!
        # !!! running obstacle discovery                                    !!!
        # ---------------------------------------------------------------------
        if c[4]:
            # Show coverage of actual map based on the throughput estimation method
            print("Recalculating loss using known obstacles")
            tun:me.Tunnel = env.tunnels[0]
            
            print("\tDetermining Routing")
            for i in range(1, len(dropped_nodes)):
                dropped_nodes[i].hops = -1
            to_process = [rx_loc]
            while len(to_process) > 0:
                next_node = to_process.pop(0)
                for n in dropped_nodes:
                    if n.hops == -1:
                        loss = calc_loss((n.x, n.y, n.h), (next_node.x, next_node.y, next_node.h), obs_mine, scale, height, 2.4e9, tun.width, tun.height, tun.eps_ceil, tun.eps_wall, tun.sig_ceil, tun.sig_wall)
                        if (loss == float("-inf")):
                            thpt = 0
                        else:
                            # loss = 10*log10(loss)
                            rx_pow = tx_pow + loss
                            rx_snr = rx_pow - noise
                            thpt = bw * log2(1 + pow(10,(rx_snr / 10))) #/ pow(2, next_node.hops + 1)
                        
                        if thpt >= 350*10**6:
                            n.hops = next_node.hops + 1
                            n.rate = thpt / pow(2, next_node.hops + 1)
                            to_process.append(n)
            for i in range(1, len(dropped_nodes)):
                if dropped_nodes[i].hops == -1:
                    dropped_nodes[i].hops = float("inf")

            print("\tCalculating Coverage")
            cmap = np.zeros(np.shape(obs_mine))
            for x in range(0,len(cmap)):
                for y in range(0,len(cmap[x])):
                    if (obs_mine[x][y] == 1):
                        best_thpt = 0
                        for n in dropped_nodes:
                            loss = calc_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), obs_mine, scale, height, 2.4e9, tun.width, tun.height, tun.eps_ceil, tun.eps_wall, tun.sig_ceil, tun.sig_wall)
                            # loss = rt.raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), 2.4e9, xml_env)
                            if (loss == float("-inf")):
                                thpt = 0
                            else:
                                # loss = 10*log10(loss)
                                rx_pow = tx_pow + loss
                                rx_snr = rx_pow - noise
                                thpt = bw * log2(1 + pow(10,(rx_snr / 10))) / pow(2, n.hops)
                            if thpt > best_thpt:
                                best_thpt = thpt
                        cmap[x][y] = best_thpt
        # ---------------------------------------------------------------------
        
        print("Saving average throughput")
        if os.path.isfile(avg_thpt_fname):
            f = open(avg_thpt_fname, "a")
        else:
            f = open(avg_thpt_fname, "w")
        
        # Run name, # nodes, avg throughput
        f.write("{},nodes:{}nodes,discovery:{},{}\n".format(c[0],c[3],c[4], np.sum(cmap * obs_mine) / np.sum(obs_mine)))
        f.close()


        print("Saving average time")
        if os.path.isfile(avg_time_fname):
            f = open(avg_time_fname, "a")
        else:
            f = open(avg_time_fname, "w")
        
        # Run name, # nodes, average time, total time
        f.write("{},nodes:{},discovery:{},{},{}\n".format(c[0],c[3],c[4], nd_time / nd_attempts, nd_time_total))
        f.close()
        
        # !!! Not needed when obstacle is fully known. Uncomment before     !!!
        # !!! running obstacle discovery                                    !!!
        # ---------------------------------------------------------------------
        if c[4]:
            cmap = cmap / np.max(cmap)
            cmap_inv = np.multiply(1 - cmap, obs_mine)

            image = np.dstack((cmap_inv, cmap, np.zeros(np.shape(cmap))))

            for n in dropped_nodes:
                image[floor(n.x / scale),floor(n.y / scale),0] = 0
                image[floor(n.x / scale),floor(n.y / scale),1] = 0
                image[floor(n.x / scale),floor(n.y / scale),2] = 1


            image *= 255

            img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
            print("Saving actual map with raytracing coverage for {}, {} nodes.png".format(c[0], c[3]))
            img_gen.save("results/raytracing_actual_{}-{}nodes.png".format(c[0], c[3]))
        # ---------------------------------------------------------------------
