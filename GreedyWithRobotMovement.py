from PIL import Image
from math import pi, sqrt, log2, log10, floor
import numpy as np
import gen_map
from matplotlib import pyplot as plt
import loadmap as me
import rt_slope
import raytrace as rt
import sys
import time

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
    
    # def connect(self, n:Node):
        


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
def calc_loss(n1, n2, rf_prop:tuple, map, scale, env:me.Environment):

    rt.raytrace_loss(n1, n2, 2.4e9. env)

    if has_los(n1, n2, map, scale):
        dist = euclidean(n1, n2)
        return dist * rf_prop[0] + rf_prop[1]
    else:
        dist = manhattan(n1, n2)
        loss = dist * rf_prop[0] + rf_prop[1]
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

def pathfind_bfs(env_map, start, dest):
    env = np.zeros(np.shape(env_map))
    for x in range(len(env)):
        for y in range(len(env[x])):
            env[x,y] = env_map[x,y]
    queue = [[start.copy(), []]]
    while (len(queue) > 0):
        item = queue.pop(0)
        loc = item[0]
        path = item[1].copy()
        path.append(loc)
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
                            # if not(o[0] == parent[0] and o[1] == parent[1]):
                            queue.append([o, path.copy()])
            
            # new_env = np.zeros(np.shape(env))
            # for i in range(len(env)):
            #     for j in range(len(env[i])):
            #         new_env[i,j] = env[i,j]
            env[loc[0], loc[1]] = 0
    return [], 0


def pathfind(env, start, dest, parent):
    
    # if start[0] == dest[0] and start[1] == dest[1]:
    #     return [dest], 1
    # else:
    #     env[start[0],start[1]] = 1
    
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
            # if s[0] < len(env) and s[1] < len(env[0]) and env[s[0], s[1]]  == 1:
            #     found += 1
            if s[0] == dest[0] and s[1] == dest[1]:
                return [s, start], 1
            # else:
            path, pathlen = pathfind(new_env, s, dest, start)
            if pathlen > 0 and (best_length < 0 or pathlen < best_length):
                # found += 1
                best_path = path
                best_length = pathlen

        if len(best_path) > 0:
            best_path.append(start)
        # if length(best_path) == 0:
        #     return [], -1
        # else:
        return best_path, best_length + 1
    else:
        return [], -1
    
    # if found > 0:
    #     best_path.append(start)
    #     return best_path, best_length + 1
    # else:
    #     return [], -1

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
                # dy = py + location[1]
                if (dy < len(obs_env[0])):
                    if obs_env[dx, dy] == 0:
                        og_env[dx, dy] = 0
                        break
        elif y == location[1]:
            dy = y
            for dx in range(location[0], x, int(sign(x - location[0]))):
                # dx = px + location[0]
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

    # if not c_map:
    c_map = np.zeros(np.shape(map))
    for x in range(0,len(c_map)):
        for y in range(0,len(c_map[x])):
            if (mine[x][y] == 1):
                best_thpt = 0
                for n in dropped:
                    loss = rt.raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), 2.4e9, xml_env)
                    if (loss <= 0):
                        thpt = 0
                    else:
                        loss = 10*log10(loss)
                    # loss = calc_loss((loc.x, loc.y, loc.h), (loc.x, loc.y, loc.h), (slope, intercept), mine, scale)
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

    
    # Image.fromarray(np.uint8((c_map / np.max(c_map)).transpose() * 255), 'L').show()

    # print("Node {}".format(num_dropped + 1))
    area = np.sum(c_map)
    c_map_best = np.zeros(np.shape(c_map))
    for x in range(0,len(map)):
        for y in range(0,len(map[x])):
            if map[x,y] == 1:
                best_sig = 0
                best_node = None
                for n in dropped:
                    loss = rt.raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), 2.4e9, xml_env)
                    # loss = calc_loss(((x+0.5)*scale, (y+0.5)*scale, height), (n.x, n.y, n.h), rf_prop, map, scale)
                    if (loss <= 0):
                        thpt = 0
                    else:
                        loss = 10*log10(loss)
                    # loss = calc_loss((loc.x, loc.y, loc.h), (loc.x, loc.y, loc.h), (slope, intercept), mine, scale)
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
                                loss = rt.raytrace_loss(((x2+0.5)*scale, (y2+0.5)*scale, height), ((x+0.5)*scale, (y+0.5)*scale, height), 2.4e9, xml_env)
                                # loss = calc_loss(((x2+0.5)*scale, (y2+0.5)*scale, height), ((x+0.5)*scale, (y+0.5)*scale, height), rf_prop, map, scale)
                                if (loss <= 0):
                                    thpt = 0
                                else:
                                    loss = 10*log10(loss)
                                # loss = calc_loss((loc.x, loc.y, loc.h), (loc.x, loc.y, loc.h), (slope, intercept), mine, scale)
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
    # d_loc.append(new_loc)
    # for x2 in range(len(c_map_new)):
    #     for y2 in range(len(c_map_new[0])):
    #         c_map[x2,y2] = c_map_best[x2,y2]
    # dropped += 1
    # Image.fromarray(np.uint8((c_map / np.max(c_map)).transpose() * 255), 'L').show()
    return new_loc, c_map_best


c1 = (7.5,7.5)
c2 = (5+2+4*24,5+2+3*24)
scale = 4

env = me.Environment()
env.load("minexml_test.xml")
mine = env.draw_basic_bitmap(scale)
unmod = np.copy(mine)
obs_mine = env.draw_obstacle_bitmap(scale)

# Empty the obstacle array so the algorithm works only
#     with discovered obstacles
# env.obstacles = []

# sl = (22,357)
# el = ((53 * scale + 2), 357 + ((48 * scale + 2)))

# path = []
# for i in range(sl[1], el[1]):
#     path.append((sl[0], i))
# for i in range(sl[0], el[0]):
#     path.append((i, el[1]))
# for i in range(el[1], sl[1], -1):
#     path.append((el[0], i))
# for i in range(el[0], sl[0], -1):
#     path.append((i, sl[1]))

# pathmap = np.dstack((mine, mine, mine))
# for p in path:
#     pathmap[floor(p[0] / scale), floor(p[1] / scale), 0] = 0
#     pathmap[floor(p[0] / scale), floor(p[1] / scale), 1] = 0
#     pathmap[floor(p[0] / scale), floor(p[1] / scale), 2] = 1

# pathmap *= 255
# img_path = Image.fromarray(np.uint8(pathmap.swapaxes(0,1)))
# img_path.show()

# og = np.copy(mine)

# for p in path:
#     mine = detect_obstacle((floor(p[0] / scale), floor(p[1] / scale)), 50, mine, obs_mine)

# obs_locs = []
# for x in range(len(mine)):
#     for y in range(len(mine[x])):
#         if og[x,y] - mine[x,y] == 1:
#             obs_locs.append((x,y))
#         # mine[x,y] = new_env[x,y]

# obs_dif = np.dstack((mine, obs_mine * mine, obs_mine * mine))
# for p in path:
#     obs_dif[floor(p[0] / scale), floor(p[1] / scale), 0] = 0
#     obs_dif[floor(p[0] / scale), floor(p[1] / scale), 1] = 0
#     obs_dif[floor(p[0] / scale), floor(p[1] / scale), 2] = 1

# obs_dif *= 255
# img_dif = Image.fromarray(np.uint8(obs_dif.swapaxes(0,1)))
# img_dif.show()

# for o in obs_locs:
#     tun_width = env.tunnels[0].width / 2
#     x_map = (o[0] + 0.5) * scale
#     y_map = (o[1] + 0.5) * scale
#     new_obs = me.Obstacle([
#         me.Point(x_map - tun_width, y_map - tun_width),
#         me.Point(x_map + tun_width, y_map + tun_width),
#         me.Point(x_map - tun_width, y_map + tun_width),
#         me.Point(x_map + tun_width, y_map - tun_width)
#     ])

#     env.add_obstacle(new_obs)

# pathmap = np.dstack((mine, mine, mine))

# pathmap *= 255
# img_path = Image.fromarray(np.uint8(pathmap.swapaxes(0,1)))
# img_path.show()

# Image.fromarray(np.uint8(env.draw_obstacle_bitmap(scale).transpose() * 255), 'L').show()

env.obstacles = []

# for x in range(len(mine)):
#     for y in range(len(mine[0])):
#         mine[x,y] = unmod[x,y]

# c_img = np.zeros(np.shape(mine))
# draw_circle([floor(118 / scale), floor(405 / scale)], 10, c_img)
# Image.fromarray(np.uint8(c_img.transpose() * 255), 'L').show()

# detect_obstacle([floor(118 / scale), floor(405 / scale)], 10, mine, obs_mine)
# Image.fromarray(np.uint8(mine.transpose() * 255), 'L').show()

# for i in range(40):
#     detect_obstacle([floor(94 / scale), floor(355 / scale) + i], 10, mine, obs_mine)
    
# Image.fromarray(np.uint8(mine.transpose() * 255), 'L').show()

# Define environment and channel parameters
noise = -120 # dBm
tx_pow = 3 # dBm
bw = 20 * 10**6 # Hz

# robot_loc = (22, 357)
rx_loc = Node((22,357,1), 0, float("inf"))
robot_loc = [floor(rx_loc.x / scale), floor(rx_loc.y / scale)]
height = 1
num_nodes = 4


# path, length = pathfind(mine, robot_loc, [robot_loc[0] + 12, robot_loc[1] + 6], [-1, -1])
# path = pathfind_bfs(mine, robot_loc, [robot_loc[0] + 12, robot_loc[1] + 6])

# path_img = np.repeat(mine[:,:,np.newaxis], 3, axis=2)
# for l in range(len(path)):
#     path_img[path[l][0], path[l][1], 0] = 0
#     path_img[path[l][0], path[l][1], 1] = 0
#     path_img[path[l][0], path[l][1], 2] = 1
# path_img *= 255
# path_img = Image.fromarray(np.uint8(path_img.swapaxes(0,1)))
# path_img.show()

dropped_nodes = [rx_loc]

cmap = np.zeros(np.shape(mine))

timer1 = time.time()

for i in range(num_nodes):
    timer2 = timer3 = time.time()
    
    while True:
        sum_before_path = np.sum(mine)
        new_node, cmap = predict_node(mine, cmap, scale, dropped_nodes, 3, noise, 350*10**6, env)
        print("\t\tTook {}s to predict node location".format(time.time() - timer3))
        path, length = pathfind_bfs(mine, robot_loc, [floor(new_node.x / scale), floor(new_node.y / scale)])
        # path_img = np.repeat(mine[:,:,np.newaxis], 3, axis=2)
        # for l in range(len(path)):
        #     path_img[path[l][0], path[l][1], 0] = 0
        #     path_img[path[l][0], path[l][1], 1] = 0
        #     path_img[path[l][0], path[l][1], 2] = 1
        # path_img *= 255
        # path_img = Image.fromarray(np.uint8(path_img.swapaxes(0,1)))
        # path_img.show()

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
                        new_node, cmap = predict_node(mine, cmap, scale, dropped_nodes, 3, noise, 350*10**6, env)
                        path, length = pathfind_bfs(mine, robot_loc, [floor(new_node.x / scale), floor(new_node.y / scale)])

                        # path_img = np.repeat(mine[:,:,np.newaxis], 3, axis=2)
                        # for l in range(len(path)):
                        #     path_img[path[l][0], path[l][1], 0] = 0
                        #     path_img[path[l][0], path[l][1], 1] = 0
                        #     path_img[path[l][0], path[l][1], 2] = 1
                        # path_img *= 255
                        # path_img = Image.fromarray(np.uint8(path_img.swapaxes(0,1)))
                        # path_img.show()

                        l = -1
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

print("Took {}s to drop all nodes".format(time.time() - timer1))

# node_locs, cmap = predict_node(mine, scale, rx_loc, num_nodes, 3, noise, 350*10**6, env)

cmap = cmap / np.max(cmap)
cmap_inv = np.multiply(1 - cmap, mine)

image = np.dstack((cmap_inv, cmap, np.zeros(np.shape(cmap))))

# obs_map = np.dstack((unmod, mine * unmod, unmod))
# obs_map *= 255
# img_obs = Image.fromarray(np.uint8(obs_map.swapaxes(0,1)))
# img_obs.show()

print("Node Locations: ")
for n in dropped_nodes:
    print("\t({x},{y})".format(x=n.x, y=n.y))
    image[floor(n.x / scale),floor(n.y / scale),0] = 0
    image[floor(n.x / scale),floor(n.y / scale),1] = 0
    image[floor(n.x / scale),floor(n.y / scale),2] = 1


image *= 255

img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
img_gen.show()
img_gen.save("GreedyFull_FinerResolution2.png")