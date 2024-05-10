from PIL import Image
from math import pi, sqrt, log2, log10, floor
import numpy as np
import gen_map
from matplotlib import pyplot as plt
import loadmap as me
import rt_slope
import raytrace as rt

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

def connected_corner_cvg(map, scale, loc:Node, num_nodes, tx_pow, noise, onehop_thpt, xml_env:me.Environment):
    dropped = 0
    d_loc = [loc]

    c_map = np.zeros(np.shape(map))
    for x in range(0,len(c_map)):
        for y in range(0,len(c_map[x])):
            if (mine[x][y] == 1):
                loss = rt.raytrace_loss(((x+0.5)*scale, (y+0.5)*scale, height), (loc.x, loc.y, loc.h), 2.4e9, xml_env)
                if (loss <= 0):
                    thpt = 0
                else:
                    loss = 10*log10(loss)
                # loss = calc_loss((loc.x, loc.y, loc.h), (loc.x, loc.y, loc.h), (slope, intercept), mine, scale)
                    rx_pow = tx_pow + loss
                    rx_snr = rx_pow - noise
                    thpt = bw * log2(1 + pow(10,(rx_snr / 10)))

                c_map[x][y] = thpt

    Image.fromarray(np.uint8((c_map / np.max(c_map)).transpose() * 255), 'L').show()

    while dropped < num_nodes:
        print("Node {}".format(dropped + 1))
        area = np.sum(c_map)
        c_map_best = np.zeros(np.shape(c_map))
        for x in range(0,len(map)):
            for y in range(0,len(map[x])):
                if map[x,y] == 1:
                    best_sig = 0
                    best_node = None
                    for n in d_loc:
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
        d_loc.append(new_loc)
        for x2 in range(len(c_map_new)):
            for y2 in range(len(c_map_new[0])):
                c_map[x2,y2] = c_map_best[x2,y2]
        dropped += 1
        Image.fromarray(np.uint8((c_map / np.max(c_map)).transpose() * 255), 'L').show()
    return d_loc, c_map


c1 = (7.5,7.5)
c2 = (5+2+4*24,5+2+3*24)
scale = 4

env = me.Environment()
env.load("minexml_test.xml")
mine = env.draw_1DTunnel_bitmap(scale)

Image.fromarray(np.uint8(mine.transpose() * 255), 'L').show()

# Define environment and channel parameters
noise = -120 # dBm
tx_pow = 3 # dBm
bw = 20 * 10**6 # Hz

# robot_loc = (22, 357)
rx_loc = Node((22,357,1), 0, float("inf"))
height = 1
num_nodes = 4

node_locs, cmap = connected_corner_cvg(mine, scale, rx_loc, num_nodes, 3, noise, 350*10**6, env)

cmap = cmap / np.max(cmap)
cmap_inv = np.multiply(1 - cmap, mine)

image = np.dstack((cmap_inv, cmap, np.zeros(np.shape(cmap))))

print("Node Locations: ")
for n in node_locs:
    print("\t({x},{y})".format(x=n.x, y=n.y))
    image[floor(n.x / scale),floor(n.y / scale),0] = 0
    image[floor(n.x / scale),floor(n.y / scale),1] = 0
    image[floor(n.x / scale),floor(n.y / scale),2] = 1


image *= 255

img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
img_gen.show()
img_gen.save("C:\\Users\\Worman 212\\Documents\\PDuane_NIOSH\\NodeDeploymentFullAlgo\\GreedyMineThroughput\\GreedyFull_FinerResolution2.png")