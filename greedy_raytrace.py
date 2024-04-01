from PIL import Image
from math import pi, sqrt, log10
import numpy as np
import gen_map
from matplotlib import pyplot as plt
import loadmap as me
from raytrace import raytrace_loss as rtloss

# args:
#     n1_height - height of the first (transmitting) node
#     n2_height - height of the second (receiving) node
#     freq      - Center frequency
#     env       - The Environment object encoding the mine
def rt_adapter(n1, n2, args):
    n1_height = args[0]
    n2_height = args[1]
    freq = args[2]
    env = args[3]
    loss = rtloss((n1[0], n1[0], n1_height), (n2[0], n2[0], n2_height), freq, env)
    # loss = 1 - loss
    if loss < 10e-13:
        loss = 10e-13
    return 10 * log10(loss)

def greedy_connected_cvg(map, c_map, loc, num_nodes, freq, env):
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
                        if rtloss((n[0],n[1],1), (x,y,1), freq, env) > -100:
                        # if conn_alg((x, y), n, [*args[:],map]) > -120:
                            con = True
                            break
                    if not con:
                        continue
                    c_map_new = np.zeros((len(c_map),len(c_map[0])))
                    for x2 in range(len(c_map)):
                        for y2 in range(len(c_map[0])):
                            c_map_new[x2,y2] = c_map[x2,y2]
                    # c_map_new = c_map.copy()
                    for x2 in range(0,len(c_map_new)):
                        for y2 in range(0,len(c_map_new[x2])):
                            if rtloss((x,y,1), (x2,y2,1), freq, env) > -100:
                            # if conn_alg((x2, y2), (x,y), [*args[:],map]):
                                c_map_new[x2][y2] = 1
                    
                    # Image.fromarray(np.uint8(c_map_new.transpose() * 255), 'L').show()
                    # diff_map = np.subtract(c_map_new, c_map)
                    # Image.fromarray(np.uint8(diff_map.transpose() * 255), 'L').show()
                    new_coverage = np.sum(c_map_new)
                    if new_coverage > area:
                        area = new_coverage
                        new_loc = (x,y)
                        # for x2 in range(len(c_map_new)):
                        #     for y2 in range(len(c_map_new[0])):
                        #         c_map[x2,y2] = c_map_new[x2,y2]
                        # Image.fromarray(np.uint8(c_map_new.transpose() * 255), 'L').show()
        d_loc.append(new_loc)
        for x2 in range(0,len(c_map_new)):
            for y2 in range(0,len(c_map_new[x2])):
                if rtloss((new_loc[0],new_loc[1],1), (x2,y2,1), freq, env) > -100:
                # if conn_alg((x2, y2), new_loc, [*args[:],map]) > -120:
                    c_map[x2][y2] = 1
        dropped += 1
        # c_map = c_map_new
        # Image.fromarray(np.uint8(c_map.transpose() * 255), 'L').show()
    return d_loc, c_map

mine = me.Environment()
mine.load("minexml_test.xml")

scale = 1 # meters per pixel

phys_map = mine.draw_bitmap(scale)
cmap = np.zeros(phys_map.shape)

robot_loc = (7, 7)
node_height = 1
freq = 2.4e9

num_nodes = 5

# comm_dist = 100

# alg_args = [10, -80, 2.4e9]

str_map = np.empty(phys_map.shape)
str_map[:] = np.nan

for x_raw in range(0,len(cmap)):
    x = (x_raw + 1/2) * scale
    for y_raw in range(0,len(cmap[x_raw])):
        y = (y_raw + 1/2) * scale
        if phys_map[x_raw][y_raw] == 1:
            stn = rtloss((robot_loc[0], robot_loc[1], node_height), (x, y, node_height), freq, mine)
            if (stn < 10e-12):
                stn = 10e-12
            stn = 10 * log10(stn)
            str_map[x_raw][y_raw] = stn
            if (rtloss((robot_loc[0], robot_loc[1], node_height), (x, y, node_height), freq, mine)) > -120:
                cmap[x_raw][y_raw] = 1

nv = np.nanmin(str_map)
str_map = str_map - nv
xv = np.nanmax(str_map)
str_map = str_map / xv

str_img = np.repeat(phys_map[:,:,np.newaxis], 3, axis=2)
for x in range(len(str_map)):
    for y in range(len(str_map[x])):
        if not str_map[x][y] == float('nan'):
            str_img[x,y,0] = 1-str_map[x][y]
            str_img[x,y,1] = str_map[x][y]
            str_img[x,y,2] = 0

str_img *= 255

str_img_gen = Image.fromarray(np.uint8(str_img.swapaxes(0,1)))
str_img_gen.show()

# pv_img = np.repeat(phys_map[:,:,np.newaxis], 3, axis=2)
# pv_img[:,:,0] = np.subtract(pv_img[:,:,0], cmap)
# pv_img[:,:,2] = np.subtract(pv_img[:,:,2], cmap)

# pv_img *= 255

# pv_img_gen = Image.fromarray(np.uint8(pv_img.swapaxes(0,1)))
# pv_img_gen.show()

node_locs, new_map = greedy_connected_cvg(phys_map, cmap, robot_loc, num_nodes, 2.4e9, mine)

image = np.repeat(phys_map[:,:,np.newaxis], 3, axis=2)
image[:,:,0] = np.subtract(image[:,:,0], cmap)
image[:,:,2] = np.subtract(image[:,:,2], cmap)

for n in node_locs:
    image[round(n[0]),round(n[1]),0] = 1
    image[round(n[0]),round(n[1]),1] = 0
    image[round(n[0]),round(n[1]),2] = 0

image *= 255

img_gen = Image.fromarray(np.uint8(image.swapaxes(0,1)))
img_gen.show()
img_gen.save("rt_nodeplacement.png")
