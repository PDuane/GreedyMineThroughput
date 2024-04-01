import loadmap as me
import numpy as np
import raytrace as rtloss
from math import log10
from PIL import Image
import sys
from multiprocessing import Process, Value, Array, Lock
import ctypes

def calc_raytrace(tx:tuple, map:Array, omap:Array, map_size:tuple, mine:me.Environment, active, lock:Lock):
    env_map = np.frombuffer(map.get_obj()).reshape(map_size[0], map_size[1])
    otpt_map = np.frombuffer(omap.get_obj()).reshape(map_size[0], map_size[1])
    for x in range(len(env_map)):
        for y in range(len(env_map[x1])):
            if env_map[x][y] == 1:
                x1_val = int((tx[0]+0.5) / scale)
                y1_val = int((tx[1]+0.5) / scale)
                x2_val = int((x+0.5) / scale)
                y2_val = int((y+0.5) / scale)
                snr = rtloss.raytrace_loss((x1_val, y1_val,height), (x2_val,y2_val,height), 2.4e9, mine)
                sum += snr
    lock.acquire()
    otpt_map[tx[0], tx[1]] = snr
    lock.release()
    active.value -= 1


if __name__ =="__main__":
    num_processes = 8
    active_processes = Value("i", 0)
    print("Loading environment")
    mine = me.Environment()
    mine.load("minexml_test.xml")

    scale = .25 # meters per pixel
    height=1

    print("Generating bitmap")
    phys_map = mine.draw_bitmap(scale)
    cmap = np.zeros(phys_map.shape)

    mp_phys_map = Array(ctypes.c_double, len(phys_map) * len(phys_map[0]))
    tmp_p = np.frombuffer(mp_phys_map.get_obj()).reshape(len(phys_map), len(phys_map[0]))
    mp_cmap = Array(ctypes.c_double, len(cmap) * len(cmap[0]))
    tmp_c = np.frombuffer(mp_cmap.get_obj()).reshape(len(cmap), len(cmap[0]))

    for x in range(len(phys_map)):
        for y in range(len(phys_map[x])):
            tmp_p[x][y] = phys_map[x][y]
            tmp_c[x][y] = cmap[x][y]
    
    phys_map = tmp_p
    cmap = tmp_c

    mutex = Lock()

    print("Calculating Coverage:")
    for x1 in range(len(phys_map)):
        for y1 in range(len(phys_map[x1])):
            coord_str = "\r\t({x},{y})          ".format(x = x1, y = y1)
            sys.stdout.write(coord_str)
            if phys_map[x1][y1] == 1:
                sum = 0
                # while (active_processes.value > (num_processes - 1)):
                #     pass
                # p = Process(target=calc_raytrace, args=((x1, y1, 2), mp_phys_map, mp_cmap, (len(phys_map), len(phys_map[0])), mine, active_processes, mutex))
                # active_processes.value += 1
                # coord_str = "\r\t({x},{y})          ".format(x = x1, y = y1)
                # sys.stdout.write(coord_str)
                for x2 in range(len(phys_map)):
                    for y2 in range(len(phys_map[x1])):
                        if phys_map[x2][y2] == 1:
                            x1_val = int((x1+0.5) / scale)
                            y1_val = int((y1+0.5) / scale)
                            x2_val = int((x2+0.5) / scale)
                            y2_val = int((y2+0.5) / scale)
                            snr = rtloss.raytrace_loss((x1_val, y1_val,height), (x2_val,y2_val,height), 2.4e9, mine)
                            sum += snr
                cmap[x1, y1] = sum
            # else:
            #     coord_str = "\r\t({x},{y})          ".format(x = x1, y = y1)
            #     sys.stdout.write(coord_str)

    cmap_max = np.nanmax(cmap)
    cmap /= cmap_max
    cvg_map = np.repeat(phys_map[:,:,np.newaxis], 3, axis=2)
    for x in range(len(cvg_map)):
        for y in range(len(cvg_map[x])):
            if phys_map[x, y] == 1:
                cvg_map[x,y,0] = 1-cmap[x,y]
                cvg_map[x,y,1] = cmap[x,y]
                cvg_map[x,y,2] = 0

    cvg_map *= 255

    img = Image.fromarray(np.uint8(cvg_map.swapaxes(0,1)))
    img.show()
