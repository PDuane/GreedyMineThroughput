import loadmap as me
from PIL import Image
import numpy as np

scale = 4

env = me.Environment()
env.load("SimRig_Obstacle.xml")
mine = env.draw_obstacle_bitmap(scale)

img = Image.fromarray(mine.transpose() * 255)
img.show()