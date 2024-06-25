import loadmap as me
from PIL import Image
import numpy as np

scale = 6

env = me.Environment()
env.load("Mine1.xml")
mine = env.draw_basic_bitmap(scale)

img = Image.fromarray(mine.transpose() * 255)
img.show()