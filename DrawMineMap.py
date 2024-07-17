import loadmap as me
from PIL import Image
import numpy as np

scale = 4

env = me.Environment()
env.load("ComplexRP.xml")
mine = env.draw_obstacle_bitmap(scale)

mine_img = np.repeat(mine[:,:,np.newaxis], 3, axis=2)
mine_img *= 255
img_gen = Image.fromarray(np.uint8(mine_img.swapaxes(0,1)))
# img = Image.fromarray(mine_img.transpose() * 255)
img_gen.show()
img_gen.save("ComplexRP_BaseMap_small.png")