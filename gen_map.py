import numpy as np

# img_size = (128, 128)
# hall_width = 4
# pillar_width = 20
# pillar_length = 20

# pixel_size = 1

# # Top, left, right, bottom
# boundary = (5, 5, 5, 5)

# # Add additional open areas with a 4-tuple of format
# #     (x1, y1, x2, y2)
# passageways = [(5,0,9,5)]

# # Add additional blocked areas with a 4-tuple of format
# #     (x1, y1, x2, y2)
# blockages = []

def gen_map(img_size, tunnel_width, pillar_width, pillar_length, pixel_size, boundary, passageways, blockages):

    image = np.zeros((128, 128))

    for y in range(boundary[0], img_size[1] - boundary[3]):
        for x in range(boundary[1], img_size[0] - boundary[2]):
            x_clear = (x - boundary[0]) % (pillar_width + tunnel_width) < tunnel_width
            y_clear = (y - boundary[1]) % (pillar_length + tunnel_width) < tunnel_width
            if x_clear or y_clear:
                image[x // pixel_size, y // pixel_size] = 1

    for p in passageways:
        for x in range(p[0], p[2]):
            for y in range (p[1], p[3]):
                image[x // pixel_size, y // pixel_size] = 1

    for b in blockages:
        for x in range(b[0], b[2]):
            for y in range (b[1], b[3]):
                image[x // pixel_size, y // pixel_size] = 0

    return image

# img_gen = Image.fromarray(np.uint8(image), 'L')

# img_gen.show()