from __future__ import annotations
import xml.etree.ElementTree as et
from math import atan2, cos, sin, pi, sqrt, pow, log10
from matplotlib import pyplot as plot
import numpy as np
from PIL import Image
# import raytrace as rtloss
from scipy import stats

class Point:
    x:float
    y:float

    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def copy(self) -> Point:
        return Point(self.x, self.y)
    
    def distance(self, another:Point):
        return sqrt((another.x - self.x)**2 + (another.y - self.y)**2)

class Line:
    p1:Point
    p2:Point
    slope:float
    intercept:float

    def __init__(self, p1:Point, p2:Point):
        self.p1 = p1.copy()
        self.p2 = p2.copy()

        if (p2.x - p1.x) == 0:
            self.slope = float('inf')
            self.intercept = float('nan')
        else:
            self.slope = (p2.y-p1.y) / (p2.x-p1.x)
            self.intercept = p1.y - self.slope * p1.x
    
    # Find the intersection of two lines
    # Input: line to find intersection
    # Output: Coordinates of intersection if lines intersect or
    #         A Point of (inf, inf) if the lines are parallel
    #         A Point of (nan, nan) of the lines do not otherwise intersect
    def find_intersection(self, another:Line):
        if self.slope == another.slope:
            return Point(float('inf'), float('inf'))

        if np.isinf(self.slope):
            x = self.p1.x
            y = another.get_y(x)
        elif np.isinf(another.slope):
            x = another.p1.x
            y = self.get_y(x)
        else:
            x = (another.intercept - self.intercept) / (self.slope - another.slope)
            y = self.get_y(x)
            if np.isnan(y):
                x = float('nan')
        
        return Point(x, y)

    def intersects(self, another:Line):
        itsc = self.find_intersection(another)
        if self.p1.x > self.p2.x:
            sp1 = self.p2
            sp2 = self.p1
        else:
            sp1 = self.p1
            sp2 = self.p2
        if another.p1.x > another.p2.x:
            op1 = another.p2
            op2 = another.p1
        else:
            op1 = another.p1
            op2 = another.p2

        if itsc.x < sp1.x or itsc.x > sp2.x or itsc.x < op1.x or itsc.x > op2.x:
            return False
        
        if self.p1.y > self.p2.y:
            sp1 = self.p2
            sp2 = self.p1
        else:
            sp1 = self.p1
            sp2 = self.p2
        if another.p1.y > another.p2.y:
            op1 = another.p2
            op2 = another.p1
        else:
            op1 = another.p1
            op2 = another.p2
        if itsc.y < sp1.y or itsc.y > sp2.y or itsc.y < op1.y or itsc.y > op2.y:
            return False
        return not (np.isinf(itsc.x) or np.isinf(itsc.y)) and not (np.isnan(itsc.x) or np.isnan(itsc.y))
    
    def get_y(self, x:float, ignore_ends = False):
        if not ignore_ends:
            if x < min(self.p1.x, self.p2.x) or x > max(self.p1.x, self.p2.x):
                return float('nan')

        return self.slope * x + self.intercept
    
    def get_x(self, y:float, ignore_ends = False):
        if not ignore_ends:
            if y < min(self.p1.y, self.p2.y) or y > max(self.p1.y, self.p2.y):
                return float('nan')

        return (y - self.intercept) / self.slope

    def length(self):
        return self.p1.distance(self.p2)

class Tunnel:
    centerline:Line
    start:Point
    end:Point
    width:float
    height:float
    eps_wall:float
    sig_wall:float
    eps_ceil:float
    sig_ceil:float
    loss_est:tuple #Linear loss in decibels per unit distance

    angle:float
    corners:tuple
    connected:tuple

    def __init__(self, start:Point, end:Point, width:float, height:float, eps_wall:float, sig_wall:float, eps_ceil:float, sig_ceil:float):
        self.start = start
        self.end = end
        self.centerline = Line(start, end)
        self.width = width
        self.height = height
        self.eps_wall = eps_wall
        self.sig_wall = sig_wall
        self.eps_ceil = eps_ceil
        self.sig_ceil = sig_ceil

        self.connected = []

        self.angle = angle = atan2(end.y-start.y, end.x-start.x)
        c1x = self.start.x + self.width / 2 * cos(angle + pi/2)
        c1y = self.start.y + self.width / 2 * sin(angle + pi/2)
        c1 = Point(c1x, c1y)
        c2x = self.start.x + self.width / 2 * cos(angle - pi/2)
        c2y = self.start.y + self.width / 2 * sin(angle - pi/2)
        c2 = Point(c2x, c2y)
        c3x = self.end.x + self.width / 2 * cos(angle - pi/2)
        c3y = self.end.y + self.width / 2 * sin(angle - pi/2)
        c3 = Point(c3x, c3y)
        c4x = self.end.x + self.width / 2 * cos(angle + pi/2)
        c4y = self.end.y + self.width / 2 * sin(angle + pi/2)
        c4 = Point(c4x, c4y)

        self.corners = (c1, c2, c3, c4)

        loss_est = (float('nan'), float('nan'))

    # def calc_slope(self, freq, height) -> tuple:
    #     ys = np.linspace(1, 600, 1000)
    #     values = []
    #     h = -(self.height / 2) + height
    #     for y in ys:
    #         loss = rtloss.raytrace_loss_raw((0, y, h), (0, 0, h), freq, self.width, self.height, (self.eps_ceil, self.sig_ceil), (self.eps_wall, self.sig_wall))
    #         loss = 10*log10(loss)
    #         values.append(loss)

    #     slope, intercept, r_value, p_value, std_err = stats.linregress(ys,values)
    #     self.loss_est = (slope, intercept)
    #     return (slope, intercept)
        # rud = (slope * ys) + intercept
        # plt.plot(ys, values)
        # plt.plot(ys, rud)
        # plt.show()

    def get_edges(self):
        return [
                Line(self.corners[0], self.corners[1]),
                Line(self.corners[1], self.corners[2]),
                Line(self.corners[2], self.corners[3]),
                Line(self.corners[3], self.corners[0])
            ]
    
    def in_tunnel(self, p:Point):
        c = self.corners
        area_tri = 0
        for i in range(0,4):
            j = (i + 1) % 4
            area_tri += abs((c[i].x * (c[j].y - p.y) + c[j].x * (p.y - c[i].y) + p.x * (c[i].y - c[j].y)) / 2)

        dx = sqrt(pow(c[1].x - c[0].x, 2) + pow(c[1].y - c[0].y, 2))
        dy = sqrt(pow(c[3].x - c[0].x, 2) + pow(c[3].y - c[0].y, 2))

        area_rect = dx * dy

        return area_tri <= area_rect + 0.000001
    
    def find_connected(self, tunnels:list, save = True):
        conn = []
        for t in tunnels:
            if t == self:
                continue

            if self.centerline.intersects(t.centerline):
                conn.append(t)
                continue
            
            num_itsc = 0
            skip = False
            for e in self.get_edges():
                for te in t.get_edges():
                    if e.intersects(te):
                        num_itsc += 1
                        if num_itsc >= 2:
                            conn.append(t)
                            skip = True
                            break
                if skip:
                    break
        if save:
            self.connected = conn
        
        return conn

    def is_connected(self, another:Tunnel):
        if another == self:
            return True
        if another in self.connected:
            return True
        
        if self.centerline.intersects(another.centerline):
            return True
        
        num_itsc = 0
        for e in self.get_edges():
            for te in another.get_edges():
                if e.intersects(te):
                    num_itsc += 1
                    if num_itsc >= 2:
                        return True
        return False
    
    def get_relative_coordinate(self, p:Point) -> Point:
        if not self.in_tunnel(p):
            return Point(float('nan'), float('nan'))

        a = atan2(p.y - self.start.y, p.x - self.start.x)
        l = sqrt(pow(p.x - self.start.x, 2) + pow(p.y - self.start.y, 2))

        rel_angle = a - self.angle
        rel_x = l * cos(rel_angle + pi/2)
        rel_y = l * sin(rel_angle + pi/2)

        return Point(rel_x, rel_y)

    def draw_tunnel(self, show = False):
        c = self.corners
        xs = [c[0].x, c[1].x, c[2].x, c[3].x, c[0].x]
        ys = [c[0].y, c[1].y, c[2].y, c[3].y, c[0].y]
        p = plot.plot(xs, ys)
        if (show):
            plot.show()
        return p

class Obstacle:
    corners:tuple
    edges:tuple
    def __init__(self, corners:tuple):
        s_x = 0
        s_y = 0
        for c in corners:
            s_x += c.x
            s_y += c.y
        
        midpoint = Point(s_x / 4, s_y / 4)

        angles = []
        for i in range(len(corners)):
            a = atan2(midpoint.y - corners[i].y, midpoint.x - corners[i].x)
            a = (a + 2*pi) % (2*pi)
            angles.append(a)

        corners_sorted = []
        num_its = len(corners)
        for i in range(num_its):
            min_idx = 0
            min = float('inf')

            for j in range(len(angles)):
                if angles[j] < min:
                    min_idx = j
                    min = angles[j]
            corners_sorted.append(corners[min_idx])
            del corners[min_idx]
            del angles[min_idx]

        self.corners = corners_sorted
        self.edges = (
            Line(self.corners[0], self.corners[1]),
            Line(self.corners[1], self.corners[2]),
            Line(self.corners[2], self.corners[3]),
            Line(self.corners[3], self.corners[0])
        )

    def inside(self, p:Point):
        c = self.corners
        area_tri = 0
        for i in range(0,4):
            j = (i + 1) % 4
            area_tri += abs((c[i].x * (c[j].y - p.y) + c[j].x * (p.y - c[i].y) + p.x * (c[i].y - c[j].y)) / 2)

        dx = sqrt(pow(c[1].x - c[0].x, 2) + pow(c[1].y - c[0].y, 2))
        dy = sqrt(pow(c[3].x - c[0].x, 2) + pow(c[3].y - c[0].y, 2))

        area_rect = dx * dy

        return area_tri <= area_rect + 0.000001

    def line_enters(self, l:Line):
        itscs = 0
        for e in self.edges:
            if l.intersects(e):
                itscs += 1
        if itscs == 1:
            return True
        
        return False

    def line_through(self, l:Line):
        itscs = 0
        for e in self.edges:
            if l.intersects(e):
                itscs += 1
        if itscs == 2:
            return True
        
        return False

    def collides(self, l:Line):
        return self.line_enters(l) or self.line_through(l)

class Environment:
    tunnels:tuple
    obstacles:tuple
    lower_bound:Point
    upper_bound:Point
    
    def load(self, file, load_obstacles = True):
        env = []
        obs = []
        tree = et.parse(file)
        root = tree.getroot()
        if (root.tag != "mine"):
            raise Exception("Bad root tag: {}; Should be 'mine'".format(root.tag))
        
        lb_str = root.attrib['lower_bound']
        op = lb_str.find('(')
        cp = lb_str.find(')')
        lb_str = lb_str[op+1:cp]
        lbx, lby = lb_str.split(',')
        self.lower_bound = Point(float(lbx), float(lby))

        ub_str = root.attrib['upper_bound']
        op = ub_str.find('(')
        cp = ub_str.find(')')
        ub_str = ub_str[op+1:cp]
        ubx, uby = ub_str.split(',')
        self.upper_bound = Point(float(ubx), float(uby))

        for child in root:
            if child.tag == "tunnel":
                p1 = child.attrib['p1']
                op = p1.find('(')
                cp = p1.find(')')
                p1 = p1[op+1:cp]
                p1x, p1y = p1.split(',')
                p1x = p1x.strip()
                p1y = p1y.strip()
                p1_p = Point(float(p1x), float(p1y))
                
                p2 = child.attrib['p2']
                op = p2.find('(')
                cp = p2.find(')')
                p2 = p2[op+1:cp]
                p2x, p2y = p2.split(',')
                p2x = p2x.strip()
                p2y = p2y.strip()
                p2_p = Point(float(p2x), float(p2y))

                w = float(child.attrib["width"])
                h = float(child.attrib["height"])

                ew = float(child.attrib["eps_wall"])
                sw = float(child.attrib["sig_wall"])
                ec = float(child.attrib["eps_ceil"])
                sc = float(child.attrib["sig_ceil"])

                env.append(Tunnel(p1_p, p2_p, w, h, ew, sw, ec, sc))
            elif load_obstacles and child.tag == "obstacle":
                corners = []
                for sc in child:
                    if sc.tag == "corner":
                        p = sc.attrib['p']
                        op = p.find('(')
                        cp = p.find(')')
                        p = p[op+1:cp]
                        px, py = p.split(',')
                        px = px.strip()
                        py = py.strip()
                        corners.append(Point(float(px), float(py)))
                if not len(corners) == 4:
                    print("Bad Obstacle, {} points".format(len(corners)))
                else:
                    obs.append(Obstacle(corners))
        
        self.tunnels = env
        self.obstacles = obs
        return self
    
    def add_obstacle(self, obs:Obstacle):
        self.obstacles.append(obs)

    def connect_tunnels(self):
        for t in self.tunnels:
            t.find_connected(self.tunnels, save=True)

    def t_los(self, line:Line, t:Tunnel, visited:list):
        p = line.p2
        visited.append(t)

        if t.in_tunnel(p):
            return True, visited

        connected = []

        if t.in_tunnel(p):
            return True, visited
        edges_in_connected = False
        for e in t.get_edges():
            if line.intersects(e):
                itsc = line.find_intersection(e)
                for tc in t.connected:
                    if tc.in_tunnel(itsc) and not tc in visited:
                        edges_in_connected = True
                        connected.append(tc)
        if not edges_in_connected:
            return False, visited
        
        for t2 in connected:
            lt, vis = self.t_los(line, t2, visited)
            visited = vis
            if lt:
                return True, visited
        
        return False, visited

    # For now, making naive assumption that two points only have line of sight
    #     if they are in the same tunnel
    def has_los(self, p1:Point, p2:Point):
        los_line = Line(p1, p2)

        for o in self.obstacles:
            if o.collides(los_line):
                return False
        
        t1s = []
        for t in self.tunnels:
            if (t.in_tunnel(p1) and t.in_tunnel(p2)):
                return True
            elif t.in_tunnel(p1):
                t1s.append(t)
                t1:Tunnel = t
            elif t.in_tunnel(p2):
                t2:Tunnel = t

        # t = t1
        visited = []
        for t in t1s:
            los, visited = self.t_los(los_line, t, visited)
            if los:
                # You can get the starting tunnel for the LOS condition with this:
                # t1 = t
                return True
        
        return False
    
    # returns:
    #     d_tot --> the total Manhattan distance between the points
    #     s_dist1 --> The distance between the first point and the corner point
    #     s_dist2 --> The distance between the corner point and the second point
    #     t_fin1 --> The tunnel containing the first point
    #     t_fin2 --> The tunnel containing the second point
    def get_distance(self, p1:Point, p2:Point):
        if (self.has_los(p1,p2)):
            tun = None
            for t in self.tunnels:
                if t.in_tunnel(p1) and t.in_tunnel(p2):
                    tun = t
                    break
            return p1.distance(p2), float('nan'), float('nan'), tun, tun

        t1s = []
        t2s = []
        for t in self.tunnels:
            if t.in_tunnel(p1):
                t1s.append(t)
            if t.in_tunnel(p2):
                t2s.append(t)
        
        c_pairs = []
        for t1 in t1s:
            for t2 in t2s:
                if t1.is_connected(t2):
                    c_pairs.append((t1, t2))
        
        if len(c_pairs) == 0:
            return float('inf'), float('nan'), float('nan'), None, None

        s_dist = [float('inf'), float('inf'), float('inf')]
        t_fin = [None, None]

        for c in c_pairs:
            t1 = c[0]
            t2 = c[1]

            # If the center lines intersect
            if (t1.centerline.intersects(t2.centerline)):
                itsc = t1.centerline.find_intersection(t2.centerline)

                p1_rel = t1.get_relative_coordinate(p1)
                t1_itc = t1.get_relative_coordinate(itsc)
                p2_rel = t2.get_relative_coordinate(p2)
                t2_itc = t2.get_relative_coordinate(itsc)

                d_p1 = abs(p1_rel.distance(t1_itc))
                d_p2 = abs(p2_rel.distance(t2_itc))

                if not self.has_los(p1, itsc):
                    d_p1 = float('inf')
                if not self.has_los(p2, itsc):
                    d_p2 = float('inf')

                if (d_p1 == float('inf') or d_p2 == float('inf')):
                    d_tot = float('inf')
                else:
                    d_tot = d_p1 + d_p2

                if (d_tot < s_dist[0]):
                    s_dist[0] = d_tot
                    s_dist[1] = d_p1
                    s_dist[2] = d_p2
                    t_fin[0] = t1
                    t_fin[1] = t2
            
                # If the center lines don't intersect;
                if (not t1.centerline.slope == 0):
                    norm_itcp_s = t2.start.y + t2.start.x * 1/t1.centerline.slope
                    norm_itcp_e = t2.end.y + t2.end.x * 1/t1.centerline.slope

                    norm_start = Line(Point(0, norm_itcp_s), Point(t2.start.x, t2.start.y))
                    norm_end = Line(Point(0, norm_itcp_e), Point(t2.end.x, t2.end.y))

                    norm_s_t1_cross = norm_start.find_intersection(t1.centerline)
                    norm_e_t1_cross = norm_end.find_intersection(t1.centerline)

                    norm_start = Line(norm_s_t1_cross, Point(t2.start.x, t2.start.y))
                    norm_end = Line(norm_e_t1_cross, Point(t2.end.x, t2.end.y))
                else:
                    norm_start = Line(Point(t2.start.x, t1.start.y), Point(t2.start.x, t2.start.y))
                    norm_end = Line(Point(t2.end.x, t1.start.y), Point(t2.end.x, t2.end.y))

                if norm_end.length() < norm_start.length():
                    norm = norm_end
                else:
                    norm = norm_start

                p_mid = Point((norm.p1.x + norm.p2.x) / 2, (norm.p1.y + norm.p2.y) / 2)

                d_p1 = p1.distance(p_mid)
                d_p2 = p2.distance(p_mid)

                if not self.has_los(p1, p_mid):
                    d_p1 = float('inf')
                if not self.has_los(p2, p_mid):
                    d_p2 = float('inf')

                if (d_p1 == float('inf') or d_p2 == float('inf')):
                    d_tot = float('inf')
                else:
                    d_tot = d_p1 + d_p2

                d_tot = d_p1 + d_p2
                if (d_tot < s_dist[0]):
                    s_dist[0] = d_tot
                    s_dist[1] = d_p1
                    s_dist[2] = d_p2
                    t_fin[0] = t1
                    t_fin[1] = t2

        return s_dist[0], s_dist[1], s_dist[2], t_fin[0], t_fin[1]

    def dist_in_tunnels(self, p1:Point, p2:Point):
        for t in self.tunnels:
            if (t.in_tunnel(p1) and t.in_tunnel(p2)):
                p1_rel = t.get_relative_coordinate(p1)
                p2_rel = t.get_relative_coordinate(p2)
                return Point(p2_rel.x - p1_rel.x, p2_rel.y - p1_rel.y)

        p1_tunnels = []
        p2_tunnels = []

        for t in self.tunnels:
            if t.in_tunnel(p1):
                p1_tunnels.append(t)
            if t.in_tunnel(p2):
                p2_tunnels.append(t)
        
        t_dists = []
        for t1 in p1_tunnels:
            for t2 in p2_tunnels:
                if t1.is_connected(t2):
                    d_tot, d_p1, d_p2 = self.get_distance(p1, p2, t1, t2)


    def draw_bitmap(self, resolution):
        nx = round((self.upper_bound.x - self.lower_bound.x) / resolution)
        ny = round((self.upper_bound.y - self.lower_bound.y) / resolution)

        bmp = np.zeros((nx,ny))

        for x in range(0, nx):
            for y in range(0, ny):
                lx = resolution * (x + 0.5) + self.lower_bound.x
                ly = resolution * (y + 0.5) + self.lower_bound.y

                # in_tun = False
                in_obs = False
                p = Point(lx, ly)
                for o in self.obstacles:
                    if o.inside(p):
                        in_obs = True
                        break
                if in_obs:
                    continue

                for tun in self.tunnels:
                    if tun.in_tunnel(Point(lx, ly)):
                        bmp[x,y] = 1
                        break
        
        return bmp
    
    def draw_basic_bitmap(self, resolution):
        nx = round((self.upper_bound.x - self.lower_bound.x) / resolution)
        ny = round((self.upper_bound.y - self.lower_bound.y) / resolution)

        bmp = np.zeros((nx,ny))

        for x in range(0, nx):
            for y in range(0, ny):
                lx = resolution * (x + 0.5) + self.lower_bound.x
                ly = resolution * (y + 0.5) + self.lower_bound.y

                # in_tun = False
                # in_obs = False
                # p = Point(lx, ly)
                # for o in self.obstacles:
                #     if o.inside(p):
                #         in_obs = True
                #         break
                
                for tun in self.tunnels:
                    if tun.in_tunnel(Point(lx, ly)):
                        relative = tun.get_relative_coordinate(Point(lx, ly))
                        if (abs(relative.x) < (resolution / 2) + 0.0001):
                            bmp[x,y] = 1
                            break
        
        return bmp
    
    def draw_obstacle_bitmap(self, resolution):
        nx = round((self.upper_bound.x - self.lower_bound.x) / resolution)
        ny = round((self.upper_bound.y - self.lower_bound.y) / resolution)

        bmp = np.zeros((nx,ny))

        for x in range(0, nx):
            for y in range(0, ny):
                lx = resolution * (x + 0.5) + self.lower_bound.x
                ly = resolution * (y + 0.5) + self.lower_bound.y

                # in_tun = False
                in_obs = False
                p = Point(lx, ly)
                for o in self.obstacles:
                    if o.inside(p):
                        in_obs = True
                        break
                if in_obs:
                    continue
                
                for tun in self.tunnels:
                    if tun.in_tunnel(Point(lx, ly)):
                        relative = tun.get_relative_coordinate(Point(lx, ly))
                        if (abs(relative.x) < (resolution / 2) + 0.0001):
                            bmp[x,y] = 1
                            break
        
        return bmp

# mine = Environment().load("minexml_test.xml")
# mine.connect_tunnels()

# p1 = Point(70,103)
# p2 = Point(140,103)
# print("Line of Sight for ({},{}) to ({},{}):".format(p1.x, p1.y, p2.x, p2.y))
# if mine.has_los(p1, p2):
#     print("True")

#     d = mine.get_distance(p1,p2)
#     print("\nDistance:  {}".format(d[0]))
#     print("Euclidean: {}".format(p1.distance(p2)))
# else:
#     print("False")
#     d = mine.get_distance(p1,p2)
#     if d[0] == float('inf'):
#         print("No Connection")
#     else:
#         print("\nDistance:  {}".format(d[0]))
#         print("Segments: {}, {}".format(d[1], d[2]))

# print("\n")

# bmp = mine.draw_bitmap(1) * 255
# i = Image.fromarray(bmp.transpose())
# i = i.convert("RGB")
# i.putpixel((p1.x, p1.y), (255,0,0))
# i.putpixel((p2.x, p2.y), (255,0,0))
# i.save("TestMap.png")