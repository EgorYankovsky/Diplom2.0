import os

class Rib:
    
    def get_antinormal(self):
        rib_len = ((self.x1 - self.x0) ** 2 + (self.y1 - self.y0) ** 2 + (self.z1 - self.z0) ** 2) ** 0.5
        return ((self.x1 - self.x0) / rib_len, 
                (self.y1 - self.y0) / rib_len,
                (self.z1 - self.z0) / rib_len)

    def get_middle_point(self):
        return ((self.x1 + self.x0) * 0.5, 
                (self.y1 + self.y0) * 0.5,
                (self.z1 + self.z0) * 0.5)

    def get_values_to_draw(self):
        return (self.x0, self.y0, self.z0, self.dx, self.dy, self.dz)

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __init__(self, x0, y0, z0, dx, dy, dz):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.dx = dx
        self.dy = dy
        self.dz = dz


def is_unique(arr, elem) -> bool:
    for item in arr:
        if item == elem:
            return False
    return True

max_values = []

x = []
y = []
z = []
func = []
ribs_arr = []
it = 0

def read_points(path):
    f = open(path, "r")
    f.readline()
    for line in f.read().split("\n"):
        if line == '': break
        coord = line.split()
        ri = float(line.split()[0])
        zi = float(line.split()[1])
        if is_unique(r, ri):
            r.append(ri)
        if is_unique(z, zi):
            z.append(zi)
        f.close()

def read_data(path, r, z, func) -> (list, list, list): # type: ignore
    f = open(path, 'r')
    r_curr = []
    z_curr = []
    for line in f.read().split('\n'):
        if line == '': break
        info = line.split(" ")
        r_curr.append(float(info[0]))
        z_curr.append(float(info[1]))
        func.append(abs(float(info[2])))
    r = list(dict.fromkeys(r_curr).keys())
    z = list(dict.fromkeys(z_curr).keys())
    return (r, z, func)

def read_ribs(ribs_path):
    f = open(ribs_path, 'r')
    for line in f.read().split('\n'):
        if line == '': break
        ribs.append(rib(line.split(" ")))
    f.close()

def read_vectors(vectors_path):
    ribs_arr.clear()
    i = 0
    for text in open(vectors_path, 'r').read().split("\n"):
        if text == "":
           break
        text_vector = text.split(" ")
        if i == 0:
            max_values.append(float(text_vector[0]))
            max_values.append(float(text_vector[1]))
            max_values.append(float(text_vector[2]))
            max_values.append(float(text_vector[3]))
            max_values.append(float(text_vector[4]))
            max_values.append(float(text_vector[5]))
            i += 1
        else:
            ribs_arr.append(Rib(float(text_vector[0]), float(text_vector[1]),float(text_vector[2]),
                                float(text_vector[3]), float(text_vector[4]),float(text_vector[5])))
