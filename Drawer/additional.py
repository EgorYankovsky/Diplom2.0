import os

class rib:

    def __init__(self, begin, end):
        self.begin = begin
        self.end = end

    def __init__(self, x0, y0, z0, x1, y1, z1):
        self.x0 = float(x0)
        self.y0 = float(y0)
        self.z0 = float(z0)

        self.x1 = float(x1)
        self.y1 = float(y1)
        self.z1 = float(z1)

    def __init__(self, arr):
        if len(arr) > 6:
            self.x0 = 0.0
            self.y0 = 0.0
            self.z0 = 0.0

            self.x1 = 0.0
            self.y1 = 0.0
            self.z1 = 0.0
        elif len(arr) == 6:
            self.x0 = float(arr[0])
            self.y0 = float(arr[1])
            self.z0 = float(arr[2])

            self.x1 = float(arr[3])
            self.y1 = float(arr[4])
            self.z1 = float(arr[5])
        else:
            self.begin = int(arr[0])
            self.end = int(arr[1])

    def get_length(self):
        return ((self.x1 - self.x0) ** 2 + (self.y1 - self.y0) ** 2 + (self.z1 - self.z0) ** 2) ** 0.5



def is_unique(arr, elem) -> bool:
    for item in arr:
        if item == elem:
            return False
    return True

x = []
y = []
r = []
z = []
basis = []
ribs = []

def read_points(path):
    f = open(path, "r")
    f.readline()
    for line in f.read().split("\n"):
        if line == '': break
        coord = line.split()
        if len(coord == 3):
            xi = float(line.split()[0])
            yi = float(line.split()[1])
            zi = float(line.split()[2])
            if is_unique(x, xi):
                x.append(xi)
            if is_unique(y, yi):
                y.append(yi)
            if is_unique(z, zi):
                z.append(zi)
        if len(coord == 2):
            ri = float(line.split()[0])
            zi = float(line.split()[1])
            if is_unique(r, ri):
                r.append(ri)
            if is_unique(z, zi):
                z.append(zi)
        f.close()


def read_basis(path_basis):
    f = open(path_basis, 'r')
    text = f.read().split("\n")
    i = 0
    it = 0
    for ri in r:
        i = 0
        basis.append([])
        for zi in z:
            basis[it].append(float(text[it * len(z) + i]))
            i+=1
        it+=1
    f.close()


def read_ribs(ribs_path):
    f = open(ribs_path, 'r')
    for line in f.read().split('\n'):
        if line == '': break
        ribs.append(rib(line.split(" ")))
    f.close()

