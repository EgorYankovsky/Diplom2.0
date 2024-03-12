import os

def isUnique(arr, elem) -> bool:
    for item in arr:
        if item == elem:
            return False
    return True

r = []
z = []
basis = []

def read_points(path):
    f = open(path, "r")
    f.readline()
    for line in f.read().split("\n"):
        if line == '':
            break
        ri = float(line.split()[0])
        zi = float(line.split()[1])
        if isUnique(r, ri):
            r.append(ri)
        if isUnique(z, zi):
            z.append(zi)
    f.close()


def read_basis(pathBasis):
    f = open(pathBasis, "r")
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