import numpy as np
from pathlib import Path


def read_mesh(path):
    path = Path(path)

    with path.open() as f:
        V, F = read_buffer(f)
    return V, F


def read_buffer(f):
    V = []
    N = []
    TC = []
    F = []
    FTC = []
    FN = []
    for line_no, line in enumerate(f):
        ltype = line[0:2]
        if ltype == 'v ':
            x = line[2:].split()
            count = len(x)
            if count != 3 and count != 4:
                raise Exception("Error: readOBJ() vertex on line %d should have 3 or 4 coordinates", line_no)
            x = np.asarray(x, dtype=np.float64)
            V.append(x)
        elif ltype == 'vn':
            x = line[2:].split()
            count = len(x)
            if count != 3:
                raise Exception("Error: readOBJ() normal on line %d should have 3 coordinates", line_no)
            x = np.asarray(x, dtype=np.float64)
            N.append(x)
        elif ltype == 'vt':
            x = line[2:].split()
            count = len(x)
            if count != 2 and count != 3:
                raise Exception("Error: readOBJ() texture coords on line %d should have 2 or 3 coordinates (%d)",
                                line_no, count)
            x = np.asarray(x, dtype=np.float64)
            TC.append(x)
        elif ltype == 'f ':
            x = line[1:].split()
            count = len(x)
            if count < 3:
                raise Exception("Error: readOBJ() face on line %d should have 3 or more vertices (%d)", line_no,
                                count)
            vertices = []
            textures = []
            normals = []
            for v in x:
                v = v.split("/")
                count = len(v)
                if count == 1:  # f v1 v2 v3 ...
                    vertices.append(v[0])
                elif count == 2:  # f v1/vt1 v2/vt2 v3/vt3 ...
                    vertices.append(v[0])
                    textures.append(v[1])
                elif count == 3:  # f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
                    vertices.append(v[0])
                    if v[1] != '':  # f v1//vn1 v2//vn2 v3//vn3
                        textures.append(v[1])
                    normals.append(v[2])

            if ((len(vertices) == len(textures)) or len(textures) == 0) and \
                    ((len(vertices) == len(normals)) or len(normals) == 0):
                F.append(vertices)
                FTC.append(textures)
                FN.append(normals)
            else:
                raise Exception("Error: readOBJ() face on line %d has invalid format\n", line_no)

        elif ltype == 'l':
            x = line[1:].split()
            vertices = [int(v) for v in x]
            F.append(vertices)

        elif ltype == '#':
            pass

    V = np.array(V, dtype=np.float64)

    # allow for diff size faces (inefficient)
    nmax = np.max([len(x) for x in F])
    F = list(map(lambda x: x + [-1]*(nmax-len(x)), F))
    F = np.array(F, dtype=np.int32)
    # obj starts at 1
    F -= 1
    return V, F
