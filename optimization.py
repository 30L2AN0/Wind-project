import numpy as np
import scipy.optimize
import math


def velocity_frame(pos, vel, size):
    frame = np.zeros(size * size * 2)
    for p in pos:
        frame[(p[1] * size + p[0]) * 2] = vel[0]
        frame[(p[1] * size + p[0]) * 2 + 1] = vel[1]

    return frame

def normals_frame(pos, normals, size):
    frame = np.zeros(size * size * 2)
    for p, n in zip(pos, normals):
        frame[(p[1] * size + p[0]) * 2] = n[0]
        frame[(p[1] * size + p[0]) * 2 + 1] = n[1]

    return frame


Velocities = np.array([])
Normals = np.array([])
BodyNormals = []

file = input()
print("Opening file:", file)
print("Reading...")
with open(file) as acc_file:
    N_boundary_cells = int(acc_file.readline())
    N_inner_cells = int(acc_file.readline())
    N_curves = int(acc_file.readline())
    grid_res = int(acc_file.readline())
    print("curves:", N_curves)
    print("boundary cells:", N_boundary_cells)
    print("inner cells:", N_inner_cells)
    print("grid size: ", grid_res)

    # read normals
    for _ in range(N_boundary_cells):
        BodyNormals.append(list(map(float, acc_file.readline().split())))
    
    for _ in range(N_curves):
        N_points = int(acc_file.readline())
        print("size of curve:", N_points)
        for i in range(N_points):
            boundary_pos = []
            inner_pos = []

            # read body pos
            for _ in range(N_boundary_cells):
                boundary_pos.append(list(map(int, acc_file.readline().split())))

            for _ in range(N_inner_cells):
                inner_pos.append(list(map(int, acc_file.readline().split())))

            # read velocity (and acceleration that is not used)
            vel = list(map(float, acc_file.readline().split()))
            acc = list(map(float, acc_file.readline().split()))

            if not math.isnan(vel[0]):
                Velocities = np.concatenate((Velocities, velocity_frame(boundary_pos + inner_pos, vel, grid_res)))
                Normals = np.concatenate((Normals, normals_frame(boundary_pos, BodyNormals, grid_res)))

                # if i == 20:
                    # frame = normals_frame(boundary_pos, BodyNormals, grid_res)
                    # for i in range(grid_res):
                    #     for j in range(grid_res):
                    #         print(frame[(i * grid_res + j) * 2 + 1], end=' ')
                    #     print()

print('Finished reading, total:')
print('Velocities shape:', Velocities.shape)
print('Normals shape:', Normals.shape)


# #######################################################

# nu = 10 ** -3
# R = 10 ** -3
# m = 10 ** -3


# Points = []
# Velocities = []
# Accelerations = []
# U_constraints = []

# with open("accelerations") as acc_file:
#     N_curves = int(acc_file.readline())
#     print("curves:", N_curves)
#     for _ in range(N_curves):
#         P = []
#         V = []
#         A = []
#         U = []
#         N_points = int(acc_file.readline())
#         print("size of curve:", N_points)

#         for _ in range(N_points):
#             line =  acc_file.readline().rstrip().split(';')

#             # deal with curves in 2d
#             P.append(list(map(float, line[0].split()[:2])))
#             v = list(map(float, line[1].split()[:2]))
#             a = list(map(float, line[2].split()[:2]))

#             # Stock's low for the simplest case of a sphere
#             u = [a[i] * m / (6. * nu * np.pi * R) - v[i] for i in range(len(a))]

#             V.append(v)
#             A.append(a)
#             U.append(u)
        
#         Points.append(np.array(P))
#         Velocities.append(np.array(V))
#         Accelerations.append(np.array(A))
#         U_constraints.append(np.array(U))

# print(Velocities)
# print(Accelerations)

# V = N # given as such from C++ program

# #######################################################

# nU = U.size // 3
# A = np.zeros([nU, 3 * nU])
# B = np.zeros([nU])
# for row in range(nU):
#     A[row, 3*row : 3*(row+1)] = N[3*row : 3*(row+1)]
#     B[row] = np.dot(V[3*row:3*(row+1)], N[3*row:3*(row+1)])

# print(A @ U)
# print(B)

# cons = ({'type': 'eq', 'fun': lambda x:  A @ x - B})

# def f(x):
#     return np.sum(x**2)

# scipy.optimize.minimize(f, U, constraints=cons)