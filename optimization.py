import numpy as np
import scipy.optimize

nu = 10 ** -3
R = 10 ** -3
m = 10 ** -3


Points = []
Velocities = []
Accelerations = []
U_constraints = []

with open("accelerations") as acc_file:
    N_curves = int(acc_file.readline())
    print("curves:", N_curves)
    for _ in range(N_curves):
        P = []
        V = []
        A = []
        U = []
        N_points = int(acc_file.readline())
        print("size of curve:", N_points)

        for _ in range(N_points):
            line =  acc_file.readline().rstrip().split(';')

            # deal with curves in 2d
            P.append(list(map(float, line[0].split()[:2])))
            v = list(map(float, line[1].split()[:2]))
            a = list(map(float, line[2].split()[:2]))

            # Stock's low for the simplest case of a sphere
            u = [a[i] * m / (6. * nu * np.pi * R) - v[i] for i in range(len(a))]

            V.append(v)
            A.append(a)
            U.append(u)
        
        Points.append(np.array(P))
        Velocities.append(np.array(V))
        Accelerations.append(np.array(A))
        U_constraints.append(np.array(U))

print(Velocities)
print(Accelerations)

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