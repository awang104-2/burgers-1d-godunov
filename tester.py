import matplotlib.pyplot as plt
import burgers1D_godunov as bg
import numpy as np

# right moving shock as expected for positive u_l > u_r
u_l = 8; u_r = 2
L = 6; T = 0.0001; Nz = 100; Nx = 100
x = np.linspace(-L, L, Nx)
u = 0
for i in range(len(x)):
    u = bg.riemann_solver(u_l, u_r, x[i]/T)
    t_i = np.linspace(0, 1, 100)
    x_i = t_i*u - T*u + x[i]
    if i%5 == 0:
        if u == u_l:
            plt.plot(x_i,t_i,color='red',linestyle='--')
        else:
            plt.plot(x_i,t_i,color='blue',linestyle='--')
x = np.linspace(0, 5, Nx)
S = (u_l+u_r)/2
plt.plot(x,x/S,color='black')
plt.title("$U_L$ > $U_R$ and $U_L$, $U_R$ > 0")
plt.xlabel("x")
plt.ylabel("t")
plt.show()

# left moving shock as expected for negative u_l > u_r
u_l = -2; u_r = -8
L = 6; T = 0.0001; Nz = 100; Nx = 100
x = np.linspace(-L, L, Nx)
u = 0
for i in range(len(x)):
    u = bg.riemann_solver(u_l, u_r, x[i]/T)
    t_i = np.linspace(0, 1, 100)
    x_i = t_i*u - T*u + x[i]
    if i%5 == 0:
        if u == u_l:
            plt.plot(x_i,t_i,color='red',linestyle='--')
        else:
            plt.plot(x_i,t_i,color='blue',linestyle='--')
x = np.linspace(-5, 0, Nx)
S = (u_l+u_r)/2
plt.plot(x,x/S,color='black')
plt.title("$U_L$ > $U_R$ and $U_L$, $U_R$ < 0")
plt.xlabel("x")
plt.ylabel("t")
plt.show()

# right moving fan in the case of positive u_l < u_r as expected
u_l = 2; u_r = 8
L = 6; T = 1; Nz = 100; Nx = 100
x = np.linspace(-L, 2*L, Nx)
u = 0
for i in range(len(x)):
    u = bg.riemann_solver(u_l, u_r, x[i]/T)
    t_i = np.linspace(0, T, 100)
    x_i = t_i*u - T*u + x[i]
    if i%5 == 0:
        if u == u_l:
            plt.plot(x_i,t_i,color='red',linestyle='--')
        else:
            plt.plot(x_i,t_i,color='blue',linestyle='--')
plt.title("$U_L$ < $U_R$ and $U_L$, $U_R$ > 0")
plt.xlabel("x")
plt.ylabel("t")
plt.show()

# left moving fan in the case of positive u_l < u_r as expected
u_l = -8; u_r = -2
L = 6; T = 1; Nz = 100; Nx = 100
x = np.linspace(-2*L, L, Nx)
u = 0
for i in range(len(x)):
    u = bg.riemann_solver(u_l, u_r, x[i]/T)
    t_i = np.linspace(0, T, 100)
    x_i = t_i*u - T*u + x[i]
    if i%5 == 0:
        if u == u_l:
            plt.plot(x_i,t_i,color='red',linestyle='--')
        else:
            plt.plot(x_i,t_i,color='blue',linestyle='--')
plt.title("$U_L$ < $U_R$ and $U_L$, $U_R$ < 0")
plt.xlabel("x")
plt.ylabel("t")
plt.show()

# centered fan in the case of positive u_l < u_r as expected
u_l = -8; u_r = 2
L = 6; T = 1; Nz = 100; Nx = 100
x = np.linspace(-2*L, L, Nx)
u = 0
for i in range(len(x)):
    u = bg.riemann_solver(u_l, u_r, x[i]/T)
    t_i = np.linspace(0, T, 100)
    x_i = t_i*u - T*u + x[i]
    if i%5 == 0:
        if u == u_l:
            plt.plot(x_i,t_i,color='red',linestyle='--')
        else:
            plt.plot(x_i,t_i,color='blue',linestyle='--')
plt.title("$U_L$ < $U_R$ and $U_L$ < 0, $U_R$ > 0")
plt.xlabel("x")
plt.ylabel("t")
plt.show()

cells = [[0, 0.1], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.5, 0.6], [0.6, 0.7], [0.7, 0.8], [0.8, 0.9],
         [0.9, 1]]
N = [1, 2, 3, 4, 5]  # number of points
avg = []
actual_avg = []
error = []

print("All errors are in percent.\n")

for j in range(5):
    K = N[j]
    for i in range(10):
        bounds = cells[i]
        a = bounds[0]
        b = bounds[1]
        avg.append(np.round(bg.trapezoidal_avg(a, b, bg.f, K), 6))
        actual_avg.append(np.round((bg.big_f(b) - bg.big_f(a)) / (b - a), 6))
        err = np.round(np.abs((bg.trapezoidal_avg(a, b, bg.f, K) - ((bg.big_f(b) - bg.big_f(a)) / (b - a))) / ((bg.big_f(b) - bg.big_f(a)) / (b - a))), 5) * 100
        error.append(err)
    print("Estimated Average:", avg)
    print("Analytical Average:", actual_avg)
    print("Error:", error)
    print('\n')
    avg = []
    actual_avg = []
    error = []

cells = [[0, 0.1], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.5, 0.6], [0.6, 0.7], [0.7, 0.8], [0.8, 0.9],
         [0.9, 1]]
N = [1, 2, 3, 4, 5]  # number of points
avg = []

print("All errors are in percent.\n")

for j in range(5):
    K = N[j]
    for i in range(10):
        bounds = cells[i]
        a = bounds[0]
        b = bounds[1]
        avg.append(np.round(bg.Gauss_Legendre_Avg(a, b, bg.f, K), 6))
        actual_avg.append(np.round((bg.big_f(b) - bg.big_f(a)) / (b - a), 6))
        err = np.round(np.abs((bg.Gauss_Legendre_Avg(a, b, bg.f, K) - ((bg.big_f(b) - bg.big_f(a)) / (b - a))) / ((bg.big_f(b) - bg.big_f(a)) / (b - a))), 4) * 100
        error.append(err)
    print("Estimated Average:", avg)
    print("Analytical Average:", actual_avg)
    print("Error:", error)
    print('\n')
    avg = []
    actual_avg = []
    error = []

domain = [0, 1]
Nx = 500
T = 1.5
cells = []
for i in range(Nx):
    cells.append([i * domain[1] / Nx, (i + 1) * domain[1] / Nx])
dx = (domain[1] - domain[0]) / Nx

# Q4.1 Initializing simulation with initial conditions. Uses the Gaussian-Legendre quadrature average.
avg = []
K = 5
for i in range(Nx):
    bounds = cells[i]
    a = bounds[0]
    b = bounds[1]
    avg.append(bg.Gauss_Legendre_Avg(a, b, bg.u, K))

U_0 = avg
U_FV = [U_0]

t = 0
CFL = 0.5
dt = CFL * dx / np.max(np.abs(U_0))

while t <= T:
    Un = bg.burgers1d_fv_godunov(U_FV[-1], dt, dx)
    U_FV.append(Un)
    t += dt

x = np.linspace(0, 1, Nx)
for i in range(16):
    plt.plot(x, U_FV[i * 100])
plt.title('U(x) from t=0 to t=1.5 by 0.1 increments.')
plt.xlabel('x position [m]')
plt.ylabel('u(x) velocity [m/s]')
plt.show()

