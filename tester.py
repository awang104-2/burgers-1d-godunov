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
