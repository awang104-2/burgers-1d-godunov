#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from scipy.linalg import circulant
from scipy.integrate import solve_ivp
from scipy.linalg import eig, eigvals
from scipy.linalg import inv
import matplotlib.pyplot as plt

def Riemann_Solver(u_l, u_r, z):
    S = 0.5*(u_l+u_r)
    if u_l > u_r:
        if z <= S:
            return u_l
        else:
            return u_r
    else:
        if z <= u_l:
            return u_l
        elif z > u_r:
            return u_r
        else:
            return z
        
# right moving shock as expected for positive u_l > u_r
u_l = 8; u_r = 2

L = 6; T = 0.0001; Nz = 100; Nx = 100
x = np.linspace(-L, L, Nx)
u = 0
for i in range(len(x)):
    u = Riemann_Solver(u_l, u_r, x[i]/T)
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
    u = Riemann_Solver(u_l, u_r, x[i]/T)
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
    u = Riemann_Solver(u_l, u_r, x[i]/T)
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
    u = Riemann_Solver(u_l, u_r, x[i]/T)
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
    u = Riemann_Solver(u_l, u_r, x[i]/T)
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

def Trapezoidal_Avg(a,b,F,K):
    sum = 0
    for i in range(K):
        sum = sum + F(a+(b-a)*i/K) + F(a+(b-a)*(i+1)/K)
    sum = sum/(2*K)
    return sum

def f(x):
    return np.sin(10*x)

def bigF(x):
    return -np.cos(10*x)/10

cells = [[0,0.1],[0.1,0.2],[0.2,0.3],[0.3,0.4],[0.4,0.5],[0.5,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1]]
N = [1,2,3,4,5] # number of points
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
        avg.append(np.round(Trapezoidal_Avg(a,b,f,K),6))
        actual_avg.append(np.round((bigF(b)-bigF(a))/(b-a),6))
        err = np.round(np.abs((Trapezoidal_Avg(a,b,f,K)-((bigF(b)-bigF(a))/(b-a)))/((bigF(b)-bigF(a))/(b-a))),5)*100
        error.append(err)
    print("Estimated Average:", avg)
    print("Analytical Average:", actual_avg)
    print("Error:", error)
    print('\n')
    avg = []
    actual_avg = []
    error = []
    

def Gauss_Legendre_Avg(a,b,F,K):
    x_k, w_k = np.polynomial.legendre.leggauss(K)
    sum = w_k*F(a+(b-a)/2*(x_k+1))
    return np.sum(sum)/2

def f(x):
    return np.sin(10*x)

cells = [[0,0.1],[0.1,0.2],[0.2,0.3],[0.3,0.4],[0.4,0.5],[0.5,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1]]
N = [1,2,3,4,5] # number of points
avg = []

print("All errors are in percent.\n")

for j in range(5):
    K = N[j] 
    for i in range(10):
        bounds = cells[i]
        a = bounds[0]
        b = bounds[1]
        avg.append(np.round(Gauss_Legendre_Avg(a,b,f,K),6))
        actual_avg.append(np.round((bigF(b)-bigF(a))/(b-a),6))
        err = np.round(np.abs((Gauss_Legendre_Avg(a,b,f,K)-((bigF(b)-bigF(a))/(b-a)))/((bigF(b)-bigF(a))/(b-a))),4)*100
        error.append(err)
    print("Estimated Average:", avg)
    print("Analytical Average:", actual_avg)
    print("Error:", error)
    print('\n')
    avg = []
    actual_avg = []
    error = []
    
def Flux(ustar):
    return ustar**2/2
    
def Burgers1D_FV_Godunov(U,dt,dx):
    Nx = len(U)
    F = np.zeros(Nx+1)
    Un = np.zeros(Nx)
    
    # CFL
    CFL = dt*np.max(np.abs(U))/dx
    if (CFL > 1/2): 
        raise Exception("CFL=" + str(CFL) + " is greater than 0.5.")
        
    # Solve the Riemann problem and flux for the discrete Un, except boundary cases.
    for i in range(Nx-1):
        ustar = Riemann_Solver(U[i],U[i+1],0)
        F[i+1] = Flux(ustar)
        
    # Left boundary
    if U[0] < 0:
        ustar = Riemann_Solver(2*U[0]-U[1],U[0],0)
    else:
        ustar = Riemann_Solver(U[0],U[0],0)
    F[0] = Flux(ustar)
    
    # Right boundary
    if U[Nx-1] > 0:
        ustar = Riemann_Solver(2*U[Nx-1]-U[Nx-2],U[Nx-1],0)
    else:
        ustar = Riemann_Solver(U[Nx-1],U[Nx-1],0)
    F[Nx] = Flux(ustar)
    
    for i in range(Nx):
        Un[i] = U[i] - (F[i+1]-F[i])*(dt/dx)
    return Un

def u(x):
    sigma = 3/40
    return 1/10 + np.exp(-(x-1/4)**2/(2*sigma**2))

domain = [0,1]
Nx = 500
T = 1.5
cells = []
for i in range(Nx):
    cells.append([i*domain[1]/Nx,(i+1)*domain[1]/Nx])
dx = (domain[1]-domain[0])/Nx

#Q4.1 Initializing simulation with initial conditions. Uses the Gaussian-Legendre quadriture average.
avg = []
K = 5
for i in range(Nx):
    bounds = cells[i]
    a = bounds[0]
    b = bounds[1]
    avg.append(Gauss_Legendre_Avg(a,b,u,K))    
    
U_0 = avg
U_FV = [U_0]

t = 0
CFL = 0.5
dt = CFL*dx/np.max(np.abs(U_0))


while t <= T:
    Un = Burgers1D_FV_Godunov(U_FV[-1],dt,dx)
    U_FV.append(Un)
    t += dt
    
x = np.linspace(0,1,Nx)
for i in range(16):
    plt.plot(x,U_FV[i*100])
plt.title('U(x) from t=0 to t=1.5 by 0.1 increments.')
plt.xlabel('x position [m]')
plt.ylabel('u(x) velocity [m/s]')
plt.show()

