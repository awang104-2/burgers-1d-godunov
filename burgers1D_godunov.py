#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


def riemann_solver(u_l, u_r, z):
    s = 0.5*(u_l+u_r)
    if u_l > u_r:
        if z <= s:
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
        

def trapezoidal_avg(a, b, F, K):
    total = 0
    for i in range(K):
        total = total + F(a+(b-a)*i/K) + F(a+(b-a)*(i+1)/K)
    average = total/(2*K)
    return average


def f(x):
    return np.sin(10*x)


def big_f(x):
    return -np.cos(10*x)/10


cells = [[0, 0.1], [0.1, 0.2], [0.2, 0.3],[0.3,0.4],[0.4,0.5],[0.5,0.6],[0.6,0.7],[0.7,0.8],[0.8,0.9],[0.9,1]]
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
        avg.append(np.round(trapezoidal_avg(a,b,f,K),6))
        actual_avg.append(np.round((big_f(b)-big_f(a))/(b-a),6))
        err = np.round(np.abs((trapezoidal_avg(a,b,f,K)-((big_f(b)-big_f(a))/(b-a)))/((big_f(b)-big_f(a))/(b-a))),5)*100
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
        actual_avg.append(np.round((big_f(b)-big_f(a))/(b-a),6))
        err = np.round(np.abs((Gauss_Legendre_Avg(a,b,f,K)-((big_f(b)-big_f(a))/(b-a)))/((big_f(b)-big_f(a))/(b-a))),4)*100
        error.append(err)
    print("Estimated Average:", avg)
    print("Analytical Average:", actual_avg)
    print("Error:", error)
    print('\n')
    avg = []
    actual_avg = []
    error = []
    

def flux(ustar):
    return ustar**2/2
    
    
def burgers1d_fv_godunov(U,dt,dx):
    Nx = len(U)
    F = np.zeros(Nx+1)
    Un = np.zeros(Nx)
    
    # CFL
    CFL = dt*np.max(np.abs(U))/dx
    if (CFL > 1/2): 
        raise Exception("CFL=" + str(CFL) + " is greater than 0.5.")
        
    # Solve the Riemann problem and flux for the discrete Un, except boundary cases.
    for i in range(Nx-1):
        ustar = riemann_solver(U[i],U[i+1],0)
        F[i+1] = flux(ustar)
        
    # Left boundary
    if U[0] < 0:
        ustar = riemann_solver(2*U[0]-U[1],U[0],0)
    else:
        ustar = riemann_solver(U[0],U[0],0)
    F[0] = flux(ustar)
    
    # Right boundary
    if U[Nx-1] > 0:
        ustar = riemann_solver(2*U[Nx-1]-U[Nx-2],U[Nx-1],0)
    else:
        ustar = riemann_solver(U[Nx-1],U[Nx-1],0)
    F[Nx] = flux(ustar)
    
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
    Un = burgers1d_fv_godunov(U_FV[-1],dt,dx)
    U_FV.append(Un)
    t += dt
    
x = np.linspace(0,1,Nx)
for i in range(16):
    plt.plot(x,U_FV[i*100])
plt.title('U(x) from t=0 to t=1.5 by 0.1 increments.')
plt.xlabel('x position [m]')
plt.ylabel('u(x) velocity [m/s]')
plt.show()

