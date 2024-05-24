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


def Gauss_Legendre_Avg(a, b, fun, k):
    x_k, w_k = np.polynomial.legendre.leggauss(k)
    sum = w_k*fun(a+(b-a)/2*(x_k+1))
    return np.sum(sum)/2


def f(x):
    return np.sin(10*x)


def flux(ustar):
    return ustar ** 2 / 2


def burgers1d_fv_godunov(U, dt, dx):
    Nx = len(U)
    F = np.zeros(Nx + 1)
    Un = np.zeros(Nx)

    # CFL
    CFL = dt * np.max(np.abs(U)) / dx
    if CFL > 1/2:
        raise Exception("CFL=" + str(CFL) + " is greater than 0.5.")

    # Solve the Riemann problem and flux for the discrete Un, except boundary cases.
    for i in range(Nx - 1):
        ustar = riemann_solver(U[i], U[i + 1], 0)
        F[i + 1] = flux(ustar)

    # Left boundary
    if U[0] < 0:
        ustar = riemann_solver(2 * U[0] - U[1], U[0], 0)
    else:
        ustar = riemann_solver(U[0], U[0], 0)
    F[0] = flux(ustar)

    # Right boundary
    if U[Nx - 1] > 0:
        ustar = riemann_solver(2 * U[Nx - 1] - U[Nx - 2], U[Nx - 1], 0)
    else:
        ustar = riemann_solver(U[Nx - 1], U[Nx - 1], 0)
    F[Nx] = flux(ustar)

    for i in range(Nx):
        Un[i] = U[i] - (F[i + 1] - F[i]) * (dt / dx)
    return Un


def u(x):
    sigma = 3 / 40
    return 1 / 10 + np.exp(-(x - 1 / 4) ** 2 / (2 * sigma ** 2))
