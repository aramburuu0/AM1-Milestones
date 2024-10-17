from numpy import array, zeros, linspace, concatenate
from numpy.linalg import norm
from scipy.optimize import newton

import matplotlib.pyplot as plt

def Kepler(U, t):

    r = U[0: 2]
    rdot = U[2: 4]

    F = concatenate((rdot, -r/norm(r)))

    return F

def Cauchy(Esquema, U0, F, t):

    U = zeros([N+1, len(U0)])
    U[0, :] = U0

    for n in range(0, N):
        
        U[n+1, :] = Esquema(U[n, :], t[n+1] - t[n], t[n], F)
    
    return U

def Euler(U, dt, t, F):

    return U + dt*F(U,t)

U0 = array([1, 0, 0, 1])
t0 = 0
tf = 20
N = 400

dt = (tf - t0)/N

t = linspace(t0, tf, N+1)

U = zeros([N+1, len(U0)])
U[0, :] = U0

    # Inputs:  
    #        Temporal_schem 
    #        U0 : initial condition at t=0
    #        F(U,t) : Function dU/dt = F(U,t); el problema que queremos resolver
    #        t : time partition t (vector of length N+1)
    #        q : order of the scheme (optional) 
    #        Tolerance: Error tolerance (optional)

def Euler_implicito(U, dt, t, F):

    def G(X):

        return X - U - dt*F(X,t)

    return newton(G, U, maxiter=500)

def Crank_Nicholson(U, dt, t, F):

    def G(X):

        return X - U - 1/2*dt*(F(U, t)+F(X, t))
    
    return newton(G, U)

U_euler = Cauchy(Euler, U0, Kepler, t)
U_implicito = Cauchy(Euler_implicito, U0, Kepler, t)
U_crank = Cauchy(Crank_Nicholson, U0, Kepler, t)

plt.plot(U_euler[:, 0],U_euler[:, 1], '-y', label="Euler Explícito")
plt.plot(U_implicito[:, 0],U_implicito[:, 1], '--g', label="Euler Implícito")
plt.plot(U_crank[:, 0],U_crank[:, 1], ':r', label="Crank-Nicholson")
plt.axis("equal")
plt.legend(loc='upper right',fontsize='small')
plt.grid()
plt.show()
