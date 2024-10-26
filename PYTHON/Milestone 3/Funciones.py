from numpy import concatenate, zeros, array, linspace
from numpy.linalg import norm
from scipy.optimize import newton

######### PROBLEMA TERMPORAL #########

######### KEPLER #########

def Kepler(U, t):

    r = U[0: 2]
    rdot = U[2: 4]

    F = concatenate((rdot, -r/norm(r)))

    return F

######### OSCILADOR #########

def Oscilador(U, t):

    x = U[0]
    xdot = U[1]

    F = array([xdot, -x])

    return F

######### ESQUEMAS NUMÉRICOS #########

######### EULER EXPLÍCITO #########

def Euler(U, dt, t, F):

    return U + dt*F(U,t)

######### EULER IMPLICITO #########

def Euler_implicito(U, dt, t, F):

    def G(X):

        return X - U - dt*F(X,t)

    return newton(G, U, maxiter=1500)

######### RUNGE-KUTTA #########

def Crank_nicholson(U, dt, t, F):

    def G(X):

        return X - U - (dt/2) *(F(U,t) + F(X,t))
    
    return newton(G, U, maxiter=1500)

######### RUNGE-KUTTA 2 (RK2) #########

def RK2(U, dt, t, F):

    k1 = F(U, t)
    k2 = F(U + dt*k1, t + dt)

    return U + dt/2 * (k1 + k2)

######### RUNGE-KUTTA 3 (RK3) #########

def RK3(U, dt, t, F):

    k1 = F(U, t)
    k2 = F(U + dt*k1/3, t + dt/3)
    k3 = F(U + 2*dt*k2/3, t + 2*dt/3)

    return U + dt/4 * (k1 + 3*k3)

######### RUNGE-KUTTA 4 (RK4) #########

def RK4(U, dt, t, F):

    k1 = F(U, t)
    k2 = F(U + dt*k1/2, t + dt/2)
    k3 = F(U + dt*k2/2, t + dt/2)
    k4 = F(U + dt*k3, t + dt)

    return U + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

######### METODOS DE RESOLUCIÓN #########

######### CAUCHY #########

def Cauchy(Esquema, U0, F, t):

    U = zeros([len(t), len(U0)])
    U[0, :] = U0

    for n in range(0, len(t)-1):
        
        U[n+1, :] = Esquema(U[n, :], t[n+1] - t[n], t[n], F)
    
    return U

######### MÉTODO DE RICHARDSON #########

def Cauchy_error(Esquema, U0, F, t, q):

    N = len(t) - 1
    t1 = t
    t2 = linspace(t[0], t[-1], 2*N+1)    # hace una partición equispaciada desde t0 hasta tf con el doble de nodos de t, si t1 tiene N+1 nodos, t2 tiene 2*N+1
    
    Error = zeros([len(t1), len(U0)])
    U1 = Cauchy(Esquema, U0, F, t1)         # solución de Cauchy al esquema para la malla original
    U2 = Cauchy(Esquema, U0, F, t2)         # solución de Cauchy al esquema para la malla refinada

    for n in range(0, N):

        Error[n, :] = (U2[2*n, :]-U1[n, :])/(1-1/2**q)

    return Error

