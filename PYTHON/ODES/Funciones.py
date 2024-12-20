from numpy import concatenate, zeros, array, linspace, exp, log10
from numpy.linalg import norm
from scipy.optimize import newton
from cmath import sqrt

######### PROBLEMA TERMPORAL #########

######### KEPLER #########

def Kepler(U, t):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal

    '''''''''''

    r = U[0: 2]
    rdot = U[2: 4]

    F = concatenate((rdot, -r/norm(r)))

    return F

######### OSCILADOR #########

def Oscilador(U, t):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal

    '''''''''''

    x = U[0]
    xdot = U[1]

    F = array([xdot, -x])

    return F

######### ESQUEMAS NUMÉRICOS #########

######### EULER EXPLÍCITO #########

def Euler(U, dt, t, F):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal
    dt: Salto temporal en cada paso dt=T/N
    F: Función a resolver

    '''''''''''

    return U + dt*F(U,t)

######### EULER IMPLICITO #########

def Euler_implicito(U, dt, t, F):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal
    dt: Salto temporal en cada paso dt=T/N
    F: Función a resolver

    '''''''''''

    def G(X):

        return X - U - dt*F(X,t)

    return newton(G, U, maxiter=3000)

######### RUNGE-KUTTA #########

def Crank_Nicholson(U, dt, t, F):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal
    dt: Salto temporal en cada paso dt=T/N
    F: Función a resolver

    '''''''''''

    def G(X):

        return X - U - (dt/2) *(F(U,t) + F(X,t))
    
    return newton(G, U, maxiter=1500)

######### RUNGE-KUTTA 2 (RK2) #########

def RK2(U, dt, t, F):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal
    dt: Salto temporal en cada paso dt=T/N
    F: Función a resolver

    '''''''''''

    k1 = F(U, t)
    k2 = F(U + dt*k1, t + dt)

    return U + dt/2 * (k1 + k2)

######### RUNGE-KUTTA 3 (RK3) #########

def RK3(U, dt, t, F):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal
    dt: Salto temporal en cada paso dt=T/N
    F: Función a resolver

    '''''''''''

    k1 = F(U, t)
    k2 = F(U + dt*k1/3, t + dt/3)
    k3 = F(U + 2*dt*k2/3, t + 2*dt/3)

    return U + dt/4 * (k1 + 3*k3)

######### RUNGE-KUTTA 4 (RK4) #########

def RK4(U, dt, t, F):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    t: Partición temporal
    dt: Salto temporal en cada paso dt=T/N
    F: Función a resolver

    '''''''''''

    k1 = F(U, t)
    k2 = F(U + dt*k1/2, t + dt/2)
    k3 = F(U + dt*k2/2, t + dt/2)
    k4 = F(U + dt*k3, t + dt)

    return U + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

######### LEAP-FROG (LP) #########

def Leap_Frog(U, Uant, dt, t, F):

    '''''''''''
    --INPUTS--

    U: Vector de estado (posición y velocidad)
    Uant: Vector de estado de N-1 (Euler)
    t: Partición temporal
    F: Función a resolver

    '''''''''''

    return Uant + 2*dt*F(U,t)

######### METODOS DE RESOLUCIÓN #########

######### CAUCHY #########

def Cauchy(Esquema, U0, F, t):

    '''''''''''
    --INPUTS--

    Esquema (U, dt, t, F): Esquema númerico para resolver el problema
    U0: Vector de condiciones iniciales del problema a resolver
    F(U, t): Función a resolver
    t: Partición temporal

    '''''''''''

    U = zeros([len(t), len(U0)])
    U[0, :] = U0
    dt = t[1] - t[0]

    if Esquema == Leap_Frog:

        U[1, :] = U[0, :] + dt * F(U[0, :], t[0])

        for n in range(1, len(t)-1):

            U[n+1,:] = Leap_Frog(U[n, :], U[n-1, :], dt ,t[n], F)

        return U

    else:

        for n in range(0, len(t)-1):

            U[n+1,:] = Esquema(U[n, :], t[n+1] - t[n], t[n], F)

        return U

######### MÉTODO DE RICHARDSON #########

def Cauchy_error(Esquema, U0, F, t, q):

    '''''''''''
    --INPUTS--

    Esquema (U, dt, t, F): Esquema númerico para resolver el problema
    U0: Vector de condiciones iniciales dl problema a resolver
    F(U, t): Función a resolver
    t: Partición temporal
    q: Orden del esquema númerico empleado 

    '''''''''''

    N = len(t) - 1
    t1 = t
    t2 = linspace(t[0], t[-1], 2*N+1)    # hace una partición equispaciada desde t0 hasta tf con el doble de nodos de t, si t1 tiene N+1 nodos, t2 tiene 2*N+1
    
    Error = zeros([len(t1), len(U0)])
    U1 = Cauchy(Esquema, U0, F, t1)         # solución de Cauchy al esquema para la malla original
    U2 = Cauchy(Esquema, U0, F, t2)         # solución de Cauchy al esquema para la malla refinada

    for n in range(0, N+1):

        Error[n, :] = (U2[2*n, :]-U1[n, :])/(1-1/2**q)

    return Error

def Cauchy_error2(Esquema, U0, F, t):

    '''''''''''
    --INPUTS--

    Esquema (U, dt, t, F): Esquema númerico para resolver el problema
    U0: Vector de condiciones iniciales dl problema a resolver
    F(U, t): Función a resolver
    t: Partición temporal 

    '''''''''''

    N = len(t) - 1
    t1 = t
    t2 = linspace(t[0], t[-1], 2*N+1)    # hace una partición equispaciada desde t0 hasta tf con el doble de nodos de t, si t1 tiene N+1 nodos, t2 tiene 2*N+1
    
    Error2 = zeros([len(t1), len(U0)])
    U1 = Cauchy(Esquema, U0, F, t1)         # solución de Cauchy al esquema para la malla original
    U2 = Cauchy(Esquema, U0, F, t2)         # solución de Cauchy al esquema para la malla refinada

    for n in range(0, N+1):

        Error2[n, :] = (U2[2*n, :]-U1[n, :])

    return Error2

######### CONVERGENCIA #########

def Convergencia(Esquema, U0, F, t, Error, Cauchy, Ptosgraf):

    '''''''''''
    Inputs:

        Esquema: Método numérico
        U0: Vector de estado inicial
        F: Función a resolver
        t: particion temporal
        Error(Esquema, U0, F, t): Función que devuelve un vector con el error de un esquema en cada paso temporal
        Cauchy: Problema del valor inicial Cauchy

    '''''''''''

    logE = zeros(Ptosgraf)
    logN = zeros(Ptosgraf)
    t1 = t
    N = len(t)-1
    
    for n in range(Ptosgraf):
        
        E = Error(Esquema, U0, F, t1)
        logE[n] = log10(norm(E[-1, :]))
        logN[n] = log10(N)

        N=2*N
        t1 = linspace(t[0], t[-1], N+1)

    return logN, logE

######### REGIONES DE ESTABILIDAD #########

def Reg_estabilidad(Esquema, x0, xf, y0, yf, N):

    '''''''''''
    Inputs:

        Esquema: Esquema númerico del que se obtiene la región de estabilidad
        x0, xf: Puntos inicial y final de la malla Re
        y0, yf: Puntos inicial y final de la malla Im
        N: Número de particiones de la malla
        
    '''''''''''

    x = linspace(x0, xf, N)
    y = linspace(y0, yf, N)

    rho = zeros((N,N))

    for i in range(N):

        for j in range(N):

            w = complex(x[i], y[j])

            if Esquema == Leap_Frog:

                r1 = (2*w + sqrt(4*w**2+4))/2
                r2 = (2*w - sqrt(4*w**2+4))/2
                r = max(abs(r1), abs(r2))

            else:
                
                r = Esquema(1, 1, 0, lambda u, t : u*w)

            rho[i,j] = abs(r)

    return x, y, rho

######### NEWTON #########

def Newton(F, x_0, Fprima = None, tol = 1e-8, maxiter=50):

    '''''''''''
    --INPUTS--

        F: Función escalar de la que sacar las raíces
        x_0: Punto inicial del eje x en el que se comienza la iteración  
        Fprima: Derivada de F. (Si no se introduce se calcula dentro de la función)
        tol: Tolerancia (por defecto es 10e-8)
        maxiter: Número máximo de iteraciones

    '''''''''''
    
    xn = x_0
    Error = tol + 1
    iter = 0

    def Fp(x):

        if Fp==None:

            delta = 1e-4

            return (F(x+delta)-F(x-delta))/(2*delta)
        
        else:

            return Fprima(x)

    while Error > tol and iter < maxiter:

        xn1 = xn - F(xn)/Fp(xn)
        Error = abs(xn-xn1)
        xn = xn1

        iter += 1

        if iter >= maxiter:

            print(f'Max number of iterations reached ({maxiter}). Returning...')
            
            return xn

    print('Número de iterciones =', iter)

    return xn

def jorge(x):

    return exp(x)-2*x-2 # Al usar la exp de numpy, se aceptan inputs tanto como escalar o vectorial (lista)

def dif_jorge(x):

    return exp(x)-2
