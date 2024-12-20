from numpy import array, zeros, concatenate, linspace
from numpy.linalg import norm
from scipy.optimize import newton
import matplotlib.pyplot as plt

######### PROBLEMAS A RESOLVER #########

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

######### CAUCHY #########

def Cauchy(Esquema, U0, F, t):

    U = zeros([N+1, len(U0)])
    U[0, :] = U0

    for n in range(0, N):
        
        U[n+1, :] = Esquema(U[n, :], t[n+1] - t[n], t[n], F)
    
    return U

######### CÓDIGO PROGRAMA #########

######### CONDICIONES INICIALES #########

while True:

    try:

        t0 = float(input("Introduce el valor del instante inicial t0: "))
        tf = float(input("Introduce el valor del instante final tf: "))
        N = int(input("Introduce el número de pasos N: "))
        
        if t0 >= tf:

            raise ValueError("El valor de t0 debe ser menor que el valor de tf")

        if tf <= 0:

            raise ValueError("El valor de tf debe ser mayor que 0")

        if N < 200:

            raise ValueError("El método de Newton no converge")

        break

    except ValueError as e:

        print(f"Error: {e}.")

dt = (tf-t0)/N
print(f"El valor de dt es: {dt}")

t = linspace(t0, tf, N+1)

######### KEPLER #########

x_kepler = 1
y_kepler = 0
xd_kepler = 0
yd_kepler = 1

######### OSCILADOR #########

x_osc = 1
xd_osc = 0

######### RESOLUCIÓN #########

# Problema = Oscilador
problema = input("¿Qué problema deseas resolver? ")

if problema == "Kepler":

    Problema = Kepler
    U0 = array([x_kepler, y_kepler, xd_kepler, yd_kepler])

elif problema == "Oscilador":
     
    Problema = Oscilador
    U0 = array([x_osc, xd_osc])
    
else:
    
    raise ValueError("Problema no especificado correctamente")

U_euler = Cauchy(Euler, U0, Problema, t)
U_euler_implicito = Cauchy(Euler_implicito, U0, Problema, t)
U_crank_nicholson = Cauchy(Crank_nicholson, U0, Problema, t)
U_rk2 = Cauchy(RK2, U0, Problema, t)
U_rk3 = Cauchy(RK3, U0, Problema, t)
U_rk4 = Cauchy(RK4, U0, Problema, t)

######### GRÁFICAS #########

plt.figure(figsize=(10, 7))
plt.axis('equal')

if Problema == Kepler:

    plt.title( r'{} con ($\Delta$t={})'.format(Problema.__name__,round(dt,4)) )
    plt.plot(U_euler[:,0], U_euler[:, 1], '-b', lw = 1, label ="Euler explícito")
    plt.plot(U_euler_implicito[:,0], U_euler_implicito[:, 1], linestyle='--', color='peru', lw = 1, label ="Euler implícito")
    plt.plot(U_crank_nicholson[:,0], U_crank_nicholson[:, 1], linestyle='-.', color='indigo', lw = 1, label ="Crank-Nicholson")
    plt.plot(U_rk2[:,0], U_rk2[:, 1], '-r', lw = 1, label ="Runge-Kutta 2")
    plt.plot(U_rk3[:,0], U_rk3[:, 1], '--g', lw = 1, label ="Runge-Kutta 3")
    plt.plot(U_rk4[:,0], U_rk4[:, 1], '-.y', lw = 1, label ="Runge-Kutta 4")
    plt.xlabel( 'Coordenada x' )
    plt.ylabel( 'Coordenada y' )
    plt.legend(loc='upper right', fontsize='small')
    plt.grid()
    plt.show()

    plt.title(r'Evolución de Vx con $\Delta$t={}'.format(round(dt,4)))
    plt.plot(t, U_euler[:,0], '-b', lw = 1, label ="Euler explícito")
    plt.plot(t, U_euler_implicito[:, 0], linestyle='--', color='peru', lw = 1, label ="Euler implícito")
    plt.plot(t, U_crank_nicholson[:, 0], linestyle='-.', color='indigo', lw = 1, label ="Crank-Nicholson")
    plt.plot(t, U_rk2[:,0], '-r', lw = 1, label ="Runge-Kutta 2")
    plt.plot(t, U_rk3[:,0], '--g', lw = 1, label ="Runge-Kutta 3")
    plt.plot(t, U_rk4[:,0], '-.y', lw = 1, label ="Runge-Kutta 4")
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Velocidad horizontal Vx' )
    plt.grid()
    plt.legend(loc='upper right', fontsize='small')
    plt.show()
    
    
    plt.title(r'Evolución de Vy con $\Delta$t={}'.format(round(dt,2)))
    plt.plot(t, U_euler[:,1], '-b', lw = 1, label ="Euler explícito")
    plt.plot(t, U_euler_implicito[:, 1], linestyle='--', color='peru', lw = 1, label ="Euler implícito")
    plt.plot(t, U_crank_nicholson[:, 1], linestyle='-.', color='indigo', lw = 1, label ="Crank-Nicholson")
    plt.plot(t, U_rk2[:,1], '-r', lw = 1, label ="Runge-Kutta 2")
    plt.plot(t, U_rk3[:,1], '--g', lw = 1, label ="Runge-Kutta 3")
    plt.plot(t, U_rk4[:,1], '-.y', lw = 1, label ="Runge-Kutta 4")
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Velocidad horizontal Vy' )
    plt.grid()
    plt.legend(loc='upper right', fontsize='small')
    plt.show()
    
else:

    plt.title( r'{} con ($\Delta$t={})'.format(Problema.__name__,round(dt,2)) )
    plt.plot(t, U_euler[:,0], '-b', lw = 1, label ="Euler explícito")
    plt.plot(t, U_euler_implicito[:, 0], linestyle='--', color='peru', lw = 1, label ="Euler implícito")
    plt.plot(t, U_crank_nicholson[:, 0], linestyle='-.', color='indigo', lw = 1, label ="Crank-Nicholson")
    plt.plot(t, U_rk2[:,0], '-r', lw = 1, label ="Runge-Kutta 2")
    plt.plot(t, U_rk3[:,0], '--g', lw = 1, label ="Runge-Kutta 3")
    plt.plot(t, U_rk4[:,0], '-.y', lw = 1, label ="Runge-Kutta 4")
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Posición' )
    plt.legend(loc='upper right', fontsize='small')
    plt.grid()
    plt.show()


