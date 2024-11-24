from numpy import array, zeros, linspace, polyfit
from Funciones import Kepler, Oscilador, Cauchy, Euler, Euler_implicito, Crank_Nicholson, RK4, Cauchy_error, Cauchy_error2, Convergencia
import matplotlib.pyplot as plt


######### CONDICIONES INICIALES #########

t0 = 0
tf = 20
N = 250

t = linspace(t0, tf, N+1)

######### KEPLER #########

x_kepler = 1
y_kepler = 0
xd_kepler = 0
yd_kepler = 1

######### OSCILADOR #########

x_osc = 1
xd_osc = 0

problema = input("¿Qué problema deseas resolver? ")

if problema == "Kepler":

    Problema = Kepler
    U0 = array([x_kepler, y_kepler, xd_kepler, yd_kepler])

elif problema == "Oscilador":
     
    Problema = Oscilador
    U0 = array([x_osc, xd_osc])

else:
    
    raise ValueError("Problema no especificado correctamente")

esquema = input("Que esquema numérico desea emplear? ")

if esquema == "EE":

    Esquema = Euler
    q = 1
    N_pts = 10

elif esquema == "RK4":

    Esquema = RK4
    q = 4
    N_pts = 7

elif esquema == "EI":

    Esquema = Euler_implicito
    q = 1
    N_pts = 10

elif esquema == "CN":

    Esquema = Crank_Nicholson
    q = 2
    N_pts = 9

else:

    raise ValueError("Esquema no válido")

######### RESOLUCIÓN Y GRÁFICAS #########

apartado = input("Que desea caclular Error o Convergencia? ")

if apartado == "Error":

    Error = Cauchy_error(Esquema, U0, Problema, t, q)
    
    plt.figure()
    plt.axis('equal')
    plt.title( r'Error de {} con {}'.format(Problema.__name__, Esquema.__name__))
    plt.xlabel('Tiempo')
    plt.ylabel('Error')

    if Problema == Kepler:
    
        plt.plot(t, Error[:,0], '-b', label = 'X')
        plt.plot(t, Error[:,1], '--r', label = 'Y')
        plt.plot(t, Error[:,2], '-g', label = 'Vx')
        plt.plot(t, Error[:,3], '--y', label = 'Vy')

    elif Problema == Oscilador:

        plt.plot(t, Error[:,0], '-b', label = 'X')
        plt.plot(t, Error[:,1], '--r', label = 'Vx')
    
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

elif apartado == "Convergencia":

    logN, logE = Convergencia(Esquema, U0, Problema, t, Cauchy_error2, Cauchy,N_pts) # Ajustar N en función del esquema para ver únicamente la parte recta

    coef = polyfit(logN, logE, 1)
    m, b = coef
    y_fit = m * logN + b

    plt.figure()
    plt.axis('equal')
    plt.title( r'Error de {} con {}'.format(Problema.__name__, Esquema.__name__))
    plt.xlabel('logN')
    plt.ylabel('logE')
    plt.plot(logN, logE, 'bo', label='Puntos de convergencia')
    # plt.plot(logN, y_fit, '-r', label=f'Regresión Lineal: y = {m:.1f}x + {b:.1f}')
    plt.legend(loc='upper right',fontsize='small')
    plt.grid()
    plt.show()

else:

    raise ValueError("Apartado no especificado correctamente")

