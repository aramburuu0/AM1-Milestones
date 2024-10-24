from numpy import array, zeros, linspace
from Esquemas import Kepler, Oscilador, Cauchy, Euler, Euler_implicito, Crank_nicholson, RK4, Cauchy_error
import matplotlib.pyplot as plt

######### CONDICIONES INICIALES #########

x0 = 0
xf = 20
N = 1000

t = linspace(x0, xf, N+1)

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

######### RESOLUCIÓN #########

Error_EE = Cauchy_error(Euler, U0, Problema, t, q=1)
Error_EI = Cauchy_error(Euler_implicito, U0, Problema, t, q=1)
Error_CN = Cauchy_error(Crank_nicholson, U0, Problema, t, q=2)
Error_RK4 = Cauchy_error(RK4, U0, Problema, t, q=4)

######### GRÁFICAS #########

plt.figure()
plt.axis('equal')

if Problema == Kepler:

    plt.title( r'Error X de {}'.format(Problema.__name__) )
    plt.plot(t, Error_EE[:, 0], '-b', lw = 1, label = 'Error EE')
    plt.plot(t, Error_EI[:, 0], '--r', lw = 1, label = 'Error EI')
    plt.plot(t, Error_CN[:, 0], '-g', lw = 1, label = 'Error CN')
    plt.plot(t, Error_RK4[:, 0], '--y', lw = 1, label = 'Error RK4')
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Error' )
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

    plt.title( r'Error Y de {}'.format(Problema.__name__) )
    plt.plot(t, Error_EE[:, 1], '-b', lw = 1, label = 'Error EE')
    plt.plot(t, Error_EI[:, 1], '--r', lw = 1, label = 'Error EI')
    plt.plot(t, Error_CN[:, 1], '-g', lw = 1, label = 'Error CN')
    plt.plot(t, Error_RK4[:, 1], '--y', lw = 1, label = 'Error RK4')
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Error' )
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

    plt.title( r'Error Vx de {}'.format(Problema.__name__) )
    plt.plot(t, Error_EE[:, 2], '-b', lw = 1, label = 'Error EE')
    plt.plot(t, Error_EI[:, 2], '--r', lw = 1, label = 'Error EI')
    plt.plot(t, Error_CN[:, 2], '-g', lw = 1, label = 'Error CN')
    plt.plot(t, Error_RK4[:, 2], '--y', lw = 1, label = 'Error RK4')
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Error' )
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

    plt.title( r'Error Vy de {}'.format(Problema.__name__) )
    plt.plot(t, Error_EE[:, 3], '-b', lw = 1, label = 'Error EE')
    plt.plot(t, Error_EI[:, 3], '--r', lw = 1, label = 'Error EI')
    plt.plot(t, Error_CN[:, 3], '-g', lw = 1, label = 'Error CN')
    plt.plot(t, Error_RK4[:, 3], '--y', lw = 1, label = 'Error RK4')
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Error' )
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

elif Problema == Oscilador:
    
    plt.title( r'Error X de {}'.format(Problema.__name__) )
    plt.plot(t, Error_EE[:, 0], '-b', lw = 1, label = 'Error EE')
    plt.plot(t, Error_EI[:, 0], '--r', lw = 1, label = 'Error EI')
    plt.plot(t, Error_CN[:, 0], '-g', lw = 1, label = 'Error CN')
    plt.plot(t, Error_RK4[:, 0], '--y', lw = 1, label = 'Error RK4')
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Error' )
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

    plt.title( r'Error Vx de {}'.format(Problema.__name__) )
    plt.plot(t, Error_EE[:, 0], '-b', lw = 1, label = 'Error EE')
    plt.plot(t, Error_EI[:, 0], '--r', lw = 1, label = 'Error EI')
    plt.plot(t, Error_CN[:, 0], '-g', lw = 1, label = 'Error CN')
    plt.plot(t, Error_RK4[:, 0], '--y', lw = 1, label = 'Error RK4')
    plt.xlabel( 'Tiempo' )
    plt.ylabel( 'Error' )
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()


