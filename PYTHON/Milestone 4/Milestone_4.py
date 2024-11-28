from numpy import array, linspace, transpose
from Funciones import Oscilador, Cauchy, Euler, Euler_implicito, Crank_Nicholson, RK2, RK3, RK4, Leap_Frog, Reg_estabilidad
import matplotlib.pyplot as plt

######### CONDICIONES INICIALES #########

t0 = 0
tf = 10
N = 400
t = linspace(t0, tf, N+1)

######### OSCILADOR #########

x_osc = 1
xd_osc = 0

U0 = array([x_osc, xd_osc])

######### SELECCION ESQUEMA #########

esquema = input("Que esquema numérico desea emplear? ")

if esquema == "EE":

    Esquema = Euler
    q = 1

elif esquema == "RK2":

    Esquema = RK2
    q = 2

elif esquema == "RK3":

    Esquema = RK3
    q = 3

elif esquema == "RK4":

    Esquema = RK4
    q = 4

elif esquema == "EI":

    Esquema = Euler_implicito
    q = 1

elif esquema == "CN":

    Esquema = Crank_Nicholson
    q = 2

elif esquema == "LF":

    Esquema = Leap_Frog
    q = 2

else:

    raise ValueError("Esquema no válido")

######### SELECCION ESQUEMA #########

apartado = input("Que desea calcular Órbita o Región de Estabilidad? ")

if apartado == "Orbita":

    U = Cauchy(Esquema, U0, Oscilador, t)

    plt.figure()
    plt.axis('equal')
    plt.title( r'Órbita del oscilador con {}'.format(Esquema.__name__))
    plt.xlabel('Tiempo')
    plt.ylabel('Posición')
    plt.plot(t, U[:, 0], '-b', label = "X")
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

elif apartado == "Estabilidad":

    x, y, rho  = Reg_estabilidad(Esquema,-4, 2, -4, 4, N)

    plt.axis('equal')
    plt.title( r'Región de Estabilidad de {}'.format(Esquema.__name__))
    plt.xlabel('Re')
    plt.ylabel('Im')
    plt.contour( x, y, transpose(rho), linspace(0, 1, 11))  # linspace desde G=10 hasta G=1 para pintar las diferentes curvas de nivel
    plt.grid()
    plt.show()

else:

    raise ValueError("Apartado no válido")

