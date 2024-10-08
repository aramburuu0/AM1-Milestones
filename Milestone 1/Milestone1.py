from numpy import array,zeros
import matplotlib.pyplot as plot

def Euler_method():

    U = array([ 1, 0, 0, 1])

    N = 1000
    x = zeros(N)
    y = zeros(N)
    t = zeros(N)
    x[0] = U[0]
    y[0] = U[1]
    t[0] = 0

    for i in range(1,N):

        deltat = 0.01
        t[i] = deltat*i
        F = Kepler( U, t[i-1])

        U = U + deltat*F
        x[i] = U[0]
        y[i] = U[1]
    
    plot.plot(x, y)
    plot.show()


def Kepler(U, t,):
        
    x = U[0]; y = U[1]; xd = U[2]; yd = U[3]
    modulo = (x**2+y**2)**(1/2)

    return array([ xd, yd, -x/modulo**3, -y/modulo**3])
    
Euler_method()
