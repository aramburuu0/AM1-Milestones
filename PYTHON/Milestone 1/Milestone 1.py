from numpy import zeros, array, linspace
import matplotlib.pyplot as plt

######### CONDICIONES INICIALES ÓRBITA DE KEPLER #########

x_kepler = 1
y_kepler = 0
xd_kepler = 0
yd_kepler = 1

t0 = 0
tf = 20
N = 200

deltat= (tf-t0)/N

U0_kepler = array([ x_kepler, y_kepler, xd_kepler, yd_kepler ])

######### CODIGO PARA LOS DIFERENTES ESQUEMAS TEMPORALES #########

######### EULER EXPLICITO (EE) #########

F_euler = zeros([N+1,len(U0_kepler)])
U_euler = zeros([N+1,len(U0_kepler)])
t = zeros(N+1)

U_euler[0, :] = U0_kepler

for n in range(0, N):

    F_euler[n, 0] = U_euler[n, 2]
    F_euler[n, 1] = U_euler[n, 3]
    F_euler[n, 2] = -U_euler[n, 0]/(U_euler[n, 0]**2 + U_euler[n, 1]**2)**(3/2)
    F_euler[n, 3] = -U_euler[n, 1]/(U_euler[n, 0]**2 + U_euler[n, 1]**2)**(3/2)

    t[n+1] = t[n] + deltat

    U_euler[n+1, :] = U_euler[n, :] + (t[n+1] - t[n]) * F_euler[n, :]

######### RUNGE-KUTTA 2 ETAPAS (RK2) #########

U_rk2 = zeros([N+1,len(U0_kepler)])
F_rk2 = zeros([N+1,len(U0_kepler)])
k1_rk2 = zeros([N+1,len(U0_kepler)])
k2_rk2 = zeros([N+1,len(U0_kepler)])

U_rk2[0, :] = U0_kepler

for n in range(0, N):

    F_rk2[n, 0] = U_rk2[n, 2]
    F_rk2[n, 1] = U_rk2[n, 3]
    F_rk2[n, 2] = -U_rk2[n, 0]/(U_rk2[n, 0]**2 + U_rk2[n, 1]**2)**(3/2)
    F_rk2[n, 3] = -U_rk2[n, 1]/(U_rk2[n, 0]**2 + U_rk2[n, 1]**2)**(3/2)
     
    k1_rk2[n, :] = F_rk2[n, :]

    k2_rk2[n, 0] = U_rk2[n, 2] + k1_rk2[n, 2] * (t[n+1] - t[n])
    k2_rk2[n, 1] = U_rk2[n, 3] + k1_rk2[n, 3] * (t[n+1] - t[n])
    k2_rk2[n, 2] = -(U_rk2[n, 0] + k1_rk2[n, 0] * (t[n+1] - t[n])) / ((U_rk2[n, 0] + k1_rk2[n, 0] * (t[n+1] - t[n]))**2 + (U_rk2[n, 1] + k1_rk2[n, 1] * (t[n+1] - t[n]))**2)**(3/2)
    k2_rk2[n, 3] = -(U_rk2[n, 1] + k1_rk2[n, 1] * (t[n+1] - t[n])) / ((U_rk2[n, 0] + k1_rk2[n, 0] * (t[n+1] - t[n]))**2 + (U_rk2[n, 1] + k1_rk2[n, 1] * (t[n+1] - t[n]))**2)**(3/2)

    U_rk2[n+1, :] = U_rk2[n, :] + deltat/2 * (k1_rk2[n, :] +  k2_rk2[n, :])

######### RUNGE-KUTTA 4 ETAPAS (RK4) #########

U_rk4 = zeros([N+1,len(U0_kepler)])
F_rk4 = zeros([N+1,len(U0_kepler)])
k1_rk4 = zeros([N+1,len(U0_kepler)])
k2_rk4 = zeros([N+1,len(U0_kepler)])
k3_rk4 = zeros([N+1,len(U0_kepler)])
k4_rk4 = zeros([N+1,len(U0_kepler)])

U_rk4[0, :] = U0_kepler

for n in range (0, N):

    F_rk4[n, 0] = U_rk4[n, 2]
    F_rk4[n, 1] = U_rk4[n, 3]
    F_rk4[n, 2] = -U_rk4[n, 0]/(U_rk4[n, 0]**2 + U_rk4[n, 1]**2)**(3/2)
    F_rk4[n, 3] = -U_rk4[n, 1]/(U_rk4[n, 0]**2 + U_rk4[n, 1]**2)**(3/2)
     
    k1_rk4[n, :] = F_rk4[n, :]

    k2_rk4[n, 0] = U_rk4[n, 2] + k1_rk4[n, 2] * (t[n+1] - t[n])/2
    k2_rk4[n, 1] = U_rk4[n, 3] + k1_rk4[n, 3] * (t[n+1] - t[n])/2
    k2_rk4[n, 2] = -(U_rk4[n, 0] + k1_rk4[n, 0] * (t[n+1] - t[n])/2) / ((U_rk4[n, 0] + k1_rk4[n, 0] * (t[n+1] - t[n])/2)**2 + (U_rk4[n, 1] + k1_rk4[n, 1] * (t[n+1] - t[n])/2)**2)**(3/2)
    k2_rk4[n, 3] = -(U_rk4[n, 1] + k1_rk4[n, 1] * (t[n+1] - t[n])/2) / ((U_rk4[n, 0] + k1_rk4[n, 0] * (t[n+1] - t[n])/2)**2 + (U_rk4[n, 1] + k1_rk4[n, 1] * (t[n+1] - t[n])/2)**2)**(3/2)

    k3_rk4[n, 0] = U_rk4[n, 2] + k2_rk4[n, 2] * (t[n+1] - t[n])/2
    k3_rk4[n, 1] = U_rk4[n, 3] + k2_rk4[n, 3] * (t[n+1] - t[n])/2
    k3_rk4[n, 2] = -(U_rk4[n, 0] + k2_rk4[n, 0] * (t[n+1] - t[n])/2) / ((U_rk4[n, 0] + k2_rk4[n, 0] * (t[n+1] - t[n])/2)**2 + (U_rk4[n, 1] + k2_rk4[n, 1] * (t[n+1] - t[n])/2)**2)**(3/2)
    k3_rk4[n, 3] = -(U_rk4[n, 1] + k2_rk4[n, 1] * (t[n+1] - t[n])/2) / ((U_rk4[n, 0] + k2_rk4[n, 0] * (t[n+1] - t[n])/2)**2 + (U_rk4[n, 1] + k2_rk4[n, 1] * (t[n+1] - t[n])/2)**2)**(3/2)

    k4_rk4[n, 0] = U_rk4[n, 2] + k3_rk4[n, 2] * (t[n+1] - t[n])
    k4_rk4[n, 1] = U_rk4[n, 3] + k3_rk4[n, 3] * (t[n+1] - t[n])
    k4_rk4[n, 2] = -(U_rk4[n, 0] + k3_rk4[n, 0] * (t[n+1] - t[n])) / ((U_rk4[n, 0] + k3_rk4[n, 0] * (t[n+1] - t[n]))**2 + (U_rk4[n, 1] + k3_rk4[n, 1] * (t[n+1] - t[n]))**2)**(3/2)
    k4_rk4[n, 3] = -(U_rk4[n, 1] + k3_rk4[n, 1] * (t[n+1] - t[n])) / ((U_rk4[n, 0] + k3_rk4[n, 0] * (t[n+1] - t[n]))**2 + (U_rk4[n, 1] + k3_rk4[n, 1] * (t[n+1] - t[n]))**2)**(3/2)

    U_rk4[n+1, :] = U_rk4[n, :] + deltat/6 * (k1_rk4[n, :] + 2 * k2_rk4[n, :] + 2 * k3_rk4[n, :] + k4_rk4[n, :])

######### GRAFICAS #########

plt.figure(figsize=(8, 5))
plt.axis("equal")

plt.plot( U_euler[:, 0], U_euler[:,1], '-b', lw = 1, label ="Euler explícito" )
plt.plot( U_rk2[:, 0], U_rk2[:,1], '-r', lw = 1, label ="Runge-Kutta 2" )
plt.plot( U_rk4[:, 0], U_rk4[:,1], ':g', lw = 2, label ="Runge-Kutta 4" )

plt.legend()
plt.xlabel( 'Coordenada x' )
plt.ylabel( 'Coordenada y' )
plt.title( r'Órbita de Kepler ($\Delta$t={})'.format(round(deltat,2)) )
plt.grid()
plt.show()
