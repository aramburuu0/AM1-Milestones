from ODES.Funciones import RK4, Cauchy
from numpy import zeros, reshape, linspace, ceil, sqrt, pi
from numpy.linalg import norm
from random import random, randint, uniform
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def Condiciones(Nb, Nc):
    
  '''''''''''
  INPUTS:
      - Nb: número de cuerpos
      - Nc: número de coordenadas
  OUTPUTS:
      - U0: un vector largo que contiene todas las coordenadas generadas de forma aleatoria en un orden lineal

  '''''''''''
  U0 = zeros(2*Nc*Nb) # vector de dimensiones 2*Nb*Nc (dos coordenadas por cuerpo)
  U1 = reshape(U0, (Nb, Nc, 2)) # Matriz NbxNc y en cada coordenada un "vector" de 2 coordenadas # REORDENA EL VECTOR U0 EN UNA MATRIZ 3D DE Nb*Nc*2 (Nº de filas, Nº de columnas, Nº de "capas")
  r0 = reshape(U1[:, :, 0], (Nb, Nc))
  v0 = reshape(U1[:, :, 1], (Nb, Nc))
  

  for i in range(Nb):
      r0[i, :] = [uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)]
      v0[i, :] = [uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)]
  
  m = zeros(Nb)
  for i in range(Nb):
      m[i] = uniform(0.001, 1)

  return U0, m



def F_NBody(U, t, Nb, Nc, m): 
     
#   Write equations: Solution( body, coordinate, position-velocity )      
  Us  = reshape( U, (Nb, Nc, 2) )  # ES COMO U1
  F =  zeros(len(U))   # INICIALIZA EL VECTOR DE FUERZAS A 0 CON LA LONGITUD DE U (INPUT QUE SERÁ U0)
  dUs = reshape( F, (Nb, Nc, 2) )  # DE NUEVO EL RESHAPE DE LA MISMA FORMA. ESTA MATRIZ 3D GUARDARÁ LOS CAMBIOS EN POSICIÓN Y VELOCIDAD POR LA FUERZA
  
  r = reshape( Us[:, :, 0], (Nb, Nc) )     # SACA LA POSICIÓN DEL INPUT Y LO PONE COMO MATRIZ 2D
  v = reshape( Us[:, :, 1], (Nb, Nc) )     # SACA LA VELOCIDAD DEL INPUT Y LO PONE COMO MATRIZ 2D
  
  drdt = reshape( dUs[:, :, 0], (Nb, Nc) ) # SACA LA DERIVADA DE LA POSICIÓN DEL dUs Y LO PONE COMO MATRIZ 2D
  dvdt = reshape( dUs[:, :, 1], (Nb, Nc) ) # SACA LA DERIVADA DE LA VELOCIDAD DEL dUs Y LO PONE COMO MATRIZ 2D

  
  dvdt[:,:] = 0  # CONDICIÓN INICIAL DE ACELERACIÓN 0 (NADA LES AFECTA DE INICIO) 

  for i in range(Nb):   
    drdt[i,:] = v[i,:]  # la derivada de la posición es la velocidad
    for j in range(Nb): 
      if j != i:  # CALCULA LA FUERZA SOBRE LA PARTÍCULA i DE TODAS LAS DEMÁS PARTÍCULAS j
        d = r[j,:] - r[i,:] # r_i - r_j (distancia entre partículas)
        dvdt[i,:] = dvdt[i,:] + m[i] * d[:] / norm(d)**3 # F = F + d / |d|^3 (el G*m nos lo fumamos?)
  
  print('Masas:', m)

  return F


#------------------------------------------------------------------
# Orbits of N bodies 
#      U : state vector
#      r, v: position and velocity points to U     
#------------------------------------------------------------------    
def Integrate_NBP():  
    
  def F(U,t):
    
    return F_NBody(U, t, Nb, Nc, m) 

  N =  1000    # time steps 
  Nb = 5      # boodies 
  Nc = 3      # coordinates 

  t0 = 0; tf = 2 * pi 
  Time = linspace(t0, tf, N+1) # Time(0:N) 

  U0, m = Condiciones( Nb, Nc )

#  U = odeint(F_NBody, Nb, Time, m) # Esto se supone que resolvería el sistema de ecuaciones diferenciales definido por la fuerza F para obtener el vector de estado U en cada paso temporal
  U = Cauchy(RK4, U0, F, Time) 

  Us  = reshape( U, (N+1, Nb, Nc, 2) ) # Convierte U en un array 4D que contiene las posiciones y velocidades de cada uno de los cuerpos en cada uno de los pasos temporales
  r   = reshape( Us[:, :, :, 0], (N+1, Nb, Nc) ) # Saca la posición como un array 3D
  
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  for i in range(Nb):
    ax.plot3D(r[:, i, 0], r[:, i, 1], r[:, i, 2])

  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  ax.set_title('Órbitas 3D')
  ax.axis('equal')
  ax.grid(True)
  plt.show()

  for i in range(Nb):
    plt.plot(  r[:, i, 0], r[:, i, 1] )
  plt.axis('equal')
  plt.grid()
  plt.xlim(-1.5, 1.5)
  plt.ylim(-1.5, 1.5)
  plt.xlabel('X')
  plt.ylabel('Y')
  plt.show()

  for i in range(Nb):
    plt.plot(  r[:, i, 0], r[:, i, 2] )
  plt.axis('equal')
  plt.grid()
  plt.xlim(-1.5, 1.5)
  plt.ylim(-1.5, 1.5)
  plt.xlabel('X')
  plt.ylabel('Z')
  plt.show()

  for i in range(Nb):
    plt.plot(  r[:, i, 1], r[:, i, 2] )
  plt.axis('equal')
  plt.grid()
  plt.xlim(-1.5, 1.5)
  plt.ylim(-1.5, 1.5)
  plt.xlabel('Y')
  plt.ylabel('Z')
  plt.show()




Integrate_NBP()