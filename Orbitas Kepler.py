import matplotlib.pyplot as plt
from numpy import zeros

# Definición de la función F
def F(x, y, xd, yd):
    F1 = xd
    F2 = yd
    F3 = -x / (x**2 + y**2)**(3/2)
    F4 = -y / (x**2 + y**2)**(3/2)
    return F1, F2, F3, F4

# Parámetros del sistema
T = 1      # Tiempo total
N = 10   # Número de pasos (más grande para mayor precisión)
deltat = T / N

# Listas para almacenar los valores de x e y
x = zeros(N+1) # Inicializa un array de tamaño N+1
y = zeros(N+1)

# Condiciones iniciales
x[0] = 1
y[0] = 0
xd = 0
yd = 1

# Bucle para la integración numérica
for i in range(N):
    F1, F2, F3, F4 = F(x[i], y[i], xd, yd)
    x[i+1] = x[i] + deltat * F1
    y[i+1] = y[i] + deltat * F2
    xd = xd + deltat * F3
    yd = yd + deltat * F4
    

# Representación gráfica
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trayectoria en 2D')
plt.grid(True)
plt.axis('equal')  # Escala igual en ambos ejes
plt.show()