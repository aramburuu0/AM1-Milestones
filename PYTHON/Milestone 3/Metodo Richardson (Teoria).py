from numpy import zeros, linspace

def Cauchy(Esquema, U0, F, t):

    U = zeros([len(t), len(U0)])
    U[0, :] = U0

    for n in range(0, len(t)-1):
        
        U[n+1, :] = Esquema(U[n, :], t[n+1] - t[n], t[n], F)
    
    return U

def Cauchy_error(Metodo, Esquema, U0, F, t, q):

    t1 = t
    t2 = linspace(t(0), t(-1), 2*len(t1))   # el N de linspace es el número de nodos 
    
    U1 = Metodo(Esquema, U0, F, t1)
    U2 = Metodo(Esquema, U0, F, t2)
    Error = zeros([len(t1), len(U0)])

    for n in range(0, len(t1)-1):

        Error[n, :] = (U2[2*n, :]-U1[n, :])/(1-1/2**q)

    return U1, Error

def Refinamiento(t1):

    N = len(t1) - 1
    t2 = zeros(2*len(t1))

    for n in range (0, N):

        t2[2*n] = t1[n]                     # nodos pares
        t2[2*n+1] = (t1[n+1]+t1[n])/2       # nodos impares

    t2[2*N] = t1[n]

    return t2

# Hace una partición equispaciada en N trozos de un segmento de ka recta real de a, b basicamente en linspace

def Particion(a, b, N):

    t = zeros(N+1)

    for n in range(0, N):

        t[n] = a + (b-a)/N * N

    return t