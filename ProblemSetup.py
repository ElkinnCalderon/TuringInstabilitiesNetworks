import numpy as np
from scipy.linalg import block_diag

class Class_SistemaReaccionDifusion:
    def __init__(self, parametros,Lu,Lv,N,M):
        """
        parámetros esperados:
        parametros[0] = [D_u, D_v] → coeficientes de difusión
        parametros[1] = eta        → constante de escala temporal
        parametros[2:8] = a, b, c, h, a0, b0 → parámetros del sistema
        """
        self.N    =N
        self.M    =N + M
        self.p    = parametros
        self.Lu   =Lu
        self.Lv   =Lv
        self.L    =block_diag(Lu,Lv)
        self.OpDif=self.OperadorDifussion()
        self.D    =self.construir_matriz_difusion()
        self.E    =self.Equilibrios()
    
    def Reaccion(self, W, i,t=0.0):

        a, b, c, h, a0, b0 = self.p[2:8]
        eta = self.p[1]
        N=self.N
        M=self.M

        div=1

        if i < N:
            u = W[i]
            v = W[i + N]
            if i < int(N / div):
                f = eta * (a0 - u + u**2 * v)
                g = eta * (b0 - u**2 * v)
            else:
                f = eta * (u + a*v - c*u*v - v**2 * u)
                g = eta * (h*u + b*v + c*u*v + v**2 * u)
            FE = f
        elif N <= i < 2 * N:
            u = W[i - N]
            v = W[i]
            if i < N + int(N / div):
                f = eta * (a0 - u + u**2 * v)
                g = eta * (b0 - u**2 * v)
            else:
                f = eta * (u + a*v - c*u*v - v**2 * u)
                g = eta * (h*u + b*v + c*u*v + v**2 * u)
            FE = g
        else:
            v = W[i]
            g = -0.1 * eta * (v - (b0 / (a0 + b0)**2))
            FE = g
            f = 0

        return FE, f, g

    def construir_matriz_difusion(self):
        N=self.N
        M=self.M
        D = np.zeros(N + M)
        D[:N] = self.p[0][0]
        D[N:N+N] = self.p[0][1]
        D[2*N:] = self.p[0][1]
        return np.diag(D)
    
    def vector_reaccion(self, W,t=0.0):
        N=self.N
        M=self.M
        F = np.zeros(N + M)
        for i in range(N + M):
            F[i] = self.Reaccion(W, i,t)[0]
        return F
    
    def OperadorDifussion(self):
        N=self.N
        M=self.M
        OpLineal = np.zeros((N + M,N + M))
        OpLineal = self.L
        return OpLineal 
    
    def Equilibrios(self):
        N=self.N
        M=self.M
        a, b, c, h, a0, b0 = self.p[2:8]    
        div=1
        F = np.zeros(N + M)
        for i in range(N+M):    
            if i< N:
                
                if i<int(N/div):
                    F[i]=a0+b0
                else:
                    F[i]=0.0
            elif i> N-1 and i< 2*N:      
                
                if i<N+int(N/div):
                    F[i]=b0/(a0+b0)**2.0
                else:
                    F[i]=0.0
            else:
                F[i]=b0/(a0+b0)**2.0
        
        return F