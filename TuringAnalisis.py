import numpy as np
from ProblemSetup import Class_SistemaReaccionDifusion
from scipy.linalg import eigvals
import math as mt
import matplotlib.pyplot as plt

class Class_TuringAnalisis:
    def __init__(self,SisReacDif):
        self.SisReacDif          =SisReacDif
        self.J                   =self.Jacobiano(SisReacDif.E)
        self.ValP_J,self.VecP_J  =np.linalg.eig(self.J)
        
        self.ValP_L,self.VecP_L  =np.linalg.eig(SisReacDif.L)
        print(self.VecP_L)

        self.Dcrit               =self.difusionCritica()
        self.RelacionDispersion  =self.relacion_dispersion_discreta()
        self.RelacionH           =self.analisis_H()
        self.omega = None
        self.delta = None
        self.mostrar_resultados()
        self.graficaRelacionDispersion()

    def Jacobiano(self,W):
        N=self.SisReacDif.N
        M=self.SisReacDif.M
        J=np.zeros((N + M, N + M))
        h=0.00001
        H = np.zeros(N + M)
        t0=0.0
        for i in range(N+M):
            for j in range(N+M):
                H = 0*H
                H[j]=h
                J[i][j]=(self.SisReacDif.Reaccion(W+H,i,t=t0)[0]-self.SisReacDif.Reaccion(W-H,i,t=t0)[0])/(2.0*h)
        return J

    def difusionCritica(self):
        N=self.SisReacDif.N
        M=self.SisReacDif.M
        difusion=[]
        for i in range(N):
            fu=self.J[i][i]
            fv=self.J[i][N+i]
            gu=self.J[N+i][i]
            gv=self.J[N+i][N+i]
            a=gv**2
            b=2*(2*fv*gu-fu*gv)
            c=fu**2
            difusion.append((-b-mt.sqrt(b**2-4*a*c))/(2*a))
        return difusion

    def MatrizH(self,VecInd):
        N=self.SisReacDif.N
        M=self.SisReacDif.M
        H=np.zeros((len(VecInd),len(VecInd)))
        for i in range(len(VecInd)):     
            MLi = -self.ValP_L[i]* np.identity(N + M)
            Msum = np.dot(self.SisReacDif.D, MLi) + self.J
            for j in range(len(VecInd)):
                H[j][i]=np.dot(self.VecP_L[:,VecInd[j]],np.dot(Msum,self.VecP_L[:,VecInd[i]]))
        return H

    def calcular_omega(self, p):
        if self.A is None or self.valores_propios is None:
            raise ValueError("A o los valores propios no están definidos.")
        self.omega = []
        for mu in self.valores_propios:
            M = self.A - mu * p
            valores = eigvals(M)
            parte_real = max(valores.real)
            self.omega.append(parte_real)
        return np.array(self.omega)

    def buscar_difusion_critica_bisec(self,max_iter=100, tolerancia=1e-6):
        N=self.SisReacDif.N
        M=self.SisReacDif.M
        p=self.Dcrit[0] 
        # Valores iniciales para la bisección (ajusta estos valores según sea necesario)
        p_min = 0.000001        # Límite inferior del intervalo de búsqueda
        p_max = 0.9      # Límite superior del intervalo de búsqueda
        iteraciones = 0

        while iteraciones < max_iter:
            # Calcular el punto medio para el valor de p[0][0]
            p = (p_min + p_max) / 2.0

            # Recalcular D, H y omega con el valor actual de p[0][0]
            self.SisReacDif.p[0][0]=p
            self.SisReacDif.D=self.SisReacDif.construir_matriz_difusion()
            rango_indices=range(self.SisReacDif.N+self.SisReacDif.M)
            H = self.MatrizH(rango_indices)
            valores_propiosH, vectores_propiosH = np.linalg.eig(H)
            omega = max(valores_propiosH.real)

            print('búsqueda de d', iteraciones, p, omega)

            # Condición de parada: si omega es cercano a cero
            if abs(omega) < tolerancia:
                print(f"Valor encontrado: p[0][0] = {p} después de {iteraciones} iteraciones")
                break

            # Actualizar intervalo de búsqueda
            if omega < 0:
                p_max = p
            else:
                p_min = p

            iteraciones += 1

            self.delta = p
        return self.delta
    
    def relacion_dispersion_discreta(self):
        """
        Calcula el máximo valor propio real de Msum para cada valor propio discreto.
        """
        N=self.SisReacDif.N
        M=self.SisReacDif.M
        detMsum = np.zeros(len(self.ValP_L))
        print(self.ValP_L)
        for i, Vali in enumerate(self.ValP_L):
            MLi = -Vali * np.identity(N + M)
            Msum = np.dot(self.SisReacDif.D, MLi) + self.J
            valores_propiosMsum, vectores_propiosMsum = np.linalg.eig(Msum)
            detMsum[i] = max(valores_propiosMsum.real)
        
        return detMsum


    def relacion_dispersion_densa(self, delta=0.01, margena=0.2, margenb=5):
        """
        Evalúa la relación de dispersión sobre una malla densa de valores.
        """
        N=self.SisReacDif.N
        M=self.SisReacDif.M
        print(self.ValP_L) 
        ValiDenso = np.arange(min(self.ValP_L) - margena, max(self.ValP_L)+margenb, delta)
        detMsum = np.zeros(len(ValiDenso))
        
        for i, Vali in enumerate(ValiDenso):
            MLi = -Vali * np.identity(N + M)
            Msum = np.dot(self.SisReacDif.D, MLi) + self.J
            valores_propiosMsum, vectores_propiosMsum = np.linalg.eig(Msum)
            detMsum[i] = max(valores_propiosMsum.real)
        
        return ValiDenso, detMsum

    def analisis_H(self):
        rango_indices=range(self.SisReacDif.N+self.SisReacDif.M)
        H = self.MatrizH(rango_indices)
        # Calcular autovalores
        valores_propios_H, vectores_propios_H = np.linalg.eig(H)
        omega = max(valores_propios_H.real)
              
        return omega, valores_propios_H

    def mostrar_resultados(self):
        print("Difusion Critica Caso Usual:\n", self.Dcrit)
        print("Relación de dispersion:\n", self.RelacionDispersion)
        print("Omega:\n", self.RelacionH[0])

    def graficaRelacionDispersion(self):
        x, y=self.relacion_dispersion_densa()

        ydiscre=self.RelacionDispersion
        xdiscre=self.ValP_L

        fig, axs = plt.subplots(1, 1, figsize=(6,6))


        axs.set_title("Relación de dispersión")
        axs.scatter(xdiscre, ydiscre, color='g', marker='o')

        axs.plot(x, y, color='g')
        axs.plot(x, np.zeros(len(y)), color='b')        
        axs.plot(x, np.zeros(len(y))+self.RelacionH[0], color='r')

        plt.show()