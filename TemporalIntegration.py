import numpy as np
import os
from ProblemSetup import Class_SistemaReaccionDifusion


class Class_TemporalIntegration:
    def __init__(self,SisReacDif,Nit=100,Niter=1,integrator='Euler'):
        self.integrator          =integrator
        self.SisReacDif          =SisReacDif
        self.Nit                 =Nit
        self.Niter               =Niter
        self.epsilon             =0.05
        self.delta               =0.0001
        self.CN,self.difCNC0     =self.Simulate()


    def MetEulerReaccion(self,Nit,C0,D,L,t,delta=0.0001):
        C=C0
        R=np.zeros(len(C0))
        for i in range(Nit):
            R=self.SisReacDif.vector_reaccion(C,t)
            C=delta*(np.dot(np.dot(D,L),C)+R)+C
            t=t+delta
        return C,t


    def RHS(self,C, D, L,t):
        R = self.SisReacDif.vector_reaccion(C,t)
        return np.dot(np.dot(D, L), C) + R

    def MetRK5Reaccion(self,Nit, C0, D, L, t, delta=0.0001):
        C = C0.copy()
        for _ in range(Nit):
            k1 = delta * self.RHS(C, D, L,t)
            k2 = delta * self.RHS(C + k1/4, D, L,t)
            k3 = delta * self.RHS(C + (3/32)*k1 + (9/32)*k2, D, L,t)
            k4 = delta * self.RHS(C + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3, D, L,t)
            k5 = delta * self.RHS(C + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4, D, L,t)
            k6 = delta * self.RHS(C - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5, D, L,t)

            C = C + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6
            t = t + delta

        return C, t

    def Simulate(self):
        N           =self.SisReacDif.N
        M           =self.SisReacDif.M
        D           =self.SisReacDif.D
        L           =self.SisReacDif.L
        Nit         =self.Nit
        Niter       =self.Niter
        
        C0          =np.zeros(N+M)
        CN          =np.zeros(N+M)

        C0          =self.SisReacDif.E
        
        print('Equilibrium')
        print(C0)

        epsilon      =self.epsilon 
        delta0       =self.delta

        for i in range(N+M):
            num_aleatorio = np.random.uniform(-epsilon, epsilon)
            C0[i]=C0[i]+num_aleatorio
            
        print('ITER',Nit)
        print('C0',C0)
        difCNC0=[]
        t=0.0
                
        folder          = 'resultados/Simulate'
        os.makedirs(folder, exist_ok=True)
        filename        = os.path.join(folder, 'vector_iteracion_0.txt')

        with open(filename, 'w') as archivo:
            contenido = ' '.join(map(str, C0))
            archivo.write(contenido)
        
        for i in range(1,Niter+1):
            print('ITER',Nit,'salto',i,'t',t)
            if self.integrator=='Euler':
                CN,t=self.MetEulerReaccion(Nit,N,C0,D,L,t,delta0)
            elif self.integrator=='RK5':
                CN,t  =self.MetRK5Reaccion(Nit=Nit, C0=C0, D=D, L=L, t=t, delta=delta0)

            difCNC0.append(np.linalg.norm(CN-C0))
            
            print('CN',CN)
            print('norm(CN-C0)',np.linalg.norm(CN-C0))
            
            C0=CN
            filename = os.path.join(folder, f'vector_iteracion_{i}.txt')
            with open(filename, 'w') as archivo:
                contenido = ' '.join(map(str, CN))
                archivo.write(contenido)

        
        return CN,difCNC0
    