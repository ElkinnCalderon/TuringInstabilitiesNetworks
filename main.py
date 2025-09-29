"""
main.py

This script defines and runs the core components of the system.

Classes:
    - Network
    - ProblemSetup
    - TuringAnalisis
    - TemporalIntegration

Author: Elkinn Adrian Calderon Barreto
Date: 2025-06-02
"""
from TuringAnalisis      import Class_TuringAnalisis
from Network             import Class_NetworkTuring
from ProblemSetup        import Class_SistemaReaccionDifusion
from TemporalIntegration import Class_TemporalIntegration
import numpy as np
import os
import matplotlib.pyplot as plt
import networkx as nx

# Network parameters
n1, p1 = 10, 0.3  # Size and probability for the first network
n2, p2 = 0, 0.3   # Size and probability for the second network (smaller)
seed = 42         # Seed for reproducibility

GraphsUV=Class_NetworkTuring(N=n1,p=p1, M=n2,p1=p2,seed=seed)
GraphsUV.dibujar_red()

# Configuration and parameter settings of the reaction-diffusion system
eps=0.01
dc=0.0532762548739727
parametrosRD=[[dc-eps,1.0],1.0,-1.5,-1.25,0.5,1.0,0.1,1.5]
Sistema=Class_SistemaReaccionDifusion(parametros=parametrosRD,Lu=GraphsUV.Lu,Lv=GraphsUV.Lv,N=n1,M=n2)

# Turing instability analysis for the reaction-diffusion system (Sistema)

Turing=Class_TuringAnalisis(Sistema)
pcritico=Turing.buscar_difusion_critica_bisec()
print('valore critico',pcritico)
# Simulation of the reaction-diffusion system (Sistema)
#Simulation=Class_TemporalIntegration(Sistema,Nit=10000,Niter=1000,integrator='Euler',imprimir=False)


# numSimul=20
# for i in range(numSimul):
#     # Configuration and parameter settings of the reaction-diffusion system
#     eps=0.0025
#     dc=0.0532762548739727

#     a = 0.1
#     b = 1.5
#     du = dc-eps*i
#     dv = 1.0
#     eta = 1.0

#     parametrosRD=[[du,dv],eta,-1.5,-1.25,0.5,1.0,a,b]
#     Sistema=Class_SistemaReaccionDifusion(parametros=parametrosRD,Lu=GraphsUV.Lu,Lv=GraphsUV.Lv,N=n1,M=n2)
#     Simulation=Class_TemporalIntegration(Sistema,Nit=10000,Niter=5000,integrator='Euler',imprimir=False)

#     folder          = 'resultados/Patrones'
#     os.makedirs(folder, exist_ok=True)
#     filename = os.path.join(folder, f'CN_{i}.txt')
    
#     contenido = f"""a:{a} b:{b} du:{du} dv:{dv} eta:{eta} LogERROR:{np.log10(Simulation.difCNC0[-1])} 
# """
#     contenido += ' '.join(map(str, Simulation.CN))
#     with open(filename, 'w') as archivo:
#         archivo.write(contenido)





""" 
folder = "resultados/Patrones"

# Archivos ordenados
archivos = sorted(
    [f for f in os.listdir(folder) if f.startswith("CN_") and f.endswith(".txt")],
    key=lambda x: int(x.split("_")[1].split(".")[0])
)

# Leer datos CN
simulaciones_CN = []
for archivo in archivos:
    with open(os.path.join(folder, archivo), "r") as f:
        lineas = f.readlines()
        datos_CN = []
        for linea in lineas:
            if linea.startswith("a:"):
                continue
            datos_CN.extend(map(float, linea.strip().split()))
        simulaciones_CN.append(np.array(datos_CN))

# Dibujar cada simulación
for i, CN in enumerate(simulaciones_CN):
    mitad = len(CN) // 2
    CN_u = CN[:mitad]
    CN_v = CN[mitad:]

    fig, axs = plt.subplots(2, 1, figsize=(6, 6))

    # Grafo de concentración u
    axs[0].set_title(f"Simulación {i} - Red de la concentración $u$")
    nodes_u = nx.draw_networkx_nodes(
        GraphsUV.G_u, pos=GraphsUV.pos_u, ax=axs[0],
        node_size=150, node_color=CN_u, cmap='viridis'
    )
    plt.colorbar(nodes_u, ax=axs[0])

    # Grafo de concentración v
    axs[1].set_title(f"Simulación {i} - Red de la concentración $v$")
    nodes_v = nx.draw_networkx_nodes(
        GraphsUV.G_v, pos=GraphsUV.pos_v, ax=axs[1],
        node_size=150, node_color=CN_v, cmap='viridis'
    )
    plt.colorbar(nodes_v, ax=axs[1])

    plt.tight_layout()
    plt.show() """