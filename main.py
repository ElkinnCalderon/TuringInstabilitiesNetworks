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

# Network parameters
n1, p1 = 10, 0.3  # Size and probability for the first network
n2, p2 = 0, 0.3   # Size and probability for the second network (smaller)
seed = 42         # Seed for reproducibility

GraphsUV=Class_NetworkTuring(N=n1,p=p1, M=n2,p1=p2,seed=seed)
GraphsUV.dibujar_red()

# Configuration and parameter settings of the reaction-diffusion system
eps=0.01
dc=0.0532762548739727
p=[[dc-eps,1.0],1.0,-1.5,-1.25,0.5,1.0,0.1,1.5]
Sistema=Class_SistemaReaccionDifusion(parametros=p,Lu=GraphsUV.Lu,Lv=GraphsUV.Lv,N=n1,M=n2)

# Turing instability analysis for the reaction-diffusion system (Sistema)
Turing=Class_TuringAnalisis(Sistema)

# Simulation of the reaction-diffusion system (Sistema)
Simulation=Class_TemporalIntegration(Sistema,Nit=1000,Niter=5,integrator='RK5')