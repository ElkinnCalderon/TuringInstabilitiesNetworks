# Turing Instabilities on Network
This code solves a reaction-diffusion system, where diffusion is defined over a complex network using the discrete Laplacian operator. However, the main interest lies in the case where the networks on which the species concentrations are defined are different.

The system has the following form:

(U, V)_t = (du * Lu * U, dv * Lv * V) + eta * (f(U, V; alpha_0, alpha_1, alpha_2, ...), g(U, V; alpha_0, alpha_1, alpha_2, ...))

In the file main.py, the parameters of the system are configured. The configuration of the Laplacians Lu and Lv is defined in the class Network.

Although the intention is for the complete configuration of the system to be contained within the class ProblemSetup.py, some parameters are kept in main.py to facilitate the exploration of parameter regions. For this purpose, the variable is defined as:

parametrosRD = [[du, dv], eta, alpha_0, alpha_1, alpha_2, ...]

The first component, parametrosRD[0], is an array with two elements containing the diffusion coefficients:
parametrosRD[0][0] corresponds to the diffusion coefficient of species u,
parametrosRD[0][1] corresponds to that of species v.

Although it is possible to directly configure these diffusion coefficients over the network in the module ProblemSetup.py, modifying them in main.py is also allowed to facilitate numerical exploration.

# Class
## Network
## ProblemSetup
## TuringAnalisis
## TemporalIntegration

# Code Availability  

The code used in this work is available at the GitHub repository: [https://github.com/usuario/repositorio](https://github.com/usuario/repositorio).  
The project is managed with version control in Git (via GitHub) and includes tagged releases to facilitate reproducibility.  

The main dependencies are:  

```bash
Python 3.11
numpy==2.0.1
matplotlib==3.9.2
networkx==3.2.1
