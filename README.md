# Turing Instabilities on Network

This repository contains code to solve a reaction–diffusion system defined on complex networks, where diffusion is modeled using the discrete Laplacian operator.
The main focus is the case where the networks on which the species concentrations are defined are **different**.

---

## Model Description

The system is written as:

[
(U, V)_t = (d_u L_u U, ; d_v L_v V) + \eta , \big(f(U, V; \alpha_0, \alpha_1, \alpha_2, \dots), ; g(U, V; \alpha_0, \alpha_1, \alpha_2, \dots)\big)
]

* `main.py` configures the system parameters.
* The Laplacians (L_u) and (L_v) are defined in the **Network** class.

Although the intention is to centralize configuration in `ProblemSetup.py`, some parameters remain in `main.py` to facilitate exploration of parameter regions. For this purpose, the variable is defined as:

```python
parametrosRD = [[du, dv], eta, alpha_0, alpha_1, alpha_2, ...]
```

* `parametrosRD[0][0]` → diffusion coefficient of species **u**
* `parametrosRD[0][1]` → diffusion coefficient of species **v**

While it is possible to directly configure these coefficients in `ProblemSetup.py`, modifying them in `main.py` is also allowed for convenience during numerical experiments.

---

## Code Structure

* **Network** → Handles network construction and Laplacians.
* **ProblemSetup** → Defines system parameters and configuration.
* **TuringAnalysis** → Performs linear stability and Turing instability analysis.
* **TemporalIntegration** → Integrates the system over time.

---

## Code Availability

The code used in this work is available at the GitHub repository:
👉 [https://github.com/ElkinnCalderon/TuringInstabilitiesNetworks](https://github.com/ElkinnCalderon/TuringInstabilitiesNetworks)

The project is managed with Git (via GitHub) and includes tagged releases to facilitate reproducibility.

---

## Requirements

```bash
Python 3.11
numpy==2.0.1
matplotlib==3.9.2
networkx==3.2.1
```

---

## License

This project is distributed under the **MIT License**, which allows use, modification, and redistribution with proper attribution.

---

## Contact

📧 [elkinncb@gmail.com](mailto:elkinncb@gmail.com)
📧 [jlaragon@unam.mx](mailto:jlaragon@unam.mx)

