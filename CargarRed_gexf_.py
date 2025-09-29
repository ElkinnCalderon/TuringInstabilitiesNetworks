import networkx as nx
import matplotlib.pyplot as plt

# Paso 1: Cargar la red desde archivo .gexf
G = nx.read_gexf("red_con_posiciones.gexf")

# Paso 2: Extraer posiciones desde atributos 'x' y 'y'
pos = {node: (float(data['x']), float(data['y'])) for node, data in G.nodes(data=True)}

# Paso 3: Dibujar la red usando las posiciones
plt.figure(figsize=(8, 6))
nx.draw(G, pos, with_labels=True, node_size=300, node_color='lightblue', edge_color='gray')
plt.title("Red cargada con posiciones desde GEXF")
plt.axis('equal')
plt.show()