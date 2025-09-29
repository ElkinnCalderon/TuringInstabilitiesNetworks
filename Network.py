import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random

class Class_NetworkTuring:
    def __init__(self, red='', N=5,p=0.4, M=5,p1=0.5, seed=29):
        self.seed = seed
        self.N = N
        self.M = M

        if red != '':
            self.G, self.G1, self.G2 = self.unir_redes_erdos(N, p, M, p1, seed)
        else:
            self.G, self.G1, self.G2 = self.unir_redes_erdos(N, p, M, p1, seed)
            self.pos_u, self.pos_v=self.calcular_pos()

        self.G_u=self.G1
        self.G_v=self.G

        self.Lu=nx.laplacian_matrix(self.G_u).todense()
        self.Lv=nx.laplacian_matrix(self.G_v).todense()

        for node, (x, y) in self.pos_u.items():
            self.G_u.nodes[node]["x"] = float(x)
            self.G_u.nodes[node]["y"] = float(y)

        nx.write_gexf(self.G_u, "red_con_posiciones.gexf")


    def generar_red_conexa(self, n, p, seed=None):
        """Genera una red Erdős-Rényi conexa."""
        if seed is not None:
            random.seed(seed)
        while True:
            G = nx.erdos_renyi_graph(n, p, seed=random.randint(0, 1000000) if seed is not None else None)
            if nx.is_connected(G):
                return G

    def unir_redes_erdos(self, n1, p1, n2, p2, seed=None):
        """Une dos redes Erdős-Rényi (la segunda puede ser vacía)."""
        seed1 = seed if seed is None else seed + 1
        seed2 = seed if seed is None else seed + 2

        G1 = self.generar_red_conexa(n1, p1, seed1)
        if n2 > 0:
            G2 = self.generar_red_conexa(n2, p2, seed2)
            mapping = {node: node + n1 for node in G2.nodes()}
            G2 = nx.relabel_nodes(G2, mapping)
            G = nx.compose(G1, G2)
            nodo1 = list(G1.nodes())[1]
            nodo2 = list(G2.nodes())[0]
            G.add_edge(nodo1, nodo2)
        else:
            G2 = G1
            G = G1

        return G, G1, G2
    
    def calcular_pos(self):
        posu = nx.circular_layout(self.G1)
        for key in posu:
            posu[key][0] -= 1.5  # Separar la primera red hacia la izquierda
        
        # Dibujar la segunda red con disposición circular pero desplazada y rotada
        posv = nx.circular_layout(self.G2)
        theta = 0.0*np.pi  # Rotar la segunda red 45 grados
        rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        for key in posv:
            rotated_pos = np.dot(rotation_matrix, np.array([posv[key][0], posv[key][1]]))
            posv[key] = [rotated_pos[0] + 1.5, rotated_pos[1] * 0.7]
        posv=posu | posv
        return posu,posv

    def dibujar_red(self):
        fig, axs = plt.subplots(2, 1, figsize=(6,6))


        color_map_u = ['skyblue'] *self.N
        color_map_v = ['skyblue'] *self.N + ['lightcoral'] * self.M

        # Graficar la red original en forma georeferenciada
        axs[0].set_title("Red de la concentración $v$")
        
        nx.draw(self.G_u, pos=self.pos_u, ax=axs[0], with_labels=True, node_size=700, node_color=color_map_u, font_size=12, font_weight='bold', edge_color='gray')
        # Graficar la red extendida en forma georeferenciada
        axs[1].set_title("Red de la concentración $p$" )
        nx.draw(self.G_v, pos=self.pos_v, ax=axs[1], with_labels=True, node_size=700, node_color=color_map_v, font_size=12, font_weight='bold', edge_color='gray')

        plt.show()