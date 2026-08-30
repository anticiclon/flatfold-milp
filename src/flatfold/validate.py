# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 13:01:32 2026

@author: Trabajador
"""

from flatfold.exact import pendiente_admisible, reflexion, compon_isometrias, es_identidad
from flatfold.io_fold import Kind
from flatfold.pattern import vecinos_ordenados
# import numpy as np


def valida_pendientes(fold_input):
    for idx, (u, v) in enumerate(fold_input.aristas):
        p = fold_input.vertices[u]
        q = fold_input.vertices[v]
        if not pendiente_admisible(p, q):
            raise ValueError(
                f"arista {idx} ({u}->{v}, {p}->{q}): pendiente no admisible")


def vertices_interiores(G):
    """
    Input: el grafo.
    Output: el conjunto de índices de vértices que están en el interior del 
    papel.
    Un vértice es interior si ninguna de sus aristas incidentes es de
    frontera. Recorres los vértices y, para cada uno, miras el tipo de las 
    aristas que salen de él.
    """
    interiores = set()
    for v in G.nodes:
        tipos = [G[v][w]["tipo"] for w in G.neighbors(v)]
        if all(t != Kind.BOUNDARY for t in tipos):
            interiores.add(v)
    return interiores


def holonomia(G, vecinos, v):
    """
    Input: el grafo, el diccionario de vecinos ordenados, y un vértice.
    
    Output: una isometría, la composición acumulada.
    
    Recorre los vecinos de v en orden angular, se queda solo con aquellos
     unidos por un pliegue, y compone las reflexiones respecto de esas rectas. 
     Empieza desde la identidad y va acumulando.
    """
    acumulada = ((1, 0, 0, 1), (0, 0))
    for w in vecinos[v]:
        if G[v][w]["tipo"] == Kind.CREASE:
            r = reflexion(G.nodes[v]["coor"], G.nodes[w]["coor"])
            acumulada = compon_isometrias(acumulada, r)
    return acumulada



def valida_kawasaki(G):
    """
    Input: el grafo. (Internamente llama a las dos anteriores.)
    Output: nada, o ValueError diciendo qué vértice falla.
    Para cada vértice interior, calcula su holonomía y comprueba que es la identidad.
    """
    vecinos = vecinos_ordenados(G)
    for v in vertices_interiores(G):
        h = holonomia(G, vecinos, v)
        if not es_identidad(h):
            raise ValueError(f"Kawasaki falla en el vértice {v}: holonomía = {h}")
        n = sum(1 for w in vecinos[v] if G[v][w]["tipo"] == Kind.CREASE)
        if n % 2 != 0:
            raise ValueError(f"vértice {v}: grado impar ({n} pliegues)")
            
    




def main():
    
    from pathlib import Path
    from flatfold.io_fold import read_fold
    from flatfold.pattern import crea_grafo_dirigido, vecinos_ordenados


    TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

    bomb = TOY / "waterbomb4.fold"

    r = read_fold(bomb)
    G = crea_grafo_dirigido(r)
    
    # vecinos = vecinos_ordenados(G)
    
    # for i in G.edges:
    #     print(i, G.edges[i])
    valida_kawasaki(G)
    

if __name__ == "__main__":
    main()
     








