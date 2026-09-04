# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 13:01:32 2026

@author: Trabajador
"""

from flatfold.exact import (pendiente_admisible, 
                            reflexion, 
                            compon_isometrias, 
                            es_identidad, direccion)
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
            
    




def sectores(G, vecinos, v):
    """
    Input:  el grafo, el diccionario de vecinos ordenados, un vértice.
    Output: lista de (angulo, w_anterior, w_siguiente) en orden cíclico,
            con el ángulo en unidades de 45 grados (entero de 1 a 8).

    Los pliegues incidentes en v dividen el plano en sectores. Este es el
    dato del que sale la condición big-little-big: los dos pliegues que
    delimitan un sector que es mínimo local estricto deben tener etiquetas
    opuestas.

    Las aristas planas NO delimitan sector: el papel sigue recto a través
    de ellas y no hay doblez.

    El ángulo se calcula sobre las ocho direcciones admisibles, así que es
    exacto y entero: no hace falta trigonometría.
    """
    pliegues = [w for w in vecinos[v] if G[v][w]["tipo"] == Kind.CREASE]
    n = len(pliegues)
    if n < 2:
        return []

    dirs = [direccion(G, v, w) for w in pliegues]

    salida = []
    for i in range(n):
        w1, w2 = pliegues[i], pliegues[(i + 1) % n]
        d1, d2 = dirs[i], dirs[(i + 1) % n]
        ang = (d2 - d1) % 8
        if ang == 0:
            ang = 8              # una vuelta completa: vértice de grado 1
        salida.append((ang, w1, w2))
    return salida




def minimos_locales(secs):
    """
    Input:  la lista que devuelve sectores(): [(angulo, w1, w2), ...]
            en orden cíclico.
    Output: los índices de los sectores cuyo ángulo es estrictamente
            menor que el de sus dos vecinos.

    Es la condición del big-little-big: un sector "pequeño" entre dos
    "grandes". Los dos pliegues que lo delimitan deben tener etiquetas
    opuestas, porque si fueran iguales el papel se atravesaría.

    Con menos de tres sectores no hay nada que comparar. Y si todos los
    ángulos son iguales —el caso de las bombas de agua— no hay ningún
    mínimo local.
    """
    n = len(secs)
    if n < 3:
        return []
    return [i for i in range(n)
            if secs[i][0] < secs[i - 1][0]
            and secs[i][0] < secs[(i + 1) % n][0]]






# def main():
    
#     from pathlib import Path
#     from flatfold.io_fold import read_fold
#     from flatfold.pattern import crea_grafo_dirigido, vecinos_ordenados


#     TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

#     bomb = TOY / "waterbomb4.fold"

#     r = read_fold(bomb)
#     G = crea_grafo_dirigido(r)
    
#     # vecinos = vecinos_ordenados(G)
    
#     # for i in G.edges:
#     #     print(i, G.edges[i])
#     valida_kawasaki(G)
    

# if __name__ == "__main__":
#     main()
     








