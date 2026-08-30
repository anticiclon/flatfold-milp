# -*- coding: utf-8 -*-

from flatfold.exact import direccion, reflexion, aplica_isometria, compon_isometrias
from flatfold.io_fold import Kind
from collections import deque
import networkx as nx




def crea_grafo_dirigido(fold_input):
    
    G = nx.DiGraph()
    for idx, value in enumerate(fold_input.vertices):
        G.add_node(idx, coor=value)    
        
    for (u, v), tipo in zip(fold_input.aristas, fold_input.tipos):
        G.add_edge(u, v, tipo=tipo)
        G.add_edge(v, u, tipo=tipo)

    return G



def vecinos_ordenados(G):
    """Por cada vértice, da los vecinos en orden antihorario."""
    return {v: sorted(G.neighbors(v), key=lambda w: direccion(G, v, w))
            for v in G.nodes}


def proxima_arista(vecinos, u, v):
    """Siguiente arista dirigida al bordear la cara desde u -> v."""
    vecinos_v = vecinos[v]
    indice_u = vecinos_v.index(u)          # la vuelta atrás
    return (v, vecinos_v[indice_u - 1])    # el índice -1 ya da la vuelta



def traza_cara(vecinos, arista_inicial):
    """Dado el diccionario de cada vértice y sus vecinos y la arista inicial,
    da la lista de aristas que limita la cara dictada por esa arista inicial"""
    cara = []
    h = arista_inicial
    while True:
        cara.append(h)
        h = proxima_arista(vecinos, *h)
        if h == arista_inicial:
            return cara



def caras_grafo(G):
    "Da todas las caras del grafo"
    vecinos = vecinos_ordenados(G)
    visitadas = set()
    caras = []
    for arista in G.edges():
        if arista in visitadas:
            continue
        cara = traza_cara(vecinos, arista)
        visitadas.update(cara)
        caras.append(cara)
    return caras

        

def area_cara(G, cara):
    "Cálcula el área de la cara."
    vertices = [G.nodes[i[0]]["coor"] for i in cara]
    desplazados = vertices[1:] + vertices[:1] 
    doble_area = sum([i[0]*j[1] - j[0]*i[1] 
                      for i, j in zip(vertices, desplazados)])
    return doble_area
    

def caras_definitivas(G, caras):
    caras_positivas = []
    areas_positivas = []
    cara_negativa = []
    area_negativa = []
    for cara in caras:
        doble_area = area_cara(G, cara)
        if doble_area > 0:
            caras_positivas.append(cara)
            areas_positivas.append(doble_area)
        elif doble_area < 0:
            cara_negativa.append(cara)
            area_negativa.append(doble_area)        
        elif doble_area == 0:
            raise ValueError("cara de área nula")
    if len(cara_negativa) != 1:
        raise ValueError(f"se esperaba una cara exterior, hay {len(cara_negativa)}")
    if sum(areas_positivas) != abs(sum(area_negativa)):
        raise ValueError("las áreas no cuadran")
    return caras_positivas, cara_negativa[0]



def aristas_caras(G, caras):
    cara_de = {}                     
    for idx, cara in enumerate(caras):
        for arista in cara:
            cara_de[arista] = idx

    dic_aristas = {}
    for (u, v) in G.edges:
        if u < v:                     # cada arista una sola vez
            dic_aristas[(u, v)] = (cara_de[(u, v)], cara_de[(v, u)])
    return dic_aristas


def grafo_dual(G, caras, idx_exterior):
    """
    Input: el grafo, la lista de caras y el índice de la exterior.

    Output: para cada cara interior, la lista de sus vecinas junto con la 
    arista que las separa. Algo de la forma {cara: [(cara_vecina, u, v), ...]},
    donde (u,v) identifica la arista compartida.
    
    Parte de aristas_caras, que ya te da qué dos caras separa cada 
    arista, y le da la vuelta: en lugar de indexar por arista, indexa por cara. 
    Descarta las aristas que tocan la exterior.
    """
    dic_aristas = aristas_caras(G, caras)
    
    H = nx.MultiGraph()
    
    for (u, v), (f, g) in dic_aristas.items():
        if f != idx_exterior and g != idx_exterior:
            H.add_edge(f, g, arista=(u, v))
    
    return H



# def plegar(G, dual, cara_referencia):
#     """
#     Input: el grafo, el dual, y la cara desde la que empezar.
#     Output: un diccionario {cara: isometría}.
#     La cara de referencia recibe la identidad. Al pasar de una cara f ya 
#     alcanzada a una vecina g cruzando la arista (u,v):

#     si G[u][v]["tipo"] == Kind.CREASE, la isometría de g es 
#     compon_isometrias(iso[f], reflexion(coor[u], coor[v]));
#     si es FLAT, la isometría de g es la misma que la de f
#     """
#     dicc_iso = {}
#     dicc_iso[cara_referencia] = ((1, 0, 0, 1), (0, 0))
#     for i in dual.nodes:
#         if i != cara_referencia:
#             dual_path = nx.shortest_path(dual, cara_referencia, i)
#             dual_edges_path = list(zip(dual_path[:-1], dual_path[1:]))
#             edges_path = [dual.edges[(j[0], j[1], 0)]["arista"] 
#                           for j in dual_edges_path]
#             type_edges = [G.edges[j]["tipo"] for j in edges_path]
#             coor_edges = [(G.nodes[n]["coor"], G.nodes[m]["coor"]) 
#                           for n, m in edges_path]
#             reflexiones = [reflexion(n, m) for n, m in coor_edges]
#             comb = list(zip(reflexiones, type_edges))
#             acumulada = ((1, 0, 0, 1), (0, 0))
#             for w in comb:
#                 if w[1] == Kind.CREASE:
#                     acumulada = compon_isometrias(w[0], acumulada,)
#             dicc_iso[i] = acumulada
#     return dicc_iso


def plegar(G, dual, cara_referencia):
    """
    Input: el grafo, el dual, y la cara desde la que empezar.
    Output: un diccionario {cara: isometría}.

    Recorrido en anchura sobre el dual. La cara de referencia recibe la
    identidad. Al cruzar de una cara f resuelta a una vecina g por la arista
    (u,v): si es pliegue, se compone la reflexión de esa recta sobre iso[f];
    si es plana, g hereda la isometría de f.

    Convenio de orden: iso[f] es el camino desde la referencia hasta f, así
    que la reflexión nueva es el último paso y va FUERA:
    compon_isometrias(reflexion, iso[f]).
    """
    iso = {cara_referencia: ((1, 0, 0, 1), (0, 0))}
    # deque es una lista optimizada para sacar por delante
    cola = deque([cara_referencia])

    while cola:
        # popleft() saca el primer elemento, que es lo que hace que el 
        # recorrido sea en anchura. Con una lista normal y pop(0) funcionaría 
        # igual pero más lento.
        f = cola.popleft()
        for g in dual.neighbors(f):
            if g in iso:
                continue
            u, v = dual[f][g][0]["arista"]
            if G[u][v]["tipo"] == Kind.CREASE:
                r = reflexion(G.nodes[u]["coor"], G.nodes[v]["coor"])
                iso[g] = compon_isometrias(iso[f], r)
            else:
                iso[g] = iso[f]
            cola.append(g)

    return iso


def cara_plegada(G, cara, isometria):
    """
    Input: el grafo, una cara y su isometría.
    Output: la lista de vértices de la cara, ya transformados.
    Extrae los vértices del ciclo y les aplica aplica_isometria.
    """
    return [aplica_isometria(isometria, G.nodes[u]["coor"]) for u, _ in cara]




def verifica_consistencia(G, dual, iso):
    """
    Input: el grafo, el dual y el resultado de plegar.
    Output: nada, o error indicando qué arista falla.
    Para cada arista interior (no solo las que usó el recorrido) comprueba que 
    la relación entre las isometrías de sus dos caras es la que debería: si es 
    pliegue, componer la de una con la reflexión debe dar la de la otra; si es 
    plana, deben ser iguales.
    """
    for i in dual.edges:
        f = i[0]
        g = i[1]
        u, v = dual.edges[i]["arista"]
        # las dos caras deben dar la misma imagen sobre la arista común
        for w in (u, v):
            pw = aplica_isometria(iso[f], G.nodes[w]["coor"])
            if pw != aplica_isometria(iso[g], G.nodes[w]["coor"]):
                raise ValueError(
                    f"el papel se rompe en la arista {(u,v)} entre {f} y {g}")
        if G[u][v]["tipo"] == Kind.CREASE:
                    r = reflexion(G.nodes[u]["coor"], G.nodes[v]["coor"])
                    esperada = compon_isometrias(iso[f], r)
        else:
            esperada = iso[f]
        if esperada != iso[g]:
            raise ValueError(
                f"inconsistencia en la arista {(u, v)} entre las caras "
                f"{i[0]} y {i[1]}: se esperaba {esperada}, pero iso[{i[1]}] = {iso[i[1]]}")
        















def main():
    
    from pathlib import Path
    from flatfold.io_fold import read_fold
    # from flatfold.pattern import crea_grafo_dirigido, vecinos_ordenados


    TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

    bomb = TOY / "waterbomb4.fold"

    r = read_fold(bomb)
    G = crea_grafo_dirigido(r)
    
    
    caras = caras_grafo(G)
    interiores, exterior = caras_definitivas(G, caras)
    idx_exterior = caras.index(exterior)
    dual = grafo_dual(G, caras, idx_exterior)
        
    cara_referencia = 0
    iso = plegar(G, dual, cara_referencia)
    # print(dicc)
    
    verifica_consistencia(G, dual, iso)
    

if __name__ == "__main__":
    main()
     






