# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 11:16:30 2026

@author: Trabajador
"""

from flatfold.exact import  sigma
from flatfold.io_fold import Kind
from flatfold.model import literal




def etiquetas_mv(G, fi, iso, dic_aristas, idx_exterior):
    """
    Input:  el grafo, el FoldInput, las isometrías y las adyacencias.
    Output: {numero_de_arista: (par_de_caras, signo)}

    El literal de la etiqueta monte/valle de un pliegue, según la
    Proposición 7.2: para un pliegue entre las caras f y g,

        m_e = X_fg        si sigma(f) = +1
        m_e = 1 - X_fg    si sigma(f) = -1

    donde X_fg vale 1 cuando f está por encima de g. El signo devuelto
    es 0 para el literal directo y 1 para el complementado, igual que
    en literal().

    Solo se incluyen los pliegues interiores: los de frontera no
    doblan y las aristas planas no tienen etiqueta.
    """
    indice = {}
    for i, (u, v) in enumerate(fi.aristas):
        indice[(u, v) if u < v else (v, u)] = i

    etiquetas = {}
    for (u, v), (f, g) in dic_aristas.items():
        if f == idx_exterior or g == idx_exterior:
            continue
        if G[u][v]["tipo"] != Kind.CREASE:
            continue

        # el literal se escribe sobre el par ordenado (f<g)
        if f < g:
            par, s = (f, g), 0 if sigma(iso[f]) == 1 else 1
        else:
            par, s = (g, f), 0 if sigma(iso[g]) == 1 else 1

        etiquetas[indice[(u, v)]] = (par, s)

    return etiquetas

def fija_etiquetas(m, x, etiquetas, prescritas):
    """
    Input:  el modelo, las variables, el resultado de etiquetas_mv y el
            diccionario dicc_montana_valle del FoldInput.
    Output: el número de etiquetas fijadas.

    Para cada pliegue con dirección prescrita en el fichero (M o V),
    fuerza el literal correspondiente a ese valor. Los pliegues sin
    prescribir (U) quedan libres.

    Es el caso de uso "asignación prescrita" del artículo: el modelo
    resulta factible si y solo si existe un orden de capas que produzca
    exactamente esa asignación monte/valle.
    """
    n = 0
    for arista, valor in prescritas.items():
        if arista not in etiquetas:
            continue                      # pliegue de frontera, o arista plana
        (u, v), s = etiquetas[arista]
        lit = literal(x, u, v) if s == 0 else 1 - literal(x, u, v)
        m.addConstr(lit == valor)
        n += 1
    return n


def energia(x, etiquetas, objetivo, pesos=None):
    """
    Input:  las variables, el resultado de etiquetas_mv,
            objetivo — {arista: 0 o 1}, la asignación deseada,
            pesos    — {arista: w}, o None para peso 1 en todas.
    Output: expresión lineal a minimizar.
    
        sum_{e: m^_e = 0} w_e * m_e(x)  +  sum_{e: m^_e = 1} w_e * (1 - m_e(x))
    
    Cada sumando penaliza que la etiqueta de e se desvíe de la deseada.
    Con pesos unitarios, el óptimo es el mínimo número de pliegues que
    hay que invertir respecto de la asignación objetivo.
    """
    expr = 0
    for arista, deseado in objetivo.items():
        if arista not in etiquetas:
            continue
        w = 1 if pesos is None else pesos.get(arista, 1)
        if w == 0:
            continue
        (u, v), s = etiquetas[arista]
        me = literal(x, u, v) if s == 0 else 1 - literal(x, u, v)
        expr += w * me if deseado == 0 else w * (1 - me)
    return expr


