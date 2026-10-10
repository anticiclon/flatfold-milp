# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 11:16:30 2026

@author: Trabajador
"""
import gurobipy as gb

from flatfold.exact import sigma
from flatfold.io_fold import Kind
from flatfold.model import literal, construye_modelo, añade_restricciones
from flatfold.pattern import vecinos_ordenados
from flatfold.validate import (vertices_interiores, sectores,
                                   minimos_locales)



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

def restricciones_blb(G, fi, iso, dic_aristas, idx_exterior, etiquetas):
    """
    Input:  el grafo, el FoldInput, las isometrías, las adyacencias,
            el índice de la cara exterior y las etiquetas de etiquetas_mv.
    Output: lista de ((par1, s1), (par2, s2)) — literales que deben ser
            COMPLEMENTARIOS, es decir m_e1 != m_e2.

    Condición big-little-big: si un sector es mínimo local estricto, los
    dos pliegues que lo delimitan tienen etiquetas opuestas.

    Es una restricción implicada (Corolario 11.4): no cambia el conjunto
    de soluciones. Su interés es doble: puede podar el árbol de búsqueda,
    y al ser afín sobre F2 puede alimentar el presolve y reducir más
    variables.
    """

    vecinos = vecinos_ordenados(G)

    # de par de vértices a número de arista del fichero
    num_arista = {}
    for i, (u, v) in enumerate(fi.aristas):
        num_arista[(u, v) if u < v else (v, u)] = i

    salida = []
    for v in vertices_interiores(G):
        secs = sectores(G, vecinos, v)
        for i in minimos_locales(secs):
            _, w1, w2 = secs[i]
            e1 = num_arista[(v, w1) if v < w1 else (w1, v)]
            e2 = num_arista[(v, w2) if v < w2 else (w2, v)]
            if e1 not in etiquetas or e2 not in etiquetas:
                continue                  # pliegue sin etiqueta: saltar
            salida.append((etiquetas[e1], etiquetas[e2]))
    return salida


def restricciones_maekawa(G, fi, iso, dic_aristas, idx_exterior, etiquetas):
    """
    Input:  el grafo, el FoldInput, las isometrías, las adyacencias,
            el índice de la cara exterior y las etiquetas de etiquetas_mv.
    Output: lista de listas de literales, una por vértice aplicable:
            [[((u,v), s), ...]]

    Teorema de Maekawa: en un vértice interior de un patrón plegable
    plano, la diferencia entre el número de montes y el de valles es
    exactamente 2.

    Solo se emiten los vértices cuyas aristas interiores sean TODAS
    pliegues: una arista plana no es ni monte ni valle, y la condición
    no aplica en esa forma. (Ahora apenas hay aristas planas, pero con
    la subdivisión convexa aparecerán.)

    Es una restricción implicada (Corolario 11.4): no cambia el conjunto
    de soluciones. Su interés es la poda del árbol de búsqueda.
    """
    from flatfold.pattern import vecinos_ordenados
    from flatfold.validate import vertices_interiores

    vecinos = vecinos_ordenados(G)

    num_arista = {}
    for i, (u, v) in enumerate(fi.aristas):
        num_arista[(u, v) if u < v else (v, u)] = i

    salida = []
    for v in vertices_interiores(G):
        ws = vecinos[v]
        if any(G[v][w]["tipo"] != Kind.CREASE for w in ws):
            continue                      # hay aristas planas: no aplica

        lits = []
        for w in ws:
            e = num_arista[(v, w) if v < w else (w, v)]
            if e not in etiquetas:
                lits = None
                break
            lits.append(etiquetas[e])
        if lits is None:
            continue

        assert len(lits) % 2 == 0, (
            f"vértice {v}: grado impar ({len(lits)} pliegues), "
            "viola Maekawa por conteo")
        salida.append(lits)

    return salida

def etiquetas_por_pliegue_original(etiquetas, representantes):
    """
    Input:  etiquetas — {arista_actual: (par, signo)}, de etiquetas_mv
                        sobre el patrón (subdividido o no)
            representantes — {arista_del_fichero: arista_actual}
    Output: {arista_del_fichero: (par, signo)}

    La energía se define sobre los pliegues del fichero, donde viven la
    asignación objetivo y los pesos. Tras subdividir, un pliegue puede
    estar partido en varios trozos; se toma la etiqueta de su
    representante.

    Es correcto porque la familia de cruces obliga a que todos los
    trozos de un mismo pliegue reciban la misma etiqueta (véase
    comprueba_trozos_coinciden, que lo verifica sobre la solución).

    Sin subdividir, representantes es la identidad y esto no cambia nada.
    """
    salida = {}
    for original, actual in representantes.items():
        if actual in etiquetas:
            salida[original] = etiquetas[actual]
    return salida

def comprueba_trozos_coinciden(valores, etiquetas, origen):
    """
    Sobre una solución, todos los trozos de un mismo pliegue deben tener
    la misma etiqueta. Lo garantiza la familia de cruces; comprobarlo es
    un test independiente de que esa familia se generó bien.
    """
    from collections import defaultdict
    por_origen = defaultdict(set)
    for arista, ((u, v), s) in etiquetas.items():
        o = origen.get(arista)
        if o is None:
            continue
        x = valores[(u, v)]
        por_origen[o].add(x ^ s)

    malos = {o: vals for o, vals in por_origen.items() if len(vals) > 1}
    if malos:
        raise ValueError(f"trozos con etiquetas distintas: {malos}")
        


def cuenta_montes(x, etiquetas):
    """
    Input:  x — las variables de capa
            etiquetas — {arista: (par, signo)}, de etiquetas_mv o de
                        etiquetas_por_pliegue_original
    Output: expresión lineal: cuántos pliegues resultan monte.

    Cada m_e es un literal sobre una variable de capa, de modo que su
    suma es lineal en x. Con ella se prescribe el número de montes,
    sum = M, o se acota, y el número de valles es el complemento:
    len(etiquetas) menos esta suma.

    El recuento abarca los pliegues presentes en el diccionario recibido,
    de modo que la elección de cuáles entran la hace quien lo construye.
    """
    return gb.quicksum(
        literal(x, u, v) if s == 0 else 1 - literal(x, u, v)
        for (u, v), s in etiquetas.values())

def distancia_hamming(x, etiquetas, m0):
    """
    Número de pliegues en que la asignación inducida difiere de m0.

    Coincide con la energía de pesos unitarios respecto de m0, y se
    implementa llamándola: la aritmética de signos que traduce el orden
    de capas a etiquetas es delicada y conviene que exista en un solo
    sitio.
    """
    return energia(x, etiquetas, m0)




def enumera_asignaciones(a, limite=1000, tiempo=None):
    """
    Input:  a — el Analisis;  limite — tope de asignaciones a devolver.
    Output: (asignaciones, agotado) donde asignaciones es una lista de
            {arista_del_fichero: 0 o 1} y agotado dice si se encontraron
            todas o se alcanzó el límite.

    Enumera las asignaciones monte/valle DISTINTAS que el patrón admite,
    no los órdenes de capas. Son cosas diferentes: un mismo etiquetado
    puede realizarse por varios apilamientos, de modo que el número de
    asignaciones es en general menor que el de estados plegados, y es el
    que responde a la pregunta de cuántas formas hay de plegar el patrón.

    El procedimiento resuelve repetidamente, añadiendo tras cada solución
    la condición de que la siguiente diste al menos uno de ella. La
    distancia es la de Hamming sobre los pliegues del fichero, lineal en
    las variables de capa, de modo que cada corte es una sola fila.
    """

    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr + a.cr_oo, a.tt,
                        a.tto, a.oo)
    m.Params.OutputFlag = 0
    m.Params.SolutionLimit = 1
    if tiempo is not None:
        m.Params.TimeLimit = tiempo

    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)

    salida = []
    while len(salida) < limite:
        m.optimize()
        if m.SolCount == 0:
            return salida, True              # no quedan más

        asignacion = {}
        for arista, ((u, v), s) in orig.items():
            valor = int(round(x[u, v].X))
            asignacion[arista] = valor ^ s
        salida.append(asignacion)

        # la siguiente debe diferir en al menos un pliegue
        m.addConstr(distancia_hamming(x, orig, asignacion) >= 1)

    return salida, False                     # se alcanzó el límite
        
