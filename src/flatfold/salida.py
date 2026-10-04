# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 10:42:30 2026

@author: Trabajador
"""
from flatfold.exact import aplica_isometria as aplica

def exporta_estado_fold(a, valores, ruta):
    """
    Escribe el estado plegado como fichero FOLD, con el orden de capas en
    el campo faceOrders.

    En FOLD, la terna [f, g, s] se lee contra la NORMAL de g: s = +1 si f
    queda del lado al que apunta esa normal. La normal se invierte cuando
    la cara está volteada, de modo que la traducción desde X_fg, que vale
    1 cuando f está por encima de g en sentido absoluto, es

        s = (2 * X_fg - 1) * sigma(g)

    Las caras se renumeran, porque la exterior no se escribe.
    """
    import json
    from flatfold.exact import sigma

    interiores = sorted(a.plegadas)
    nuevo = {f: i for i, f in enumerate(interiores)}

    # vértices usados por las caras interiores, en coordenadas plegadas
    coords, idx_v = [], {}
    caras_fv = []
    for f in interiores:
        fila = []
        for p in a.plegadas[f]:
            if p not in idx_v:
                idx_v[p] = len(coords)
                coords.append(list(p))
            fila.append(idx_v[p])
        caras_fv.append(fila)

    ordenes = []
    for (f, g), v in valores.items():
        if f not in nuevo or g not in nuevo:
            continue
        ordenes.append([nuevo[f], nuevo[g],
                        (2 * v - 1) * sigma(a.iso[g])])

    doc = {
        "file_spec": 1.1,
        "file_creator": "flatfold-milp",
        "file_classes": ["singleModel"],
        "frame_attributes": ["2D", "convexFaces"],
        "vertices_coords": coords,
        "faces_vertices": caras_fv,
        "faceOrders": ordenes,
    }
    with open(ruta, "w", encoding="utf-8") as h:
        json.dump(doc, h, indent=1)
    return ruta




def verifica_solucion(a, valores, tol_tripletas=None):
    """
    Comprueba que un orden de capas es físicamente válido, recalculando
    las incidencias desde la geometría plegada.

    Input:  a — el Analisis;  valores — {(f,g): 0 o 1} con f < g
    Output: lista de violaciones; vacía si la solución es válida.

    Deliberadamente NO usa a.trip, a.cr, a.tt ni las demás listas: las
    recalcula con los predicados geométricos. Así comprueba un segundo
    camino, y detecta tanto errores en las funciones que generan las
    incidencias como errores de signo en la traducción al modelo, que
    son la clase de fallo que ni los recuentos de estados ni las
    comprobaciones de factibilidad pueden revelar.

    No es independiente de los predicados geométricos de geometry.py,
    que comparte con el generador.
    """
    from itertools import combinations
    from flatfold.geometry import (caja, interseca_cajas, tiene_interior,
                                   se_solapan, atraviesa, recta_soporte,
                                   solapan_colinealmente)
    from flatfold.exact import aplica_isometria
    from collections import defaultdict

    P = a.plegadas
    fallos = []

    def X(u, v):
        if u == v:
            raise ValueError("literal sobre la misma cara")
        return valores[(u, v)] if u < v else 1 - valores[(v, u)]

    def solapan(u, v):
        return se_solapan(P[u], P[v])

    # --- (T) tres caras sobre un punto común: orden total, no cíclico ---
    cajas = {f: caja(p) for f, p in P.items()}
    for f, g, h in combinations(sorted(P), 3):
        if not (solapan(f, g) and solapan(f, h) and solapan(g, h)):
            continue
        if not tiene_interior(interseca_cajas([cajas[f], cajas[g],
                                               cajas[h]])):
            continue
        s = X(f, g) + X(g, h) + X(h, f)
        if not 1 <= s <= 2:
            fallos.append(("ciclo", (f, g, h), s))

    # --- (E) una cara que atraviesa la imagen de una arista ---
    segmentos = {}
    for (u, v), (f, g) in a.dic_aristas.items():
        if f == a.idx_ext or g == a.idx_ext:
            continue
        seg = (aplica_isometria(a.iso[f], a.G.nodes[u]["coor"]),
               aplica_isometria(a.iso[f], a.G.nodes[v]["coor"]))
        segmentos[(u, v)] = (seg, f, g)
        for t in P:
            if t in (f, g):
                continue
            if atraviesa(seg, P[t]) and X(f, t) != X(g, t):
                fallos.append(("cruce", (f, g, t), (X(f, t), X(g, t))))

    # --- (S) aristas colineales cuyas imágenes se solapan ---
    grupos = defaultdict(list)
    for arista, (seg, f, g) in segmentos.items():
        grupos[recta_soporte(seg)].append((seg, f, g))

    for recta, lista in grupos.items():
        for i in range(len(lista)):
            for j in range(i + 1, len(lista)):
                s1, f1, g1 = lista[i]
                s2, f2, g2 = lista[j]
                if solapan_colinealmente(s1, s2) is None:
                    continue
                if len({f1, g1, f2, g2}) < 4:
                    continue
                # criterio de Akitaya: la clasificación la dan los
                # solapamientos entre caras, no el tipo de arista
                if solapan(f1, g1) and solapan(f1, f2) and solapan(f1, g2):
                    # taco-taco: paridad par de los cuatro literales
                    par = (X(f1, f2) + X(g1, f2)
                           + X(f1, g2) + X(g1, g2))
                    if par % 2 != 0:
                        fallos.append(("taco-taco",
                                       (f1, g1, f2, g2), par))

    return fallos



def verifica_punto_a_punto(a, valores, paso=1, max_puntos=2_000_000):
    """
    Comprueba el apilamiento punto a punto, en lugar de por pares de caras.

    Input:  a — el Analisis;  valores — {(f,g): 0 o 1} con f < g
            paso — separación de la retícula de muestreo, en unidades de
                   media coordenada del plegado
    Output: lista de violaciones; vacía si no se encuentra ninguna.

    Esta es la definición original de estado plegado, anterior a la
    reducción a comparaciones entre pares: sobre cada punto del plano, las
    caras que lo cubren deben estar totalmente ordenadas, y a lo largo de
    cada pliegue ninguna cara puede colarse entre las dos que el pliegue
    une.

    Lo que hace útil esta comprobación es que NO usa ninguna función de
    geometry.py: el único predicado es un test de punto en polígono
    convexo por signos del producto cruzado, con aritmética entera. Las
    demás validaciones del proyecto recalculan las incidencias, pero con
    los mismos predicados que el generador, de modo que un error en ellos
    les resultaría invisible.

    Es un muestreo, no una comprobación exhaustiva: cubre las celdas del
    arreglo que contienen algún punto de la retícula. Con paso pequeño y
    patrones pequeños es efectivamente exhaustivo; en patrones grandes
    puede dejar celdas estrechas sin visitar.
    """
    # se trabaja al doble de escala para que los puntos medios de la
    # retícula del plegado caigan también en coordenadas enteras
    P = {f: [(2 * x, 2 * y) for x, y in pol]
         for f, pol in a.plegadas.items()}

    def cubre(p, pol):
        """Cierto si p está en el INTERIOR de pol, convexo. El borde no
        cuenta: dos caras que solo se tocan no compiten por un punto."""
        x, y = p
        n = len(pol)
        s = 0
        for i in range(n):
            (x1, y1), (x2, y2) = pol[i], pol[(i + 1) % n]
            cr = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)
            if cr == 0:
                return False
            s += 1 if cr > 0 else -1
        return abs(s) == n

    def X(u, v):
        return valores[(u, v)] if u < v else 1 - valores[(v, u)]

    def desordenadas(caras):
        """Las claves de un orden total son 0..k-1 distintas; un ciclo
        produce claves repetidas."""
        if len(caras) < 3:
            return False
        clave = {c: sum(1 for o in caras if o != c and X(o, c) == 1)
                 for c in caras}
        return sorted(clave.values()) != list(range(len(caras)))

    fallos = []

    # --- 1) sobre cada punto de la retícula, orden total -------------------
    xs = [x for pol in P.values() for x, _ in pol]
    ys = [y for pol in P.values() for _, y in pol]
    ancho = (max(xs) - min(xs)) // paso + 1
    alto = (max(ys) - min(ys)) // paso + 1
    if ancho * alto > max_puntos:
        raise ValueError(
            f"la retícula tendría {ancho * alto} puntos; sube el paso "
            f"o el límite")

    for x in range(min(xs), max(xs) + 1, paso):
        for y in range(min(ys), max(ys) + 1, paso):
            cubren = [f for f, pol in P.items() if cubre((x, y), pol)]
            if desordenadas(cubren):
                fallos.append(("punto-ciclo", (x, y), tuple(cubren)))

    # --- 2) a lo largo de cada pliegue, nadie entre sus dos caras ----------
    for (u, v), (f, g) in a.dic_aristas.items():
        if f == a.idx_ext or g == a.idx_ext:
            continue
        # imagen del segmento, al doble de escala
        (x1, y1), (x2, y2) = (
            aplica(a.iso[f], a.G.nodes[u]["coor"]),
            aplica(a.iso[f], a.G.nodes[v]["coor"]))
        x1, y1, x2, y2 = 2 * x1, 2 * y1, 2 * x2, 2 * y2
        n = max(abs(x2 - x1), abs(y2 - y1)) // paso
        for i in range(1, n):            # interior del segmento
            p = (x1 + (x2 - x1) * i // n, y1 + (y2 - y1) * i // n)
            for t in P:
                if t in (f, g) or not cubre(p, P[t]):
                    continue
                if X(f, t) != X(g, t):
                    fallos.append(("punto-cruce", p, (f, g, t)))

    return fallos










