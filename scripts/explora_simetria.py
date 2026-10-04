# -*- coding: utf-8 -*-
"""
Determina empíricamente cómo se transporta el orden de capas bajo una
simetría del patrón.

Una isometría del papel que aplica el patrón en sí mismo induce una
permutación de las caras. La pregunta es si, para un estado plegado
válido, el orden transportado coincide con el original,

    X_{sigma(f) sigma(g)} = X_{fg},

o con su complemento. La conjetura es que el signo lo decide si la
isometría preserva o invierte la orientación, pero conviene comprobarlo
antes de escribir restricciones que lo supongan.

El procedimiento enumera todos los estados de un patrón simétrico y, para
cada estado y cada simetría, mira cuál de las dos relaciones se cumple.
"""
# from itertools import combinations
from pathlib import Path

# import gurobipy as gb
from flatfold.exact import aplica_isometria
from flatfold.model import (analiza_completo, construye_modelo,
                            añade_restricciones)

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"


def simetrias_del_cuadrado(L):
    """Las ocho isometrías del cuadrado [0,L]^2, con su nombre y el signo
    de su determinante."""
    return [
        ("identidad",      ((1, 0, 0, 1), (0, 0)),   +1),
        ("giro 90",        ((0, -1, 1, 0), (L, 0)),  +1),
        ("giro 180",       ((-1, 0, 0, -1), (L, L)), +1),
        ("giro 270",       ((0, 1, -1, 0), (0, L)),  +1),
        ("reflexion x",    ((-1, 0, 0, 1), (L, 0)),  -1),
        ("reflexion y",    ((1, 0, 0, -1), (0, L)),  -1),
        ("reflexion diag", ((0, 1, 1, 0), (0, 0)),   -1),
        ("reflexion anti", ((0, -1, -1, 0), (L, L)), -1),
    ]


def preserva_el_patron(a, iso):
    """Cierto si la isometría aplica cada arista del patrón en otra del
    mismo tipo."""
    coor = {v: a.G.nodes[v]["coor"] for v in a.G.nodes}
    puntos = {p: v for v, p in coor.items()}
    for u, v in a.G.edges:
        pu = aplica_isometria(iso, coor[u])
        pv = aplica_isometria(iso, coor[v])
        if pu not in puntos or pv not in puntos:
            return False
        if a.G[puntos[pu]][puntos[pv]]["tipo"] != a.G[u][v]["tipo"]:
            return False
    return True


def permutacion_de_caras(a, iso):
    """{cara: cara imagen}, identificando cada cara por el conjunto de sus
    vértices en el patrón sin plegar. Devuelve None si alguna imagen no
    corresponde a ninguna cara."""
    def vertices(f):
        return frozenset(a.G.nodes[h[0]]["coor"] for h in a.caras[f])

    indice = {vertices(f): f for f in a.plegadas}
    perm = {}
    for f in a.plegadas:
        img = frozenset(aplica_isometria(iso, p) for p in vertices(f))
        if img not in indice:
            return None
        perm[f] = indice[img]
    return perm


def estados(a, limite=1000):
    """Todos los órdenes de capas válidos, como lista de {(f,g): 0 o 1}."""
    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr + a.cr_oo, a.tt,
                        a.tto, a.oo)
    m.Params.OutputFlag = 0
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.optimize()
    if m.SolCount >= limite:
        raise RuntimeError("el pool se llenó")

    salida = []
    for k in range(m.SolCount):
        m.Params.SolutionNumber = k
        salida.append({par: int(round(v.Xn)) for par, v in x.items()})
    return salida


def relacion(vals, perm, pares):
    """Compara el orden con el transportado por la permutación.
    Devuelve 'igual', 'complementario' o 'ninguna'."""
    def X(d, u, v):
        return d[(u, v)] if u < v else 1 - d[(v, u)]

    iguales = complementarios = 0
    total = 0
    for f, g in pares:
        pf, pg = perm[f], perm[g]
        if pf == pg:
            continue
        par = (min(pf, pg), max(pf, pg))
        if par not in vals:
            return "ninguna"
        total += 1
        if X(vals, pf, pg) == X(vals, f, g):
            iguales += 1
        else:
            complementarios += 1
    if iguales == total:
        return "igual"
    if complementarios == total:
        return "complementario"
    return "ninguna"


def main(nombre="waterbomb4"):
    a = analiza_completo(TOY / f"{nombre}.fold")
    L = max(max(p) for p in a.fi.vertices)
    print(f"{nombre}: {len(a.plegadas)} caras, {len(a.pares)} pares, "
          f"papel [0,{L}]^2")

    sims = [(n, iso, det) for n, iso, det in simetrias_del_cuadrado(L)
            if preserva_el_patron(a, iso)]
    print(f"simetrías del patrón: {len(sims)} de 8\n")

    todos = estados(a)
    print(f"estados: {len(todos)}\n")

    for nom, iso, det in sims:
        perm = permutacion_de_caras(a, iso)
        if perm is None:
            print(f"{nom:18s} det {det:+d}: no permuta las caras")
            continue
        cuenta = {"igual": 0, "complementario": 0, "ninguna": 0}
        for vals in todos:
            cuenta[relacion(vals, perm, a.pares)] += 1
        print(f"{nom:18s} det {det:+d}: {cuenta}")


if __name__ == "__main__":
    main()