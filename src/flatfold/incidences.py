# -*- coding: utf-8 -*-

from collections import defaultdict
from itertools import combinations
from flatfold.io_fold import Kind
from flatfold.geometry import (se_solapan, atraviesa, caja, interseca_cajas, 
                               tiene_interior, recta_soporte, lado, 
                               solapan_colinealmente)
from flatfold.exact import aplica_isometria



def pares_comparables(caras_plegadas):
    """
    Input:  {cara: polígono plegado}
    Output: lista de pares (f, g) con f < g cuyas imágenes se solapan
            con área positiva.

    Cada par es una variable binaria del modelo: x_fg = 1 significa que la
    cara f está por encima de la g. Los pares que no se solapan no generan
    variable porque nunca compiten por el mismo punto del plano.
    """
    pares = []
    # Todas las posibles combinaciones de 2 caras
    for f, g in combinations(sorted(caras_plegadas), 2):
        if se_solapan(caras_plegadas[f], caras_plegadas[g]):
            pares.append((f, g))
    return pares


# def cruces(G, caras, plegadas, iso, dic_aristas, idx_exterior, pares=None):
#     """
#     Input:  G, caras, plegadas ({cara: polígono}), iso ({cara: isometría}),
#             dic_aristas ({(u,v): (f,g)}), idx_exterior,
#             pares (opcional, solo para comprobar el invariante)
#     Output: lista de (f, g, t) — la cara t atraviesa la arista que separa
#             las caras f y g.

#     Se traduce en la restricción x_ft = x_gt: una cara que cruza la línea de
#     una arista no puede colarse entre los dos trozos de papel unidos por ella.
#     Vale igual para pliegues (el doblez lo impide) y para aristas planas
#     (el papel sigue recto y las dos hojas tendrían que cruzarse).
#     """
#     restricciones = []

#     for (u, v), (f, g) in dic_aristas.items():
#         if f == idx_exterior or g == idx_exterior:
#             continue

#         # imagen plegada de la arista; da igual cuál de las dos caras se use,
#         # porque coinciden sobre ella (lo garantiza verifica_consistencia)
#         segmento = (aplica_isometria(iso[f], G.nodes[u]["coor"]),
#                     aplica_isometria(iso[f], G.nodes[v]["coor"]))

#         for t in plegadas:
#             if t == f or t == g:
#                 continue
#             if atraviesa(segmento, plegadas[t]):
#                 restricciones.append((f, g, t))
#                 if pares is not None:
#                     assert (min(f, t), max(f, t)) in pares, (f, t)
#                     assert (min(g, t), max(g, t)) in pares, (g, t)

#     return restricciones



def cruces(G, caras, plegadas, iso, dic_aristas, idx_exterior, pares=None):
    """
    Output: (taco_tortilla, tortilla_tortilla), listas de (f, g, t).

    Ambas dan la misma restricción X_ft = X_gt; la distinción es de
    clasificación, siguiendo el criterio de Akitaya: si las dos caras de
    la arista se solapan entre sí, la configuración es un taco y el caso
    es taco-tortilla; si no, es tortilla-tortilla.

    Se separan porque las tablas del artículo reportan restricciones por
    familia.
    """
    conjunto = set(pares) if pares is not None else None
    taco_tortilla, tortilla_tortilla = [], []

    for (u, v), (f, g) in dic_aristas.items():
        if f == idx_exterior or g == idx_exterior:
            continue

        segmento = (aplica_isometria(iso[f], G.nodes[u]["coor"]),
                    aplica_isometria(iso[f], G.nodes[v]["coor"]))

        # ¿las dos caras de la arista se solapan? (f1f2 de Akitaya)
        if conjunto is not None:
            es_taco = (min(f, g), max(f, g)) in conjunto
        else:
            es_taco = se_solapan(plegadas[f], plegadas[g])
        destino = taco_tortilla if es_taco else tortilla_tortilla

        for t in plegadas:
            if t == f or t == g:
                continue
            if atraviesa(segmento, plegadas[t]):
                destino.append((f, g, t))
                if conjunto is not None:
                    assert (min(f, t), max(f, t)) in conjunto, (f, t)
                    assert (min(g, t), max(g, t)) in conjunto, (g, t)

    return taco_tortilla, tortilla_tortilla




def tripletas(caras_plegadas, pares):
    """
    Input:  {cara: polígono plegado} y la lista de pares comparables.
    Output: lista de (f, g, h) con f < g < h cuyas tres imágenes tienen
            un punto interior común.

    Se traduce en la restricción de transitividad: si tres caras cubren el
    mismo punto del plano, su orden debe ser total y no cíclico.

    Solo se examinan las tripletas cuyos tres pares se solapan, es decir los
    triángulos del grafo de solapamiento; en las demás falta al menos un
    solapamiento y no puede haber punto común.
    """
    conjunto_pares = set(pares)
    cajas = {f: caja(p) for f, p in caras_plegadas.items()}

    resultado = []
    for f, g, h in combinations(sorted(caras_plegadas), 3):
        if ((f, g) not in conjunto_pares
                or (f, h) not in conjunto_pares
                or (g, h) not in conjunto_pares):
            continue
        if tiene_interior(interseca_cajas([cajas[f], cajas[g], cajas[h]])):
            resultado.append((f, g, h))
    return resultado





def _lado_de_cara(key, poligono):
    """Lado de la recta que ocupa una cara convexa apoyada en ella."""
    for v in poligono:
        s = lado(key, v)
        if s != 0:
            return s
    raise ValueError("cara contenida en la recta")


# def _lado_de_cara(key, poligono):
#     """Lado de la recta que ocupa una cara convexa apoyada en ella."""
#     signos = {lado(key, v) for v in poligono} - {0}
#     if len(signos) > 1:
#         raise ValueError(f"la cara cruza la recta {key}: {poligono}")
#     if not signos:
#         raise ValueError(f"cara contenida en la recta {key}")
#     return signos.pop()


def colineales(G, caras, plegadas, iso, dic_aristas, idx_exterior):
    """
    Output: (taco_taco, taco_tortilla, tortilla_tortilla)

    taco_taco:        listas de cuatro pares de caras; la restricción es que
                      un número PAR de los cuatro literales valga 1.
    taco_tortilla:    pares ((a,t1), (b,t1)); la restricción es la igualdad.
    tortilla_tortilla: pares ((t1,u1), (t2,u2)); igualdad.
    """
    # --- agrupar las aristas interiores por recta soporte de su imagen ---
    grupos = defaultdict(list)
    for (u, v), (f, g) in dic_aristas.items():
        if f == idx_exterior or g == idx_exterior:
            continue
        seg = (aplica_isometria(iso[f], G.nodes[u]["coor"]),
               aplica_isometria(iso[f], G.nodes[v]["coor"]))
        tipo = G[u][v]["tipo"]
        grupos[recta_soporte(seg)].append((seg, tipo, f, g))

    taco_taco, taco_tortilla, tortilla_tortilla = [], [], []

    for key, aristas in grupos.items():
        for i in range(len(aristas)):
            for j in range(i + 1, len(aristas)):
                seg1, tipo1, f1, g1 = aristas[i]
                seg2, tipo2, f2, g2 = aristas[j]
                if solapan_colinealmente(seg1, seg2) is None:
                    continue

                l1f, l1g = _lado_de_cara(key, plegadas[f1]), _lado_de_cara(key, plegadas[g1])
                l2f, l2g = _lado_de_cara(key, plegadas[f2]), _lado_de_cara(key, plegadas[g2])

                # invariante: pliegue -> mismo lado; plana -> lados opuestos
                assert (l1f == l1g) == (tipo1 == Kind.CREASE)
                assert (l2f == l2g) == (tipo2 == Kind.CREASE)

                # --- taco - taco ---
                if tipo1 == Kind.CREASE and tipo2 == Kind.CREASE:
                    if l1f != l2f:
                        continue                       # abren a lados opuestos
                    if len({f1, g1, f2, g2}) < 4:
                        continue                       # comparten cara
                    taco_taco.append([(f1, f2), (g1, f2), (f1, g2), (g1, g2)])

                # --- taco - tortilla ---
                elif tipo1 == Kind.CREASE or tipo2 == Kind.CREASE:
                    if tipo1 == Kind.CREASE:
                        (a, b), lado_taco = (f1, g1), l1f
                        t1 = f2 if l2f == lado_taco else g2
                    else:
                        (a, b), lado_taco = (f2, g2), l2f
                        t1 = f1 if l1f == lado_taco else g1
                    if t1 in (a, b):
                        continue
                    taco_tortilla.append(((a, t1), (b, t1)))

                # --- tortilla - tortilla ---
                else:
                    t1, t2 = (f1, g1) if l1f > 0 else (g1, f1)
                    u1, u2 = (f2, g2) if l2f > 0 else (g2, f2)
                    if t1 == u1 or t2 == u2:
                        continue
                    tortilla_tortilla.append(((t1, u1), (t2, u2)))

    return taco_taco, taco_tortilla, tortilla_tortilla














def main():
    from pathlib import Path
    from flatfold.io_fold import read_fold
    from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                                  caras_definitivas, grafo_dual, plegar,
                                  cara_plegada, aristas_caras)

    TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

    for nombre in ["single_fold", "square", "map_1x3",
                   "waterbomb4", "waterbomb8"]:
        print("=" * 60)
        print(nombre)

        fi = read_fold(TOY / f"{nombre}.fold")
        G = crea_grafo_dirigido(fi)
        caras = caras_grafo(G)
        _, exterior = caras_definitivas(G, caras)
        idx_exterior = caras.index(exterior)
        dual = grafo_dual(G, caras, idx_exterior)
        dic_aristas = aristas_caras(G, caras)

        if len(dual) == 0:
            print("  sin caras interiores conectadas")
            continue

        iso = plegar(G, dual, min(dual.nodes))
        plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
        pares = pares_comparables(plegadas)
        trip = tripletas(plegadas, pares)
        cr = cruces(G, caras, plegadas, iso, dic_aristas, idx_exterior,
                    set(pares))
        tt, tto, oo = colineales(G, caras, plegadas, iso, dic_aristas,
                                 idx_exterior)

        print(f"  caras: {len(plegadas)}   variables: {len(pares)}")
        print(f"  tripletas:          {len(trip)}")
        print(f"  cruces:             {len(cr)}")
        print(f"  taco-taco:          {len(tt)}")
        print(f"  taco-tortilla:      {len(tto)}")
        print(f"  tortilla-tortilla:  {len(oo)}")

        for cuatro in tt:
            print(f"    paridad par sobre {cuatro}")
        for p1, p2 in tto:
            print(f"    igualdad {p1} = {p2}")
        for p1, p2 in oo:
            print(f"    igualdad {p1} = {p2}")


if __name__ == "__main__":
    main()




