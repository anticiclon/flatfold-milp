# -*- coding: utf-8 -*-
"""
Comprobación manual de la transitividad perezosa.

Compara, sobre unas cuantas instancias, la formulación estática con la
separación por callback: mismo resultado, y cuántas tripletas hicieron
falta de verdad.
"""
from pathlib import Path

from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.model import (construye_modelo, añade_restricciones,
                            resuelve_con_lazy, separa_tripletas)

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
GRIDS = RAIZ / "instances" / "flatfolder" / "grids"

INSTANCIAS = [TOY / "waterbomb4.fold",
              TOY / "map_1x5.fold",
              GRIDS / "004_brown_3x3_A3.fold",
              GRIDS / "182_brown_4x4_C41.fold"]


def monta(ruta):
    fi = read_fold(ruta)
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, ext = caras_definitivas(G, caras)
    idx = caras.index(ext)
    dual = grafo_dual(G, caras, idx)
    dic = aristas_caras(G, caras)
    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    cr = cruces(G, caras, plegadas, iso, dic, idx)
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)
    return pares, trip, cr, tt, tto, oo


def prueba(ruta):
    pares, trip, cr, tt, tto, oo = monta(ruta)
    ind = {p: i for i, p in enumerate(pares)}
    print(f"\n{ruta.name}: {len(pares)} variables, {len(trip)} tripletas")

    # --- estática ---
    m1, x1 = construye_modelo(pares)
    añade_restricciones(m1, x1, ind, trip, cr, tt, tto, oo)
    m1.Params.OutputFlag = 0
    m1.optimize()
    print(f"  estática: estado {m1.Status}, "
          f"{m1.NumConstrs} restricciones, {m1.Runtime:.2f} s")

    # --- perezosa ---
    m2, x2 = construye_modelo(pares)
    añade_restricciones(m2, x2, ind, trip, cr, tt, tto, oo,
                        transitividad_lazy=True)
    m2.Params.OutputFlag = 0
    resuelve_con_lazy(m2)
    print(f"  perezosa: estado {m2.Status}, "
          f"{m2.NumConstrs} restricciones iniciales, {m2.Runtime:.2f} s")
    print(f"            separadas {m2._n_lazy // 2} de {len(trip)} tripletas")

    # --- ¿la solución cumple TODAS las tripletas? ---
    if m2.SolCount:
        valores = {par: int(round(v.X)) for par, v in x2.items()}
        violadas = separa_tripletas(valores, trip)
        print(f"            violadas en la solución: {len(violadas)}"
              f"{'  *** MAL ***' if violadas else '  (correcto)'}")
    else:
        print("            sin solución")


def main():
    for ruta in INSTANCIAS:
        if ruta.exists():
            prueba(ruta)
        else:
            print(f"{ruta.name}: no está")


if __name__ == "__main__":
    main()