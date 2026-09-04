# -*- coding: utf-8 -*-
"""
Mide cómo degeneran las tripletas al sustituir por raíces del presolve.

La clave del conteo es (número de literales constantes, raíces repetidas).
El caso sano es (0, 0): tres literales sobre tres variables supervivientes
distintas. Cualquier otro exige tratamiento especial, y este script dice
cuáles aparecen de verdad.
"""
from pathlib import Path

from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.presolve import (literales_a_ecuaciones, Infactible, 
                               reduce_completo)

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
EJEMPLOS = RAIZ / "instances" / "flatfolder"

# --- configuración ---
CARPETA = "grids"        # grids | instagram | unsatisfiable
DETALLE = ["003_brown_3x3_A2.fold",
            "182_brown_4x4_C41.fold",
            "226_brown_4x4_C85.fold"]
# DETALLE = ["125_brown_seamless4x4.fold", "357_chan_8x8_Checkerboard.fold"]
# ---------------------


def monta(ruta):
    fi = read_fold(ruta)
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, ext = caras_definitivas(G, caras)
    idx = caras.index(ext)
    dual = grafo_dual(G, caras, idx)
    if len(dual) == 0:
        return None
    dic = aristas_caras(G, caras)
    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    cr = cruces(G, caras, plegadas, iso, dic, idx)
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)
    return pares, trip, cr, tt, tto, oo


def analiza(ruta):
    datos = monta(ruta)
    if datos is None:
        return None
    pares, trip, cr, tt, tto, oo = datos
    ig, par, fij = literales_a_ecuaciones(pares, cr, tto, oo, tt)
    # mapeo, restantes, _ = reduce(len(pares), ig, par, fij)
    mapeo, restantes, vivas, n_iter = reduce_completo(
    len(pares), ig, par, fij, trip, pares)

    supervivientes = {e[1] for e in mapeo.values() if e[0] == "var"}
    determinadas = sum(1 for e in mapeo.values() if e[0] == "const")

    return {
        "vars": len(pares),
        "supervivientes": len(supervivientes),
        "determinadas": determinadas,
        "paridades_restantes": len(restantes),
        "tripletas": len(trip),
        "iteraciones": n_iter,
        "tripletas_vivas": len(vivas),
    }


def informe(nombre, r):
    if r is None:
        print(f"{nombre}: sin pliegues")
        return
    reduccion = (100 * (1 - r["supervivientes"] / r["vars"])
                 if r["vars"] else 0)
    print(f"{nombre}")
    print(f"  variables {r['vars']} -> supervivientes {r['supervivientes']} "
          f"({reduccion:.0f}% menos), determinadas {r['determinadas']}")
    print(f"  paridades restantes: {r['paridades_restantes']}")
    print(f"  tripletas {r['tripletas']} -> vivas {r['tripletas_vivas']} "
          f"({100 * (1 - r['tripletas_vivas'] / r['tripletas']):.0f}% menos)"
          if r['tripletas'] else "  sin tripletas")
    print(f"  iteraciones del punto fijo: {r['iteraciones']}")
    

def main():
    carpeta = EJEMPLOS / CARPETA
    if not carpeta.is_dir():
        raise SystemExit(f"no existe {carpeta}")

    print("=== instancias de juguete ===")
    for nombre in ["map_1x3", "map_1x5", "waterbomb4", "waterbomb8"]:
        ruta = TOY / f"{nombre}.fold"
        if ruta.exists():
            informe(nombre, analiza(ruta))

    print(f"\n=== detalle de {CARPETA} ===")
    for nombre in DETALLE:
        ruta = carpeta / nombre
        if ruta.exists():
            informe(nombre, analiza(ruta))
        else:
            print(f"{nombre}: no está en {CARPETA}")

    print(f"\n=== resumen sobre {CARPETA} ===")
    fallos = 0
    n_vars = n_sup = n_trip = n_vivas = 0
    peor = (0, None)
    for ruta in sorted(carpeta.glob("*.fold")):
        try:
            r = analiza(ruta)
        except (ValueError, Infactible):
            fallos += 1
            continue
        if not r:
            continue
        n_vars += r["vars"]
        n_sup += r["supervivientes"]
        n_trip += r["tripletas"]
        n_vivas += r["tripletas_vivas"]
        if r["vars"] > peor[0]:
            peor = (r["vars"], ruta.name)

    if n_vars:
        print(f"  variables {n_vars} -> supervivientes {n_sup} "
              f"({100 * (1 - n_sup / n_vars):.0f}% menos)")
    if n_trip:
        print(f"  tripletas {n_trip} -> vivas {n_vivas} "
              f"({100 * (1 - n_vivas / n_trip):.0f}% menos)")
    print(f"  instancias descartadas: {fallos}")
    print(f"  mayor: {peor[1]} con {peor[0]} variables")
    
    
if __name__ == "__main__":
    main()