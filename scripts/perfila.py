# -*- coding: utf-8 -*-
"""
Mide dónde se va el tiempo en una instancia real.
"""
import time
from pathlib import Path
# import argparse

from flatfold.io_fold import read_fold
from flatfold.validate import valida_pendientes, valida_kawasaki
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              verifica_consistencia, cara_plegada,
                              aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)

# RAIZ = Path(__file__).resolve().parent.parent
# # GRIDS = RAIZ / "instances" / "flatfolder" / "grids"
# GRIDS = RAIZ / "instances" / "flatfolder" / "instagram"


RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"



def cronometra(etiqueta, funcion):
    t = time.perf_counter()
    resultado = funcion()
    dt = time.perf_counter() - t
    n = len(resultado) if hasattr(resultado, "__len__") else "-"
    print(f"  {etiqueta:22s} {str(n):>10}  {dt:8.2f} s")
    return resultado


def perfila(ruta):
    print("=" * 55)
    print(ruta.name)

    fi = read_fold(ruta)
    valida_pendientes(fi)
    G = crea_grafo_dirigido(fi)
    valida_kawasaki(G)

    caras = cronometra("caras_grafo", lambda: caras_grafo(G))
    interiores, exterior = caras_definitivas(G, caras)
    idx_ext = caras.index(exterior)
    dual = grafo_dual(G, caras, idx_ext)
    if len(dual) == 0:
        print("  (sin pliegues: nada que perfilar)")
        return
    dic_aristas = aristas_caras(G, caras)

    iso = cronometra("plegar", lambda: plegar(G, dual, min(dual.nodes)))
    cronometra("verifica_consistencia",
               lambda: verifica_consistencia(G, dual, iso) or [])

    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
    print(f"  {'caras interiores':22s} {len(plegadas):>10}")

    pares = cronometra("pares_comparables", lambda: pares_comparables(plegadas))
    trip = cronometra("tripletas", lambda: tripletas(plegadas, pares))
    cr = cronometra("cruces",
                    lambda: cruces(G, caras, plegadas, iso, dic_aristas,
                                   idx_ext))
    col = cronometra("colineales",
                     lambda: colineales(G, caras, plegadas, iso, dic_aristas,
                                        idx_ext))
    tt, tto, oo = col
    print(f"  {'taco-taco':22s} {len(tt):>10}")


# def main():
#     # # de menor a mayor, para poder cortar si se atasca
#     # for nombre in ["003_brown_3x3_A2.fold",
#     #                "182_brown_4x4_C41.fold",
#     #                "226_brown_4x4_C85.fold"]:
#     #     perfila(nombre)
#     for nombre in ["125_brown_seamless4x4.fold",        # la mediana, 81 caras
#                    "357_chan_8x8_Checkerboard.fold"]:   # la grande, 1368
#         perfila(nombre)

# --- configuración: edita estas dos líneas ---
CARPETA = "instagram"
FICHEROS = ["357_chan_8x8_Checkerboard.fold"]
# ---------------------------------------------


def main():
    carpeta = EJEMPLOS / CARPETA
    for nombre in FICHEROS:
        perfila(carpeta / nombre)


if __name__ == "__main__":
    main()


    
    
    
#  Para lanzar: !python scripts/perfila.py instagram 357_chan_8x8_Checkerboard.fold