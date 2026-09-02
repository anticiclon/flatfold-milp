# -*- coding: utf-8 -*-
"""
Recorre los ejemplos de Flat-Folder y clasifica hasta dónde llega la cadena:
lectura, validación de hipótesis, topología, plegado y convexidad.

No resuelve ningún modelo. Para eso está resuelve_conjunto.py.
"""
import json
from pathlib import Path
from collections import Counter
import argparse

from flatfold.io_fold import read_fold
from flatfold.validate import valida_pendientes, valida_kawasaki
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              verifica_consistencia, cara_plegada)
from flatfold.geometry import es_convexa






RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"


def carpeta_desde_argumentos(por_defecto="grids"):
    p = argparse.ArgumentParser()
    p.add_argument("carpeta", nargs="?", default=por_defecto,
                   help="subcarpeta de instances/flatfolder (grids, "
                        "instagram, unsatisfiable)")
    args = p.parse_args()
    ruta = EJEMPLOS / args.carpeta
    if not ruta.is_dir():
        p.error(f"no existe {ruta}")
    return ruta











def clasifica(mensaje):
    """Etiqueta un ValueError según su texto."""
    if "retícula" in mensaje:
        return "fuera de retícula"
    if "pendiente" in mensaje:
        return "pendiente no admisible"
    if "Kawasaki" in mensaje or "grado impar" in mensaje:
        return "kawasaki"
    if "área nula" in mensaje or "cara exterior" in mensaje \
            or "áreas no cuadran" in mensaje:
        return "topología"
    if "se rompe" in mensaje or "inconsistencia" in mensaje:
        return "plegado inconsistente"
    return "otro"


def main():
    resultados = Counter()
    discrepancias = []
    sin_convexidad = []
    otros = []
    tamanos = []
    
    carpeta = carpeta_desde_argumentos("unsatisfiable")

    for ruta in sorted(carpeta.glob("*.fold")):
        try:
            fi = read_fold(ruta)
            valida_pendientes(fi)
            G = crea_grafo_dirigido(fi)
            valida_kawasaki(G)

            caras = caras_grafo(G)
            interiores, exterior = caras_definitivas(G, caras)

            n_fichero = len(json.load(open(ruta)).get("faces_vertices", []))
            if n_fichero and len(interiores) != n_fichero:
                resultados["caras discrepan"] += 1
                discrepancias.append((ruta.name, len(interiores), n_fichero))
                continue

            tamanos.append((len(interiores), ruta.name))

            idx_ext = caras.index(exterior)
            dual = grafo_dual(G, caras, idx_ext)
            if len(dual) == 0:
                resultados["sin pliegues"] += 1
                continue

            iso = plegar(G, dual, min(dual.nodes))
            verifica_consistencia(G, dual, iso)

            plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
            no_convexas = [f for f, p in plegadas.items() if not es_convexa(p)]

            if no_convexas:
                resultados["caras no convexas"] += 1
                sin_convexidad.append((ruta.name, len(no_convexas),
                                       len(plegadas)))
            else:
                resultados["plegado ok, todo convexo"] += 1

        except ValueError as e:
            etiqueta = clasifica(str(e))
            resultados[etiqueta] += 1
            if etiqueta in ("otro", "plegado inconsistente", "topología"):
                otros.append((ruta.name, str(e)))
        except Exception as e:
            resultados[f"excepción {type(e).__name__}"] += 1
            otros.append((ruta.name, f"{type(e).__name__}: {e}"))

    print(resultados)

    if discrepancias:
        print(f"\nDiscrepancias con faces_vertices ({len(discrepancias)}):")
        for nombre, mio, suyo in discrepancias[:15]:
            print(f"  {nombre}: yo {mio}, fichero {suyo}")

    if sin_convexidad:
        print(f"\nCon caras no convexas ({len(sin_convexidad)}):")
        for nombre, n, total in sin_convexidad[:15]:
            print(f"  {nombre}: {n} de {total}")

    if otros:
        print(f"\nErrores que merecen mirarse ({len(otros)}):")
        for nombre, msg in otros[:15]:
            print(f"  {nombre}: {msg}")

    if tamanos:
        tamanos.sort()
        print(f"\nCaras interiores: min {tamanos[0][0]} ({tamanos[0][1]}), "
              f"max {tamanos[-1][0]} ({tamanos[-1][1]})")
        medio = tamanos[len(tamanos) // 2]
        print(f"  mediana: {medio[0]} ({medio[1]})")


if __name__ == "__main__":
    main()