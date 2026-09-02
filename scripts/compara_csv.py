# -*- coding: utf-8 -*-
"""
Compara la generación de restricciones contra grids_data.csv, la salida
de referencia de Flat-Folder.

Columnas contrastadas: faces, variables, taco-taco, taco-tortilla,
tortilla-tortilla, transitivity.

Correspondencia de nombres:
    Flat-Folder            aquí
    faces                  caras interiores
    variables              pares_comparables
    transitivity           tripletas
    taco-tortilla          cruces
    taco-taco              colineales -> taco_taco
    tortilla-tortilla      colineales -> tortilla_tortilla
"""
import csv
from pathlib import Path
from collections import Counter

from flatfold.io_fold import read_fold
from flatfold.validate import valida_pendientes, valida_kawasaki
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)

RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"

# --- configuración ---
CARPETA = "grids"
CSV = "grids_data.csv"
MAX_CARAS = 400          # salta las instancias mayores; None para no saltar
# ---------------------

# CAMPOS = ["faces", "variables", "taco-taco", "taco-tortilla",
#           "tortilla-tortilla", "transitivity"]

CAMPOS = ["vertices", "edges", "faces", "variables", "taco-taco",
          "taco-tortilla", "tortilla-tortilla", "transitivity"]


def mide(ruta):
    """Devuelve los seis recuentos, con los nombres del CSV."""
    fi = read_fold(ruta)
    from collections import Counter
    print("tipos de arista:", Counter(k.name for k in fi.tipos))
    valida_pendientes(fi)
    G = crea_grafo_dirigido(fi)
    valida_kawasaki(G)

    caras = caras_grafo(G)
    interiores, exterior = caras_definitivas(G, caras)
    idx_ext = caras.index(exterior)
    dual = grafo_dual(G, caras, idx_ext)
    dic_aristas = aristas_caras(G, caras)

    if len(dual) == 0:
        return {"faces": len(interiores), 
                "variables": 0, 
                "vertices": 0,
                "edges": 0,
                "taco-taco": 0,
                "taco-tortilla": 0, 
                "tortilla-tortilla": 0, 
                "transitivity": 0}

    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}

    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    cr = cruces(G, caras, plegadas, iso, dic_aristas, idx_ext)
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic_aristas, idx_ext)

    # return {"faces": len(plegadas),
    #         "variables": len(pares),
    #         "taco-taco": len(tt),
    #         "taco-tortilla": len(cr),
    #         "tortilla-tortilla": len(oo),
    #         "transitivity": len(trip)}
    
    return {"vertices": len(fi.vertices),
            "edges": len(fi.aristas),
            "faces": len(plegadas),
            "variables": len(pares),
            "taco-taco": len(tt),
            "taco-tortilla": len(cr) + len(tto),
            "tortilla-tortilla": len(oo),
            "transitivity": len(trip)}


def main():
    carpeta = EJEMPLOS / CARPETA
    referencia = {r["number"]: r
                  for r in csv.DictReader(open(EJEMPLOS / CSV,
                                               encoding="utf-8"))}
    print(f"Carpeta: {CARPETA}   referencia: {CSV} ({len(referencia)} filas)")

    resultados = Counter()
    diferencias = []

    for ruta in sorted(carpeta.glob("*.fold")):
        numero = ruta.name.split("_")[0]
        esperado = referencia.get(numero)
        if esperado is None:
            resultados["sin fila en el CSV"] += 1
            continue
        if MAX_CARAS and int(esperado["faces"]) > MAX_CARAS:
            resultados["saltada por tamaño"] += 1
            continue

        try:
            mio = mide(ruta)
        except ValueError:
            resultados["rechazada"] += 1
            continue
        except Exception as e:
            resultados[f"excepción {type(e).__name__}"] += 1
            diferencias.append((ruta.name, f"{type(e).__name__}: {e}"))
            continue

        fallos = [(c, mio[c], int(esperado[c]))
                  for c in CAMPOS if mio[c] != int(esperado[c])]
        if fallos:
            resultados["DIFIERE"] += 1
            diferencias.append((ruta.name, fallos))
        else:
            resultados["coincide"] += 1

    print(resultados)

    if diferencias:
        print(f"\nDiferencias ({len(diferencias)}):")
        for nombre, detalle in diferencias[:20]:
            if isinstance(detalle, str):
                print(f"  {nombre}: {detalle}")
            else:
                print(f"  {nombre}:")
                for campo, mio, suyo in detalle:
                    print(f"      {campo:20s} yo {mio:>8}  CSV {suyo:>8}")


if __name__ == "__main__":
    main()