# -*- coding: utf-8 -*-
"""
Resuelve las instancias de una carpeta y comprueba la factibilidad.

En grids/ e instagram/ todas deben salir FACTIBLES: son patrones que
Flat-Folder pliega. En unsatisfiable/ deben salir INFACTIBLES.
"""
import time
from pathlib import Path
from collections import Counter

import gurobipy as gb
from flatfold.model import analiza, construye_modelo, añade_restricciones

RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"

# --- configuración: edita esta línea ---
CARPETA = "unsatisfiable"          # grids | instagram | unsatisfiable
LIMITE_SEGUNDOS = 300
# ---------------------------------------


def factible(datos, tiempo_limite=LIMITE_SEGUNDOS):
    """¿Existe algún orden de capas válido? Para en la primera solución."""
    pares, trip, cr, tt, tto, oo = datos
    m, x = construye_modelo(pares)
    añade_restricciones(m, x, trip, cr, tt, tto, oo)
    m.Params.OutputFlag = 0
    m.Params.SolutionLimit = 1
    m.Params.TimeLimit = tiempo_limite
    m.optimize()
    return m.SolCount > 0, m.Status, m.Runtime


def main():
    carpeta = EJEMPLOS / CARPETA
    if not carpeta.is_dir():
        raise SystemExit(f"no existe {carpeta}")

    ficheros = sorted(carpeta.glob("*.fold"))
    print(f"Carpeta: {CARPETA}  ({len(ficheros)} ficheros)")

    resultados = Counter()
    problemas = []
    lentas = []

    for ruta in ficheros:
        try:
            t0 = time.perf_counter()
            datos = analiza(ruta)
            t_gen = time.perf_counter() - t0

            if not datos[0]:                       # sin variables
                resultados["sin pliegues"] += 1
                continue

            ok, estado, t_solve = factible(datos)
            if ok:
                resultados["FACTIBLE"] += 1
            elif estado == gb.GRB.INFEASIBLE:
                resultados["INFACTIBLE"] += 1
                problemas.append((ruta.name, "infactible", len(datos[0])))
            else:
                resultados[f"sin resolver (estado {estado})"] += 1
                problemas.append((ruta.name, f"estado {estado}",
                                  len(datos[0])))

            if t_gen + t_solve > 5:
                lentas.append((t_gen + t_solve, ruta.name,
                               len(datos[0]), t_gen, t_solve))

        except ValueError:
            resultados["rechazada"] += 1
        except Exception as e:
            resultados[f"excepción {type(e).__name__}"] += 1
            problemas.append((ruta.name, f"{type(e).__name__}: {e}", 0))

    print(resultados)

    if problemas:
        print(f"\nCasos a revisar ({len(problemas)}):")
        for nombre, motivo, nvars in problemas[:20]:
            print(f"  {nombre}: {motivo}  ({nvars} variables)")

    if lentas:
        lentas.sort(reverse=True)
        print("\nMás lentas:")
        for total, nombre, nvars, tg, ts in lentas[:10]:
            print(f"  {nombre}: {total:.1f}s  "
                  f"(geometría {tg:.1f}, solver {ts:.1f}, {nvars} vars)")


if __name__ == "__main__":
    main()