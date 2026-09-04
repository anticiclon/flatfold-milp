# -*- coding: utf-8 -*-
from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido)
from flatfold.pattern import vecinos_ordenados
from flatfold.validate import vertices_interiores, minimos_locales, sectores
from pathlib import Path
# from collections import Counter
from compara_presolve import geometria
from flatfold.model import construye_modelo, añade_restricciones



def cuenta_blb(ruta):
    fi = read_fold(ruta)
    G = crea_grafo_dirigido(fi)
    vecinos = vecinos_ordenados(G)
    interiores = vertices_interiores(G)
    n_min = sum(len(minimos_locales(sectores(G, vecinos, v)))
                for v in interiores)
    return len(interiores), n_min


RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"

# --- configuración ---
CARPETA = "grids"
DETALLE = ["003_brown_3x3_A2.fold",
           "182_brown_4x4_C41.fold",
           "226_brown_4x4_C85.fold"]
# ---------------------


def main():
    carpeta = EJEMPLOS / CARPETA
    print(f"Carpeta: {CARPETA}")

    print("\n=== detalle ===")
    for nombre in DETALLE:
        ruta = carpeta / nombre
        if ruta.exists():
            n_int, n_min = cuenta_blb(ruta)
            print(f"  {nombre}: {n_int} vértices interiores, "
                  f"{n_min} mínimos locales")

    print(f"\n=== resumen sobre {CARPETA} ===")
    total_int = total_min = 0
    con_blb = sin_blb = fallos = 0
    for ruta in sorted(carpeta.glob("*.fold")):
        try:
            n_int, n_min = cuenta_blb(ruta)
        except (ValueError, KeyError, AssertionError):
            fallos += 1
            continue
        total_int += n_int
        total_min += n_min
        if n_min:
            con_blb += 1
        else:
            sin_blb += 1

    print(f"  vértices interiores: {total_int}")
    print(f"  mínimos locales:     {total_min}")
    print(f"  instancias con al menos uno: {con_blb}")
    print(f"  instancias sin ninguno:      {sin_blb}")
    print(f"  descartadas: {fallos}")
    
    # print("\n=== candidatas para tests: pocas variables y algún BLB ===")
    # candidatas = []
    # for ruta in sorted(carpeta.glob("*.fold")):
    #     try:
    #         _, n_min = cuenta_blb(ruta)
    #         if n_min == 0:
    #             continue
    #         datos = geometria(ruta)
    #         if datos:
    #             candidatas.append((len(datos[1]), n_min, ruta.name))
    #     except Exception:
    #         pass
    # candidatas.sort()
    # for n_vars, n_min, nombre in candidatas[:10]:
    #     print(f"  {nombre}: {n_vars} variables, {n_min} BLB")
        
    # from compara_presolve import geometria
    
    print('-------------')
    
    datos = geometria(carpeta / "004_brown_3x3_A3.fold")
    fi, pares, trip, cr, tt, tto, oo, etiq = datos
    m, x = construye_modelo(pares)
    indice_pares = {par: i for i, par in enumerate(pares)}
    añade_restricciones(m, x, indice_pares, trip, cr, tt, tto, oo)
    m.Params.OutputFlag = 0
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = 100000
    m.optimize()
    print("estados:", m.SolCount)


if __name__ == "__main__":
    main()

