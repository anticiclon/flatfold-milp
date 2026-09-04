# -*- coding: utf-8 -*-
"""
Compara resolver el modelo completo frente al reducido por el presolve.

Mide, por instancia y por modo:
    tiempos de geometría, presolve y solver, con y sin reducción,
    tamaños antes y después, y número de componentes.

Comprueba además que ambos caminos dan el mismo resultado: si difieren,
el presolve altera el conjunto de soluciones y hay un error.
"""
import time
import csv
from pathlib import Path

import gurobipy as gb
from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.objective import etiquetas_mv, energia
from flatfold.model import (construye_modelo, añade_restricciones,
                            resuelve_por_componentes)
from flatfold.presolve import (literales_a_ecuaciones, reduce_completo,
                               componentes, energia_reducida,
                               reparte_por_componentes, Infactible)

RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"

# --- configuración ---
CARPETA = "instagram"
MODO = "optimizacion"      # factibilidad | optimizacion
MAX_VARS = None          # salta las mayores; None para no saltar
TIEMPO = 2*3600
PRESOLVE_GUROBI = 0     # -1 automático, 0 desactivado, 1 conservador, 2 agresivo
# SOLO = [
#     "321_maeng_Dicranocephalus_wallichii.fold",
#     "343_maeng_Orchid_Mantis.fold",
#     "351_ku_Three_Headed_Dragon.fold",
#     "357_chan_8x8_Checkerboard.fold",
#     "313_maeng_Shrimp.fold",
# ]
SOLO = ["321_maeng_Dicranocephalus_wallichii.fold",
        "343_maeng_Orchid_Mantis.fold",
        "351_ku_Three_Headed_Dragon.fold",
        "357_chan_8x8_Checkerboard.fold"]
SALIDA = RAIZ / "experiments" / f"presolve_{CARPETA}_{MODO}.csv"
# ---------------------


def geometria(ruta):
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
    etiq = etiquetas_mv(G, fi, iso, dic, idx)
    return fi, pares, trip, cr, tt, tto, oo, etiq


def sin_presolve(fi, pares, trip, cr, tt, tto, oo, etiq):
    """Modelo completo. Devuelve (valor, tiempo, estado)."""
    t = time.perf_counter()
    m, x = construye_modelo(pares)
    añade_restricciones(m, x, trip, cr, tt, tto, oo)
    m.Params.OutputFlag = 0
    m.Params.IntegralityFocus = 1
    m.Params.MIPGap = 0
    m.Params.TimeLimit = TIEMPO
    m.Params.Presolve = PRESOLVE_GUROBI
    if MODO == "optimizacion":
        invertida = {a: 1 - v for a, v in fi.dicc_montana_valle.items()}
        m.setObjective(energia(x, etiq, invertida), gb.GRB.MINIMIZE)
    else:
        m.Params.SolutionLimit = 1
    m.optimize()
    dt = time.perf_counter() - t

    if m.SolCount == 0:
        # return None, dt, m.Status
        return None, dt, m.Status
    valor = int(round(m.ObjVal)) if MODO == "optimizacion" else 1
    # return valor, dt, m.Status
    return valor, dt, m.Status, m.NodeCount



def con_presolve(fi, pares, trip, cr, tt, tto, oo, etiq):
    """Modelo reducido. Devuelve (valor, t_presolve, t_solver, n_sup, n_comp)."""
    t = time.perf_counter()
    ig, par, fij = literales_a_ecuaciones(pares, cr, tto, oo, tt)
    mapeo, restantes, vivas, _ = reduce_completo(
        len(pares), ig, par, fij, trip, pares)
    comps = componentes(mapeo, restantes, trip, pares)

    if MODO == "optimizacion":
        invertida = {a: 1 - v for a, v in fi.dicc_montana_valle.items()}
        terminos, const = energia_reducida(mapeo, pares, etiq, invertida)
    else:
        terminos, const = {}, 0

    datos = reparte_por_componentes(comps, vivas, restantes, terminos)
    t_pre = time.perf_counter() - t

    n_sup = sum(len(d["raices"]) for d in datos)

    t = time.perf_counter()
    valor, _, nodes_con = resuelve_por_componentes(datos, const, modo=MODO, tiempo=TIEMPO,
                                        presolve_gurobi=PRESOLVE_GUROBI)
    t_sol = time.perf_counter() - t

    if MODO == "factibilidad":
        valor = 1 if valor else None
    elif valor is not None:
        valor = int(round(valor))
    return valor, t_pre, t_sol, n_sup, len(comps), nodes_con


def main():
    carpeta = EJEMPLOS / CARPETA
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    print(f"Carpeta: {CARPETA}   modo: {MODO}")

    filas = []
    discrepancias = []
    sin_terminar = []

    for ruta in sorted(carpeta.glob("*.fold")):
        if SOLO and ruta.name not in SOLO:
            continue
        try:
            t = time.perf_counter()
            datos = geometria(ruta)
            t_geo = time.perf_counter() - t
        except (ValueError, Infactible):
            continue
        if datos is None:
            continue

        fi, pares, trip, cr, tt, tto, oo, etiq = datos
        if MAX_VARS and len(pares) > MAX_VARS:
            continue

        # v_sin, t_sin, estado, nodes_sin = sin_presolve(*datos)
        try:
            v_sin, t_sin, estado, nodes_sin = sin_presolve(*datos)
        except gb.GurobiError as e:
            print(f"  {ruta.name}: sin presolve falló ({e})")
            continue
        v_con, t_pre, t_con, n_sup, n_comp, nodes_con = con_presolve(*datos)

        # if v_sin != v_con:
        #     discrepancias.append((ruta.name, v_sin, v_con))
        if v_sin is not None and v_con is not None and v_sin != v_con:
            discrepancias.append((ruta.name, v_sin, v_con))
        elif (v_sin is None) != (v_con is None):
            sin_terminar.append((ruta.name, v_sin, v_con))

        filas.append({
            "instancia": ruta.name,
            "vars": len(pares),
            "supervivientes": n_sup,
            "componentes": n_comp,
            "tripletas": len(trip),
            "nodes_sin": int(nodes_sin),
            "nodes_con": int(nodes_con),
            "t_geometria": round(t_geo, 3),
            "t_presolve": round(t_pre, 3),
            "t_solver_sin": round(t_sin, 3),
            "t_solver_con": round(t_con, 3),
            "valor_sin": v_sin,
            "valor_con": v_con,
        })
        # print(f"  {ruta.name:34s} {len(pares):6d} -> {n_sup:5d} vars, "
        #       f"{n_comp:3d} comps | solver {t_sin:7.2f} vs "
        #       f"{t_pre:6.2f}+{t_con:7.2f}")
        print(
        f"  {ruta.name:34s} "
        f"{len(pares):6d} -> {n_sup:5d} vars, "
        f"{n_comp:3d} comps | "
        f"nodes {nodes_sin:.0f} vs {nodes_con:.0f} | "
        f"solver {t_sin:7.2f} vs {t_pre:6.2f}+{t_con:7.2f}"
        )


    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)

    n = len(filas)
    sin = sum(r["t_solver_sin"] for r in filas)
    con = sum(r["t_presolve"] + r["t_solver_con"] for r in filas)
    print(f"\n{n} instancias")
    print(f"  sin presolve: {sin:.1f} s de solver")
    print(f"  con presolve: {con:.1f} s (presolve + solver)")
    print(f"  {'gana' if con < sin else 'pierde'} el presolve "
          f"por un factor {sin / con:.2f}" if con else "")
    if discrepancias:
        print(f"\n*** DISCREPANCIAS ({len(discrepancias)}) ***")
        for nombre, a, b in discrepancias[:10]:
            print(f"  {nombre}: sin={a}, con={b}")
    if sin_terminar:
        print(f"\nSin terminar por tiempo ({len(sin_terminar)}):")
        for nombre, a, b in sin_terminar:
            print(f"  {nombre}: sin={a}, con={b}")
    


if __name__ == "__main__":
    main()