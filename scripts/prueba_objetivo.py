# -*- coding: utf-8 -*-
"""
Valida fija_etiquetas y energia sobre el mapa 1x3.
"""
from pathlib import Path
import gurobipy as gb

from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.objective import etiquetas_mv, fija_etiquetas, energia
from flatfold.model import construye_modelo, añade_restricciones

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"


def prepara(nombre):
    fi = read_fold(TOY / f"{nombre}.fold")
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, exterior = caras_definitivas(G, caras)
    idx_ext = caras.index(exterior)
    dual = grafo_dual(G, caras, idx_ext)
    dic = aristas_caras(G, caras)
    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}

    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    cr = cruces(G, caras, plegadas, iso, dic, idx_ext)
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx_ext)

    m, x = construye_modelo(pares)
    añade_restricciones(m, x, trip, cr, tt, tto, oo)
    m.Params.OutputFlag = 0

    etiq = etiquetas_mv(G, fi, iso, dic, idx_ext)
    
    print(f"  {nombre}: {len(pares)} vars, {len(trip)} tripletas, "
      f"{len(cr)} cruces, {len(tt)} taco-taco")
    
    return fi, m, x, etiq


def prueba_prescrita(nombre):
    fi, m, x, etiq = prepara(nombre)
    n = fija_etiquetas(m, x, etiq, fi.dicc_montana_valle)
    m.optimize()
    if m.SolCount == 0:
        print(f"  {nombre:16s} INFACTIBLE ({n} etiquetas fijadas)")
        return
    print(f"  {nombre:16s} FACTIBLE ({n} de "
          f"{len(fi.dicc_montana_valle)} fijadas)")
    print(f"     variables: {{k: int(round(v.X)) for k, v in x.items()}}")
    for arista, ((u, v), s) in sorted(etiq.items()):
        xv = int(round(x[u, v].X))
        me = xv if s == 0 else 1 - xv
        print(f"     arista {arista}: caras {(u,v)} sigma-signo {s}  "
              f"x={xv}  m_e={me}  pedido={fi.dicc_montana_valle.get(arista)}")


def prueba_reparacion(nombre):
    fi, m, x, etiq = prepara(nombre)
    m.setObjective(energia(x, etiq, fi.dicc_montana_valle), gb.GRB.MINIMIZE)
    m.optimize()
    if m.SolCount == 0:
        print(f"  {nombre:16s} reparación: INFACTIBLE")
    else:
        print(f"  {nombre:16s} reparación: {int(round(m.ObjVal))} "
              f"pliegues a invertir")




def main():
    print("etiquetas prescritas (debe ser factible la MV, infactible la MM):")
    prueba_prescrita("map_1x3")
    prueba_prescrita("map_1x3_MM")

    print("\nreparación (mínimo de inversiones):")
    prueba_reparacion("map_1x3")
    prueba_reparacion("map_1x3_MM")


if __name__ == "__main__":
    main()