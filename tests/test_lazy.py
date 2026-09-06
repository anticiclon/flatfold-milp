# -*- coding: utf-8 -*-
"""
Tests de la transitividad perezosa.

AVISO: el solution pool NO es válido aquí. Con restricciones perezosas
Gurobi puede enumerar soluciones que violan restricciones aún no
añadidas, así que el criterio de "mismo recuento de estados" que se usa
en el resto del proyecto no sirve. Se valida por óptimo y comprobando
que la solución devuelta cumple TODAS las tripletas.
"""
import pytest
from pathlib import Path

import gurobipy as gb
# from flatfold.io_fold import read_fold
# from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
#                               caras_definitivas, grafo_dual, plegar,
#                               cara_plegada, aristas_caras)
# from flatfold.incidences import (pares_comparables, tripletas, cruces,
#                                  colineales)
from flatfold.objective import etiquetas_mv, energia
from flatfold.model import (construye_modelo, añade_restricciones,
                            resuelve_con_lazy, 
                            separa_tripletas, analiza_completo)

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
GRIDS = RAIZ / "instances" / "flatfolder" / "grids"

PEQUEÑA = GRIDS / "004_brown_3x3_A3.fold"
MEDIANA = GRIDS / "182_brown_4x4_C41.fold"


# def monta(ruta):
#     fi = read_fold(ruta)
#     G = crea_grafo_dirigido(fi)
#     caras = caras_grafo(G)
#     _, ext = caras_definitivas(G, caras)
#     idx = caras.index(ext)
#     dual = grafo_dual(G, caras, idx)
#     dic = aristas_caras(G, caras)
#     iso = plegar(G, dual, min(dual.nodes))
#     plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
#     pares = pares_comparables(plegadas)
#     trip = tripletas(plegadas, pares)
#     cr = cruces(G, caras, plegadas, iso, dic, idx)
#     tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)
#     etiq = etiquetas_mv(G, fi, iso, dic, idx)
#     return fi, pares, trip, cr, tt, tto, oo, etiq
def monta(ruta):
    a = analiza_completo(ruta)
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    # return a.fi, a.pares, a.trip, a.cr, a.tt, a.tto, a.oo, etiq
    return a.fi, a.pares, a.trip, a.cr + a.cr_oo, a.tt, a.tto, a.oo, etiq

def modelo(datos, lazy, objetivo=None):
    fi, pares, trip, cr, tt, tto, oo, etiq = datos
    m, x = construye_modelo(pares)
    ind = {p: i for i, p in enumerate(pares)}
    añade_restricciones(m, x, ind, trip, cr, tt, tto, oo,
                        transitividad_lazy=lazy)
    m.Params.OutputFlag = 0
    if objetivo is not None:
        m.setObjective(energia(x, etiq, objetivo), gb.GRB.MINIMIZE)
    if lazy:
        resuelve_con_lazy(m)
    else:
        m.optimize()
    return m, x


# ===== corrección ==========================================================

@pytest.mark.parametrize("ruta", [TOY / "waterbomb4.fold",
                                  TOY / "map_1x5.fold",
                                  PEQUEÑA, MEDIANA])
def test_la_solucion_cumple_todas_las_tripletas(ruta):
    """
    Lo que valida la separación: la solución final no puede violar
    ninguna tripleta, tampoco las que nunca se añadieron.
    """
    datos = monta(ruta)
    m, x = modelo(datos, lazy=True)
    assert m.SolCount > 0
    valores = {par: int(round(v.X)) for par, v in x.items()}
    assert separa_tripletas(valores, datos[2]) == []


@pytest.mark.parametrize("ruta", [TOY / "waterbomb4.fold", PEQUEÑA, MEDIANA])
def test_mismo_optimo_que_la_estatica(ruta):
    """Con objetivo, el valor óptimo debe coincidir."""
    datos = monta(ruta)
    fi = datos[0]
    invertida = {a: 1 - v for a, v in fi.dicc_montana_valle.items()}

    m1, _ = modelo(datos, lazy=False, objetivo=invertida)
    m2, _ = modelo(datos, lazy=True, objetivo=invertida)
    assert m1.SolCount > 0 and m2.SolCount > 0
    assert int(round(m1.ObjVal)) == int(round(m2.ObjVal))


def test_infactible_sigue_siendo_infactible():
    """waterbomb8 con su asignación M/V no es realizable."""
    from flatfold.objective import fija_etiquetas
    datos = monta(TOY / "waterbomb8.fold")
    fi, pares, trip, cr, tt, tto, oo, etiq = datos
    m, x = construye_modelo(pares)
    ind = {p: i for i, p in enumerate(pares)}
    añade_restricciones(m, x, ind, trip, cr, tt, tto, oo,
                        transitividad_lazy=True)
    fija_etiquetas(m, x, etiq, fi.dicc_montana_valle)
    m.Params.OutputFlag = 0
    resuelve_con_lazy(m)
    assert m.SolCount == 0


# ===== el modelo es más pequeño ============================================

def test_menos_restricciones_iniciales():
    datos = monta(MEDIANA)
    trip = datos[2]
    m1, _ = modelo(datos, lazy=False)
    m2, _ = modelo(datos, lazy=True)
    assert m1.NumConstrs - m2.NumConstrs == 2 * len(trip)


def test_solo_se_separa_una_parte():
    """El interés de la variante: la mayoría de tripletas no hace falta."""
    datos = monta(MEDIANA)
    m, _ = modelo(datos, lazy=True)
    separadas = m._n_lazy // 2
    assert 0 < separadas < len(datos[2])


# ===== el contador =========================================================

def test_contador_empieza_en_cero():
    datos = monta(TOY / "waterbomb4.fold")
    fi, pares, trip, cr, tt, tto, oo, etiq = datos
    m, x = construye_modelo(pares)
    ind = {p: i for i, p in enumerate(pares)}
    añade_restricciones(m, x, ind, trip, cr, tt, tto, oo,
                        transitividad_lazy=True)
    assert m._n_lazy == 0
    assert m._tripletas is trip


def test_sin_tripletas_colgadas_falla():
    """resuelve_con_lazy sobre un modelo estático debe protestar."""
    datos = monta(TOY / "waterbomb4.fold")
    fi, pares, trip, cr, tt, tto, oo, etiq = datos
    m, x = construye_modelo(pares)
    ind = {p: i for i, p in enumerate(pares)}
    añade_restricciones(m, x, ind, trip, cr, tt, tto, oo)
    with pytest.raises(AssertionError):
        resuelve_con_lazy(m)
        
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

        