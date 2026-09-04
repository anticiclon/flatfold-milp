# -*- coding: utf-8 -*-
"""
Tests de la formulación por facetas del politopo de paridad.

Las dos formulaciones —variable de acarreo y facetas— describen el mismo
conjunto de puntos enteros. Lo que cambia es la relajación lineal, no las
soluciones. Ese es el criterio que las valida.
"""
import pytest
from itertools import product
from pathlib import Path

import gurobipy as gb
from flatfold.model import (facetas_paridad, construye_modelo,
                            añade_restricciones, añade_taco_taco,
                            resuelve_por_componentes)
from flatfold.presolve import (literales_a_ecuaciones, reduce_completo,
                               componentes, reparte_por_componentes)

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
GRIDS = RAIZ / "instances" / "flatfolder" / "grids"


# ===== facetas_paridad, aislada ===========================================

def test_numero_de_facetas():
    """2^(n-1): una por cada punto de paridad impar."""
    for n in (3, 4, 5):
        assert len(facetas_paridad(n)) == 2 ** (n - 1)


def test_todas_de_cardinal_impar():
    for n in (3, 4, 5):
        assert all(len(s) % 2 == 1 for s in facetas_paridad(n))


def test_facetas_de_cuatro():
    assert set(facetas_paridad(4)) == {
        (0,), (1,), (2,), (3,),
        (0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)}


def test_facetas_de_tres():
    assert set(facetas_paridad(3)) == {(0,), (1,), (2,), (0, 1, 2)}


def test_una_faceta_por_punto_impar():
    """Cada S excluye exactamente un punto: el que vale 1 en S y 0 fuera."""
    n = 4
    impares = {p for p in product((0, 1), repeat=n) if sum(p) % 2 == 1}
    excluidos = {tuple(1 if i in S else 0 for i in range(n))
                 for S in facetas_paridad(n)}
    assert excluidos == impares


# ===== la formulación deja exactamente los puntos pares ===================

def _cuenta_paridad_aislada(n, facetas):
    """Una sola restricción de paridad sobre n literales, nada más."""
    pares = [(0, i + 1) for i in range(n)]
    m, x = construye_modelo(pares)
    añade_taco_taco(m, x, [pares], facetas=facetas)
    m.Params.OutputFlag = 0
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = 100
    m.optimize()
    return m.SolCount


@pytest.mark.parametrize("n", [3, 4])
def test_ambas_formulaciones_dejan_los_mismos_puntos(n):
    """De 2^n asignaciones sobreviven las 2^(n-1) de paridad par."""
    esperado = 2 ** (n - 1)
    assert _cuenta_paridad_aislada(n, facetas=False) == esperado
    assert _cuenta_paridad_aislada(n, facetas=True) == esperado


# ===== sobre instancias ===================================================

def monta(ruta):
    from flatfold.io_fold import read_fold
    from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                                  caras_definitivas, grafo_dual, plegar,
                                  cara_plegada, aristas_caras)
    from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                     colineales)
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


def cuenta_directo(datos, facetas, limite=50000):
    pares, trip, cr, tt, tto, oo = datos
    m, x = construye_modelo(pares)
    indice = {par: i for i, par in enumerate(pares)}
    añade_restricciones(m, x, indice, trip, cr, tt, tto, oo,
                        facetas_paridad_flag=facetas)
    m.Params.OutputFlag = 0
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.optimize()
    assert m.SolCount < limite, "el pool se llenó"
    return m.SolCount


def cuenta_presolve(datos, facetas, limite=50000):
    pares, trip, cr, tt, tto, oo = datos
    ig, par, fij = literales_a_ecuaciones(pares, cr, tto, oo, tt)
    mapeo, restantes, vivas, _ = reduce_completo(
        len(pares), ig, par, fij, trip, pares)
    comps = componentes(mapeo, restantes, trip, pares)
    d = reparte_por_componentes(comps, vivas, restantes, {})
    n, _, _ = resuelve_por_componentes(d, modo="recuento",
                                       limite_pool=limite, facetas=facetas)
    return n


@pytest.mark.parametrize("nombre", ["waterbomb4", "map_1x5"])
def test_recuento_igual_en_el_modelo_completo(nombre):
    datos = monta(TOY / f"{nombre}.fold")
    assert cuenta_directo(datos, facetas=True) == \
           cuenta_directo(datos, facetas=False)


@pytest.mark.parametrize("nombre", ["waterbomb4", "map_1x5"])
def test_recuento_igual_con_presolve(nombre):
    datos = monta(TOY / f"{nombre}.fold")
    assert cuenta_presolve(datos, facetas=True) == \
           cuenta_presolve(datos, facetas=False)


def test_instancia_real():
    """004: 104 variables, 13122 estados."""
    datos = monta(GRIDS / "004_brown_3x3_A3.fold")
    assert cuenta_directo(datos, facetas=False) == 13122
    assert cuenta_directo(datos, facetas=True) == 13122


def test_facetas_no_añaden_variables():
    """La ventaja de la formulación: cero variables enteras nuevas."""
    datos = monta(TOY / "waterbomb4.fold")
    pares, trip, cr, tt, tto, oo = datos

    m1, x1 = construye_modelo(pares)
    ind = {p: i for i, p in enumerate(pares)}
    añade_restricciones(m1, x1, ind, trip, cr, tt, tto, oo,
                        facetas_paridad_flag=False)
    m1.update()

    m2, x2 = construye_modelo(pares)
    añade_restricciones(m2, x2, ind, trip, cr, tt, tto, oo,
                        facetas_paridad_flag=True)
    m2.update()

    assert m2.NumVars == len(pares)          # ninguna auxiliar
    assert m1.NumVars == len(pares) + len(tt)
    assert m2.NumConstrs > m1.NumConstrs     # a cambio, más restricciones
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

    
    
    
    
    