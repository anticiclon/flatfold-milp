# -*- coding: utf-8 -*-
"""
Tests del teorema de Maekawa como restricción de refuerzo.

Es implicada (Corolario 11.4): no cambia el conjunto de soluciones. Ese
es el criterio que valida la formulación. A diferencia de las otras tres
variantes, añade una variable binaria por vértice.
"""
import pytest
from pathlib import Path

import gurobipy as gb
from flatfold.io_fold import read_fold, Kind
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras, vecinos_ordenados)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.validate import vertices_interiores
from flatfold.objective import (etiquetas_mv, restricciones_maekawa,
                                fija_etiquetas)
from flatfold.model import construye_modelo, añade_restricciones

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
GRIDS = RAIZ / "instances" / "flatfolder" / "grids"

PEQUEÑA = GRIDS / "004_brown_3x3_A3.fold"
ESTADOS_PEQUEÑA = 13122


def monta(ruta):
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
    etiq = etiquetas_mv(G, fi, iso, dic, idx)
    mk = restricciones_maekawa(G, fi, iso, dic, idx, etiq)
    return dict(fi=fi, G=G, pares=pares, trip=trip, cr=cr, tt=tt,
                tto=tto, oo=oo, etiq=etiq, mk=mk)


def construye(d, con_maekawa):
    m, x = construye_modelo(d["pares"])
    ind = {p: i for i, p in enumerate(d["pares"])}
    añade_restricciones(m, x, ind, d["trip"], d["cr"], d["tt"],
                        d["tto"], d["oo"], maekawa=d["mk"],
                        con_maekawa=con_maekawa)
    m.Params.OutputFlag = 0
    return m, x


def cuenta(d, con_maekawa, limite=50000):
    m, x = construye(d, con_maekawa)
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.optimize()
    assert m.SolCount < limite, "el pool se llenó"
    return m.SolCount


# ===== restricciones_maekawa ==============================================

def test_un_vertice_en_waterbomb4():
    d = monta(TOY / "waterbomb4.fold")
    assert len(d["mk"]) == 1
    assert len(d["mk"][0]) == 4          # grado 4


def test_un_vertice_en_waterbomb8():
    d = monta(TOY / "waterbomb8.fold")
    assert len(d["mk"]) == 1
    assert len(d["mk"][0]) == 8


def test_sin_vertices_interiores_no_hay_maekawa():
    for nombre in ["map_1x3", "map_1x5", "single_fold"]:
        d = monta(TOY / f"{nombre}.fold")
        assert d["mk"] == []


def test_todos_los_grados_son_pares():
    d = monta(PEQUEÑA)
    for lits in d["mk"]:
        assert len(lits) % 2 == 0


def test_formato():
    d = monta(TOY / "waterbomb4.fold")
    for lits in d["mk"]:
        for (u, v), s in lits:
            assert u < v
            assert s in (0, 1)


def test_coincide_con_los_vertices_aplicables():
    """Un vértice por cada interior cuyas aristas sean todas pliegues."""
    d = monta(PEQUEÑA)
    G = d["G"]
    vec = vecinos_ordenados(G)
    n = sum(1 for v in vertices_interiores(G)
            if all(G[v][w]["tipo"] == Kind.CREASE for w in vec[v]))
    assert len(d["mk"]) == n


# ===== implicada: no cambia el conjunto de soluciones =====================

def test_recuento_base():
    d = monta(PEQUEÑA)
    assert cuenta(d, con_maekawa=False) == ESTADOS_PEQUEÑA


def test_maekawa_no_cambia_el_recuento():
    """
    Es implicada: no puede eliminar estados. Si el recuento baja, la
    formulación de |M - V| = 2 está mal.
    """
    d = monta(PEQUEÑA)
    assert cuenta(d, con_maekawa=True) == ESTADOS_PEQUEÑA


@pytest.mark.parametrize("nombre", ["waterbomb4"])
def test_recuento_igual_en_juguete(nombre):
    d = monta(TOY / f"{nombre}.fold")
    assert cuenta(d, con_maekawa=True) == cuenta(d, con_maekawa=False)


# ===== captura la condición ===============================================

def test_waterbomb8_sigue_infactible_con_maekawa():
    """
    Su asignación M,V,M,V,M,V,M,V da |M - V| = 0, que viola Maekawa.
    Imponerla no puede volverlo factible; si saliera factible, la
    fórmula 4s - 2 estaría mal.
    """
    d = monta(TOY / "waterbomb8.fold")
    for con in (False, True):
        m, x = construye(d, con_maekawa=con)
        fija_etiquetas(m, x, d["etiq"], d["fi"].dicc_montana_valle)
        m.optimize()
        assert m.SolCount == 0


def test_waterbomb4_es_factible_con_su_asignacion():
    """Contraste: aquí sí se cumple Maekawa y debe seguir siendo factible."""
    d = monta(TOY / "waterbomb4.fold")
    m, x = construye(d, con_maekawa=True)
    fija_etiquetas(m, x, d["etiq"], d["fi"].dicc_montana_valle)
    m.optimize()
    assert m.SolCount > 0


# ===== el coste: variables nuevas =========================================

def test_maekawa_anade_una_binaria_por_vertice():
    """La única de las cuatro variantes que añade variables."""
    d = monta(PEQUEÑA)
    m1, _ = construye(d, con_maekawa=False)
    m2, _ = construye(d, con_maekawa=True)
    m1.update()
    m2.update()
    assert m2.NumVars - m1.NumVars == len(d["mk"])
    assert m2.NumConstrs - m1.NumConstrs == len(d["mk"])


# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

