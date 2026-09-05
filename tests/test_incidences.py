# -*- coding: utf-8 -*-
import pytest
from pathlib import Path

# from flatfold.io_fold import read_fold
# from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              # caras_definitivas, grafo_dual, plegar,
                              # cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.model import analiza_completo

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"

RUTAS = {
    "single": TOY / "single_fold.fold",
    "square": TOY / "square.fold",
    "map13":  TOY / "map_1x3.fold",
    "bomb":   TOY / "waterbomb4.fold",
    "bomb8":  TOY / "waterbomb8.fold",
}



# def monta(nombre):
#     """Devuelve todo lo que consumen las funciones de incidencias."""
#     fi = read_fold(RUTAS[nombre])
#     G = crea_grafo_dirigido(fi)
#     caras = caras_grafo(G)
#     _, exterior = caras_definitivas(G, caras)
#     idx_ext = caras.index(exterior)
#     dual = grafo_dual(G, caras, idx_ext)
#     dic = aristas_caras(G, caras)
#     if len(dual) == 0:
#         return G, caras, {}, {}, dic, idx_ext
#     iso = plegar(G, dual, min(dual.nodes))
#     plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
#     return G, caras, plegadas, iso, dic, idx_ext
def monta(nombre):
    a = analiza_completo(RUTAS[nombre])
    return a.G, a.caras, a.plegadas, a.iso, a.dic_aristas, a.idx_ext

# --- pares comparables -----------------------------------------------------

@pytest.mark.parametrize("nombre,n", [
    ("single", 1), ("square", 0), ("map13", 3),
    ("bomb", 6), ("bomb8", 28),
])
def test_numero_de_pares(nombre, n):
    _, _, plegadas, *_ = monta(nombre)
    assert len(pares_comparables(plegadas)) == n


def test_pares_ordenados_y_sin_repetir():
    _, _, plegadas, *_ = monta("bomb")
    pares = pares_comparables(plegadas)
    assert all(f < g for f, g in pares)
    assert len(set(pares)) == len(pares)


def test_todas_las_caras_de_un_par_existen():
    _, _, plegadas, *_ = monta("bomb8")
    for f, g in pares_comparables(plegadas):
        assert f in plegadas and g in plegadas


# --- tripletas -------------------------------------------------------------

@pytest.mark.parametrize("nombre,n", [
    ("single", 0), ("map13", 1), ("bomb", 4), ("bomb8", 56),
])
def test_numero_de_tripletas(nombre, n):
    _, _, plegadas, *_ = monta(nombre)
    pares = pares_comparables(plegadas)
    assert len(tripletas(plegadas, pares)) == n


def test_tripletas_ordenadas():
    _, _, plegadas, *_ = monta("bomb")
    pares = pares_comparables(plegadas)
    for f, g, h in tripletas(plegadas, pares):
        assert f < g < h


def test_los_tres_pares_de_una_tripleta_existen():
    """Si tres caras cubren un punto, las tres se solapan dos a dos."""
    _, _, plegadas, *_ = monta("bomb8")
    pares = set(pares_comparables(plegadas))
    for f, g, h in tripletas(plegadas, pares):
        assert (f, g) in pares and (f, h) in pares and (g, h) in pares


# --- cruces ----------------------------------------------------------------

@pytest.mark.parametrize("nombre,n", [
    ("single", 0), ("map13", 0), ("bomb", 0), ("bomb8", 0),
])
def test_numero_de_cruces(nombre, n):
    """
    En las instancias de juguete ninguna cara atraviesa la imagen de una
    arista: todas caen sobre los bordes del montón.
    """
    G, caras, plegadas, iso, dic, idx = monta(nombre)
    assert len(cruces(G, caras, plegadas, iso, dic, idx)) == n


def test_cruces_implican_solapamiento():
    """Si t atraviesa la arista de f y g, t se solapa con las dos."""
    G, caras, plegadas, iso, dic, idx = monta("bomb8")
    pares = set(pares_comparables(plegadas))
    for f, g, t in cruces(G, caras, plegadas, iso, dic, idx):
        assert (min(f, t), max(f, t)) in pares
        assert (min(g, t), max(g, t)) in pares
        assert t != f and t != g


# --- colineales ------------------------------------------------------------

@pytest.mark.parametrize("nombre,n_tt", [
    ("single", 0), ("map13", 0), ("bomb", 2), ("bomb8", 12),
])
def test_numero_de_taco_taco(nombre, n_tt):
    G, caras, plegadas, iso, dic, idx = monta(nombre)
    tt, _, _ = colineales(G, caras, plegadas, iso, dic, idx)
    assert len(tt) == n_tt


@pytest.mark.parametrize("nombre", ["single", "map13", "bomb", "bomb8"])
def test_sin_aristas_planas_no_hay_tortillas(nombre):
    """Ninguna instancia de juguete tiene aristas F."""
    G, caras, plegadas, iso, dic, idx = monta(nombre)
    _, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)
    assert tto == [] and oo == []


def test_taco_taco_usa_cuatro_caras_distintas():
    G, caras, plegadas, iso, dic, idx = monta("bomb8")
    pares = set(pares_comparables(plegadas))
    tt, _, _ = colineales(G, caras, plegadas, iso, dic, idx)
    for cuatro in tt:
        caras_usadas = {c for par in cuatro for c in par}
        assert len(caras_usadas) == 4
        for u, v in cuatro:
            assert (min(u, v), max(u, v)) in pares
            
            
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

