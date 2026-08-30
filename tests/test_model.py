# -*- coding: utf-8 -*-
import pytest
from pathlib import Path
import gurobipy as gb

from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada)
from flatfold.incidences import pares_comparables
from flatfold.model import (construye_modelo, literal, resuelve,
                            orden_de_capas)

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"

RUTAS = {
    "single": TOY / "single_fold.fold",
    "square": TOY / "square.fold",
    "map13":  TOY / "map_1x3.fold",
    "bomb":   TOY / "waterbomb4.fold",
    "bomb8":  TOY / "waterbomb8.fold",
}


def cadena(nombre):
    """Tubería completa hasta los pares comparables."""
    fi = read_fold(RUTAS[nombre])
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, exterior = caras_definitivas(G, caras)
    dual = grafo_dual(G, caras, caras.index(exterior))
    if len(dual) == 0:
        return []
    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
    return pares_comparables(plegadas)


def cuenta_soluciones(pares, limite=1000):
    m, _ = construye_modelo(pares)
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.Params.LogToConsole = 0
    m.optimize()
    return m.SolCount


# --- construcción del modelo ----------------------------------------------

@pytest.mark.parametrize("nombre,n_vars", [
    ("single", 1),
    ("square", 0),
    ("map13",  3),
])
def test_numero_de_variables(nombre, n_vars):
    m, x = construye_modelo(cadena(nombre))
    assert len(x) == n_vars
    assert m.NumVars == n_vars


def test_modelo_sin_variables_es_factible():
    """El cuadrado no tiene pliegues: nada que decidir, pero no debe fallar."""
    m, x = construye_modelo([])
    estado, valores = resuelve(m, x)
    assert valores == {}
    assert orden_de_capas([], valores) == []


# --- literal ---------------------------------------------------------------

def test_literal_devuelve_la_variable():
    m, x = construye_modelo([(0, 2)])
    assert literal(x, 0, 2) is x[0, 2]


def test_literal_antisimetrico():
    m, x = construye_modelo([(0, 2)])
    suma = literal(x, 0, 2) + literal(x, 2, 0)
    assert suma.getConstant() == 1
    assert suma.size() == 2                       # x  y  -x
    assert suma.getCoeff(0) + suma.getCoeff(1) == 0


def test_literal_sobre_la_misma_cara_falla():
    m, x = construye_modelo([(0, 2)])
    with pytest.raises(ValueError):
        literal(x, 1, 1)


def test_literal_de_par_inexistente_falla():
    """Si dos caras no se solapan no hay variable: debe saltar, no inventar."""
    m, x = construye_modelo([(0, 2)])
    with pytest.raises(KeyError):
        literal(x, 0, 5)


# --- orden de capas --------------------------------------------------------

def test_orden_de_capas_dos_caras():
    assert orden_de_capas([(0, 2)], {(0, 2): 1}) == [[0, 2]]
    assert orden_de_capas([(0, 2)], {(0, 2): 0}) == [[2, 0]]


def test_orden_de_capas_tres_caras():
    """0 encima de 1, 1 encima de 2: la pila es 0, 1, 2."""
    pares = [(0, 1), (0, 2), (1, 2)]
    valores = {(0, 1): 1, (0, 2): 1, (1, 2): 1}
    assert orden_de_capas(pares, valores) == [[0, 1, 2]]


def test_orden_de_capas_invertido():
    pares = [(0, 1), (0, 2), (1, 2)]
    valores = {(0, 1): 0, (0, 2): 0, (1, 2): 0}
    assert orden_de_capas(pares, valores) == [[2, 1, 0]]


def test_orden_de_capas_dos_grupos():
    """Caras 0-1 en un montón y 2-3 en otro: dos grupos independientes."""
    pares = [(0, 1), (2, 3)]
    valores = {(0, 1): 1, (2, 3): 0}
    assert orden_de_capas(pares, valores) == [[0, 1], [3, 2]]


# --- la cadena completa ----------------------------------------------------

def test_pliegue_unico_tiene_dos_estados():
    """
    Un papel doblado por la mitad se pliega de dos formas: monte y valle.
    Es el primer resultado verificable de toda la tubería.
    """
    assert cuenta_soluciones(cadena("single")) == 2


def test_cuadrado_tiene_un_estado():
    assert cuenta_soluciones(cadena("square")) == 1


@pytest.mark.parametrize("nombre", ["single", "square", "map13",
                                    "bomb", "bomb8"])
def test_la_cadena_produce_solucion(nombre):
    """De fichero FOLD a orden de capas, sin reventar."""
    pares = cadena(nombre)
    m, x = construye_modelo(pares)
    estado, valores = resuelve(m, x)
    assert valores is not None
    capas = orden_de_capas(pares, valores)
    caras_en_capas = sum(len(g) for g in capas)
    caras_con_variable = len({c for p in pares for c in p})
    assert caras_en_capas == caras_con_variable
    
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
