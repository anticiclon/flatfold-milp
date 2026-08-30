# -*- coding: utf-8 -*-
import pytest
from pathlib import Path
from flatfold.io_fold import read_fold
from flatfold.pattern import crea_grafo_dirigido, vecinos_ordenados
from flatfold.exact import es_identidad
from flatfold.validate import (valida_pendientes, vertices_interiores,
                               holonomia, valida_kawasaki)

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"

RUTAS = {
    "single": TOY / "single_fold.fold",
    "square": TOY / "square.fold",
    "map13":  TOY / "map_1x3.fold",
    "bomb":   TOY / "waterbomb4.fold",
}


def grafo(nombre):
    return crea_grafo_dirigido(read_fold(RUTAS[nombre]))


# --- pendientes ------------------------------------------------------------

@pytest.mark.parametrize("nombre", ["single", "square", "map13", "bomb"])
def test_pendientes_ok(nombre):
    valida_pendientes(read_fold(RUTAS[nombre]))


def test_pendiente_ilegal_se_rechaza():
    fi = read_fold(TOY / "bad_slope.fold")
    with pytest.raises(ValueError):
        valida_pendientes(fi)


# --- vértices interiores ---------------------------------------------------

@pytest.mark.parametrize("nombre", ["single", "square", "map13"])
def test_sin_vertices_interiores(nombre):
    assert vertices_interiores(grafo(nombre)) == set()


def test_vertice_interior_waterbomb():
    assert vertices_interiores(grafo("bomb")) == {4}


# --- Kawasaki --------------------------------------------------------------

def test_holonomia_es_identidad_waterbomb():
    G = grafo("bomb")
    assert es_identidad(holonomia(G, vecinos_ordenados(G), 4))


@pytest.mark.parametrize("nombre", ["single", "square", "map13", "bomb"])
def test_kawasaki_ok(nombre):
    valida_kawasaki(grafo(nombre))


def test_kawasaki_se_rechaza():
    G = crea_grafo_dirigido(read_fold(TOY / "bad_kawasaki.fold"))
    with pytest.raises(ValueError):
        valida_kawasaki(G)
        
        
        
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
