# -*- coding: utf-8 -*-
import pytest
from pathlib import Path
from flatfold.io_fold import read_fold, Kind

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"

SINGLE = TOY / "single_fold.fold"
HALF   = TOY / "single_fold_half.fold"
SQUARE = TOY / "square.fold"
MAP13  = TOY / "map_1x3.fold"


def test_lectura_pliegue_unico():
    r = read_fold(SINGLE)
    assert len(r.vertices) == 6
    assert len(r.aristas) == 7
    assert r.escala == 2
    assert r.tipos[0] == Kind.BOUNDARY
    assert r.tipos[6] == Kind.CREASE
    assert r.dicc_montana_valle == {6: 1}


def test_lectura_cuadrado():
    r = read_fold(SQUARE)
    assert len(r.vertices) == 4
    assert len(r.aristas) == 4
    assert all(k == Kind.BOUNDARY for k in r.tipos)
    assert r.dicc_montana_valle == {}


def test_lectura_map13():
    r = read_fold(MAP13)
    assert len(r.vertices) == 8
    assert len(r.aristas) == 10
    assert r.tipos[8] == Kind.CREASE
    assert r.tipos[9] == Kind.CREASE
    assert r.dicc_montana_valle == {8: 1, 9: 0}   # M y V


@pytest.mark.parametrize("ruta", [SINGLE, HALF, SQUARE, MAP13])
def test_coordenadas_pares(ruta):
    r = read_fold(ruta)
    assert all(x % 2 == 0 and y % 2 == 0 for x, y in r.vertices)


def test_escala_normaliza():
    r1, r2 = read_fold(SINGLE), read_fold(HALF)
    assert r2.escala == 4
    assert r1.vertices == r2.vertices
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
    