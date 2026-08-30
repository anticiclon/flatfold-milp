# -*- coding: utf-8 -*-
import pytest
from flatfold.geometry import (es_convexa, caja, interseca_cajas, recorta,
                               poligono_de_caja, doble_area_poligono,
                               tiene_interior, se_solapan)

# Polígonos de referencia, todos con las cuatro direcciones admisibles
CUADRADO  = [(0, 0), (4, 0), (4, 4), (0, 4)]
TRIANGULO = [(0, 0), (4, 0), (4, 4)]
OCTOGONO  = [(2, 0), (6, 0), (8, 2), (8, 6),
             (6, 8), (2, 8), (0, 6), (0, 2)]
ELE       = [(0, 0), (4, 0), (4, 2), (2, 2), (2, 4), (0, 4)]   # no convexa


# --- convexidad ------------------------------------------------------------

@pytest.mark.parametrize("p", [CUADRADO, TRIANGULO, OCTOGONO])
def test_convexos(p):
    assert es_convexa(p)


def test_no_convexo():
    assert not es_convexa(ELE)


def test_menos_de_tres_vertices():
    assert not es_convexa([(0, 0), (2, 2)])


# --- caja ------------------------------------------------------------------

def test_caja_triangulo():
    # x en [0,4], y en [0,4], x+y en [0,8], x-y en [0,4]
    assert caja(TRIANGULO) == (0, 4, 0, 4, 0, 8, 0, 4)


def test_interseca_cajas_toma_las_cotas_mas_ajustadas():
    a = (0, 10, 0, 10, 0, 20, -10, 10)
    b = (2,  8, 3,  9, 5, 15,  -5,  5)
    assert interseca_cajas([a, b]) == (2, 8, 3, 9, 5, 15, -5, 5)


# --- recorte ---------------------------------------------------------------

def test_recorte_produce_el_triangulo():
    """El cuadrado recortado por x - y >= 0 es el triángulo inferior."""
    r = recorta(CUADRADO, "d", 0, -1)
    assert set(r) == set(TRIANGULO)
    assert doble_area_poligono(r) == doble_area_poligono(TRIANGULO)


def test_recorte_que_no_toca_devuelve_el_poligono():
    r = recorta(CUADRADO, "s", 100, +1)      # x+y <= 100, no recorta nada
    assert set(r) == set(CUADRADO)


def test_recorte_que_deja_fuera_devuelve_vacio():
    assert recorta(CUADRADO, "s", -100, +1) == []   # x+y <= -100


# --- ida y vuelta: la representación por ocho números es exacta ------------

@pytest.mark.parametrize("p", [CUADRADO, TRIANGULO, OCTOGONO])
def test_caja_reproduce_el_poligono(p):
    """Ocho números describen exactamente el polígono, no una envolvente."""
    q = poligono_de_caja(caja(p))
    assert set(q) == set(p)
    assert doble_area_poligono(q) == doble_area_poligono(p)


# --- interior --------------------------------------------------------------

@pytest.mark.parametrize("p", [CUADRADO, TRIANGULO, OCTOGONO])
def test_tiene_interior(p):
    assert tiene_interior(caja(p))


def test_caja_inconsistente():
    """
    x,y en [0,10] pero x+y <= 1 y x-y >= 5: cada intervalo por separado es
    razonable, pero el sistema no tiene solución. Es el caso que justifica
    convertir a polígono en lugar de mirar los ocho números sueltos.
    """
    assert not tiene_interior((0, 10, 0, 10, 0, 1, 5, 20))


def test_caja_degenerada_es_un_segmento():
    """x fijo: la región es un segmento vertical, sin interior."""
    assert not tiene_interior((2, 2, 0, 4, 2, 6, -2, 2))


# --- solapamiento ----------------------------------------------------------

def test_solapan_de_verdad():
    a = [(0, 0), (4, 0), (4, 4), (0, 4)]
    b = [(2, 2), (6, 2), (6, 6), (2, 6)]
    assert se_solapan(a, b)


def test_pegados_por_un_lado_no_solapan():
    """Compartir una arista no es solaparse: no hay área común."""
    a = [(0, 0), (4, 0), (4, 4), (0, 4)]
    b = [(4, 0), (8, 0), (8, 4), (4, 4)]
    assert not se_solapan(a, b)


def test_tocandose_en_un_vertice_no_solapan():
    a = [(0, 0), (4, 0), (4, 4), (0, 4)]
    b = [(4, 4), (8, 4), (8, 8), (4, 8)]
    assert not se_solapan(a, b)


def test_disjuntos_no_solapan():
    a = [(0, 0), (4, 0), (4, 4), (0, 4)]
    b = [(6, 0), (10, 0), (10, 4), (6, 4)]
    assert not se_solapan(a, b)


def test_uno_dentro_de_otro_solapan():
    a = [(0, 0), (8, 0), (8, 8), (0, 8)]
    b = [(2, 2), (4, 2), (4, 4), (2, 4)]
    assert se_solapan(a, b)
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

    
