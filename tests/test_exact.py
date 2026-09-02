# -*- coding: utf-8 -*-
from fractions import Fraction
import pytest
from flatfold.exact import (pendiente_admisible, reflexion, aplica_isometria,
                            compon_isometrias, es_identidad, 
                            normaliza_a_reticula)



def test_pendientes_admisibles():
    assert pendiente_admisible((0,0), (4,0))     # horizontal
    assert pendiente_admisible((0,0), (0,4))     # vertical
    assert pendiente_admisible((0,0), (4,4))     # diagonal
    assert pendiente_admisible((0,0), (4,-4))    # antidiagonal


def test_pendientes_rechazadas():
    assert not pendiente_admisible((0,0), (6,2))   # pendiente 1/3
    assert not pendiente_admisible((0,0), (2,6))   # pendiente 3
    assert not pendiente_admisible((0,0), (0,0))   # degenerada
    

def test_reflexion_horizontal():
    r = reflexion((0, 4), (6, 4))                  # recta y = 4
    assert aplica_isometria(r, (2, 0)) == (2, 8)
    assert aplica_isometria(r, (2, 4)) == (2, 4)   # sobre la recta, fijo


def test_no_depende_de_los_extremos():
    assert reflexion((0, 4), (6, 4)) == reflexion((2, 4), (10, 4))


def test_compon_equivale_a_aplicar_dos_veces():
    f = reflexion((0,4), (6,4))
    g = reflexion((2,0), (2,6))
    p = (10, 0)
    assert (aplica_isometria(compon_isometrias(f, g), p)
            == aplica_isometria(f, aplica_isometria(g, p)))


def test_reflexion_es_involucion():
    for p, q in [((0,4),(6,4)), ((2,0),(2,6)), ((0,0),(4,4)), ((0,4),(4,0))]:
        r = reflexion(p, q)
        assert es_identidad(compon_isometrias(r, r))


def test_compon_no_conmuta():
    f = reflexion((0,0), (6,0))     # y = 0
    g = reflexion((0,0), (4,4))     # y = x
    assert compon_isometrias(f, g) != compon_isometrias(g, f)
    
    
    
    
def test_reticula_limpia():
    coords = [(0, 0), (1, 0), (1, 1)]
    salida, k = normaliza_a_reticula(coords)
    assert k == 1
    assert salida == [(0, 0), (1, 0), (1, 1)]

def test_quita_el_ruido():
    sucio = [(0.9999999999999998, 5.68e-16), (0.05000000000000004, 0.15)]
    salida, k = normaliza_a_reticula(sucio)
    assert k == 20
    assert salida == [(Fraction(1), Fraction(0)),
                      (Fraction(1, 20), Fraction(3, 20))]

def test_rechaza_lo_que_no_es_reticula():
    """
    Un patrón con ángulos arbitrarios no cae en ninguna retícula racional
    razonable y debe rechazarse en lugar de forzarse a una.
    """
    from math import sqrt
    with pytest.raises(ValueError):
        normaliza_a_reticula([(0, 0), (sqrt(2) / 2, 0)])
        
        
def test_acepta_reticulas_de_papel():
    """Sextos y séptimos aparecen en patrones 3x3 y son legítimos."""
    salida, k = normaliza_a_reticula([(0, 0),
                                      (0.1666666666666667, 0),
                                      (0.8571428571428572, 0)])
    assert k == 42                       # lcm(6, 7)
    assert salida[1][0] == Fraction(1, 6)
    assert salida[2][0] == Fraction(6, 7)
    
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
