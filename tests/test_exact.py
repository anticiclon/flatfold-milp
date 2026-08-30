# -*- coding: utf-8 -*-
from flatfold.exact import (pendiente_admisible, reflexion, aplica_isometria,
                            compon_isometrias, es_identidad)


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
    
    
    
    
    
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
