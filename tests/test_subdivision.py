# -*- coding: utf-8 -*-
"""
Tests de la descomposición trapezoidal del Lema 6.4.

Estas tres funciones operan sobre un polígono suelto, sin grafo ni
FoldInput, así que se prueban aisladas.
"""
import pytest
from pathlib import Path
from flatfold.subdivision import (franjas, aristas_en_franja,
                                  trapecios_de_cara, cortes_de_vertical, 
                                  inserta_cortes, segmentos_verticales_nuevos,
                                  necesita_subdivision, subdivide_vertical)
from flatfold.geometry import es_convexa, doble_area_poligono
from flatfold.io_fold import read_fold, Kind
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas)
from flatfold.model import analiza_completo



TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"
# hexágono en L, en coordenadas escaladas (pares)
ELE = [(0, 0), (4, 0), (4, 2), (2, 2), (2, 4), (0, 4)]
CUADRADO = [(0, 0), (4, 0), (4, 4), (0, 4)]
TRIANGULO = [(0, 0), (4, 0), (4, 4)]
# una U: dos brazos hacia arriba
U = [(0, 0), (6, 0), (6, 4), (4, 4), (4, 2), (2, 2), (2, 4), (0, 4)]


# ===== franjas =============================================================

def test_franjas_ele():
    assert franjas(ELE) == [(0, 2), (2, 4)]


def test_franjas_cuadrado():
    assert franjas(CUADRADO) == [(0, 4)]


def test_franjas_u():
    assert franjas(U) == [(0, 2), (2, 4), (4, 6)]


# ===== aristas_en_franja ===================================================

def test_numero_par_de_aristas():
    """El lema lo demuestra: cada cruce alterna dentro/fuera."""
    for p in (ELE, CUADRADO, TRIANGULO, U):
        for a, b in franjas(p):
            assert len(aristas_en_franja(p, a, b)) % 2 == 0


def test_orden_de_abajo_arriba():
    for p in (ELE, U):
        for a, b in franjas(p):
            cr = aristas_en_franja(p, a, b)
            ys = [(ya, yb) for _, ya, _, yb in cr]
            assert ys == sorted(ys)


def test_franja_del_cuadrado():
    """Dos aristas: la base y la tapa."""
    cr = aristas_en_franja(CUADRADO, 0, 4)
    assert cr == [(0, 0, 4, 0), (0, 4, 4, 4)]


def test_verticales_se_descartan():
    """Los lados verticales no cruzan ninguna franja."""
    cr = aristas_en_franja(CUADRADO, 0, 4)
    assert all(x1 != x2 or ya != yb for x1, ya, x2, yb in cr)
    assert len(cr) == 2          # y no 4


# ===== trapecios_de_cara ===================================================

@pytest.mark.parametrize("p", [ELE, CUADRADO, TRIANGULO, U])
def test_todas_las_piezas_son_convexas(p):
    for pieza in trapecios_de_cara(p):
        assert es_convexa(pieza), pieza


@pytest.mark.parametrize("p", [ELE, CUADRADO, TRIANGULO, U])
def test_se_conserva_el_area(p):
    total = sum(doble_area_poligono(q) for q in trapecios_de_cara(p))
    assert total == doble_area_poligono(p)


@pytest.mark.parametrize("p", [ELE, CUADRADO, TRIANGULO, U])
def test_orientacion_antihoraria(p):
    """Cada pieza sale con área positiva."""
    for pieza in trapecios_de_cara(p):
        assert doble_area_poligono(pieza) > 0


@pytest.mark.parametrize("p", [ELE, CUADRADO, TRIANGULO, U])
def test_vertices_en_la_reticula(p):
    """Coordenadas pares: la exactitud posterior depende de ello."""
    for pieza in trapecios_de_cara(p):
        for x, y in pieza:
            assert x % 2 == 0 and y % 2 == 0


def test_convexo_no_se_parte():
    """Un cuadrado ya es convexo: una sola pieza."""
    assert len(trapecios_de_cara(CUADRADO)) == 1


def test_ele_da_dos_piezas():
    piezas = trapecios_de_cara(ELE)
    assert len(piezas) == 2
    areas = sorted(doble_area_poligono(q) for q in piezas)
    assert areas == [8, 16]        # rectángulos 2x2 y 2x4


def test_u_da_tres_piezas():
    piezas = trapecios_de_cara(U)
    assert len(piezas) == 3
    assert sum(doble_area_poligono(q) for q in piezas) == \
           doble_area_poligono(U)


def test_triangulo_no_se_parte():
    """Ya es convexo y tiene una sola franja."""
    piezas = trapecios_de_cara(TRIANGULO)
    assert len(piezas) == 1
    assert doble_area_poligono(piezas[0]) == doble_area_poligono(TRIANGULO)




CUAD_V = [(0, 0), (4, 0), (4, 4), (0, 4)]
CUAD_A = [(0, 1), (1, 2), (2, 3), (3, 0)]


def test_corta_base_y_tapa():
    """La vertical x=2 atraviesa la base y la tapa, no los lados."""
    c = cortes_de_vertical(CUAD_V, CUAD_A, 2)
    assert c == {0: (2, 0), 2: (2, 4)}


def test_no_corta_en_los_extremos():
    """x=0 pasa por dos vértices: no hay nada que insertar."""
    assert cortes_de_vertical(CUAD_V, CUAD_A, 0) == {}


def test_fuera_del_poligono():
    assert cortes_de_vertical(CUAD_V, CUAD_A, 8) == {}


def test_corta_una_diagonal():
    v = [(0, 0), (4, 4)]
    assert cortes_de_vertical(v, [(0, 1)], 2) == {0: (2, 2)}


def test_corta_una_antidiagonal():
    v = [(0, 4), (4, 0)]
    assert cortes_de_vertical(v, [(0, 1)], 2) == {0: (2, 2)}


def test_vertical_sobre_si_misma():
    """Una arista vertical en x=2 no se parte por la recta x=2."""
    v = [(2, 0), (2, 4)]
    assert cortes_de_vertical(v, [(0, 1)], 2) == {}


V = [(0, 0), (4, 0), (4, 4), (0, 4)]
A = [(0, 1), (1, 2), (2, 3), (3, 0)]
T = [Kind.BOUNDARY] * 4


def test_una_arista_se_parte_en_dos():
    c = {0: (2, 0)}
    v2, a2, t2, o2 = inserta_cortes(V, A, T, c)
    assert len(v2) == 5
    assert len(a2) == 5
    assert (0, 4) in a2 and (4, 1) in a2      # el nuevo vértice es el 4


def test_el_tipo_se_hereda():
    tipos = [Kind.CREASE] + [Kind.BOUNDARY] * 3
    _, a2, t2, _ = inserta_cortes(V, A, tipos, {0: (2, 0)})
    creases = [t for t in t2 if t == Kind.CREASE]
    assert len(creases) == 2


def test_el_origen_se_hereda():
    _, a2, _, o2 = inserta_cortes(V, A, T, {0: (2, 0)})
    assert sum(1 for j in o2 if o2[j] == 0) == 2


def test_deduplica_puntos():
    """Dos aristas cortadas en el mismo punto comparten vértice."""
    v = [(0, 0), (4, 0), (0, 4), (4, 4)]
    a = [(0, 1), (2, 3)]
    # cortes distintos, mismo punto imposible aquí; probamos contra
    # un vértice existente:
    v2, a2, _, _ = inserta_cortes(v, a, [Kind.BOUNDARY] * 2, {})
    assert len(v2) == 4


def test_sin_cortes_no_cambia_nada():
    v2, a2, t2, o2 = inserta_cortes(V, A, T, {})
    assert v2 == V and a2 == A and t2 == T




def test_cuadrado_no_genera_nada():
    """Ya es convexo: una pieza, y sus lados son la frontera."""
    v = [(0, 0), (4, 0), (4, 4), (0, 4)]
    a = [(0, 1), (1, 2), (2, 3), (3, 0)]
    piezas = trapecios_de_cara(v)
    assert segmentos_verticales_nuevos(piezas, v, a) == []


def test_ele_genera_el_corte():
    """El hexágono en L se parte por x=2, entre y=0 e y=2."""
    v = [(0, 0), (4, 0), (4, 2), (2, 2), (2, 4), (0, 4)]
    a = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    piezas = trapecios_de_cara(v)
    nuevos = segmentos_verticales_nuevos(piezas, v, a)
    assert nuevos == [((2, 0), (2, 2))]
    
# def test_subdivision_deja_todo_convexo():
#     fi = read_fold(TOY / "dos_eles.fold")
#     G, caras, ext = ...
#     assert necesita_subdivision(G, caras, ext)      # hay no convexas
#     fi2, origen = subdivide_vertical(fi, G, caras, ext)
#     G2, caras2, ext2 = ...                          # reconstruir
#     assert necesita_subdivision(G2, caras2, ext2) == []

def topologia(fi):
    """Grafo, caras e índice de la exterior a partir de un FoldInput."""
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, exterior = caras_definitivas(G, caras)
    return G, caras, caras.index(exterior)


def test_subdivision_deja_todo_convexo():
    fi = read_fold(TOY / "dos_eles.fold")
    G, caras, ext = topologia(fi)
    assert necesita_subdivision(G, caras, ext)          # hay no convexas

    fi2, origen = subdivide_vertical(fi, G, caras, ext)
    G2, caras2, ext2 = topologia(fi2)
    assert necesita_subdivision(G2, caras2, ext2) == []


def test_subdivision_conserva_el_area():
    from flatfold.pattern import area_cara
    fi = read_fold(TOY / "dos_eles.fold")
    G, caras, ext = topologia(fi)
    area = sum(area_cara(G, c) for i, c in enumerate(caras) if i != ext)

    fi2, _ = subdivide_vertical(fi, G, caras, ext)
    G2, caras2, ext2 = topologia(fi2)
    area2 = sum(area_cara(G2, c) for i, c in enumerate(caras2) if i != ext2)
    assert area == area2


def test_subdivision_anade_aristas_planas():
    from flatfold.io_fold import Kind
    fi = read_fold(TOY / "dos_eles.fold")
    G, caras, ext = topologia(fi)
    fi2, _ = subdivide_vertical(fi, G, caras, ext)
    planas = [t for t in fi2.tipos if t == Kind.FLAT]
    assert len(planas) == 2          # los cortes en x=2 y x=6


def test_patron_convexo_no_se_toca():
    fi = read_fold(TOY / "waterbomb4.fold")
    G, caras, ext = topologia(fi)
    fi2, origen = subdivide_vertical(fi, G, caras, ext)
    assert fi2 is fi
    assert origen == {i: i for i in range(len(fi.aristas))}
    
    
    
    
def test_representantes_uno_por_pliegue():
    from flatfold.subdivision import representantes
    from flatfold.io_fold import Kind

    fi = read_fold(TOY / "dos_eles.fold")
    G, caras, ext = topologia(fi)
    fi2, origen = subdivide_vertical(fi, G, caras, ext)

    rep = representantes(origen, fi2.tipos, fi)
    n_pliegues = sum(1 for t in fi.tipos if t == Kind.CREASE)
    assert len(rep) == n_pliegues


def test_representantes_apuntan_a_pliegues():
    from flatfold.subdivision import representantes
    from flatfold.io_fold import Kind

    fi = read_fold(TOY / "dos_eles.fold")
    G, caras, ext = topologia(fi)
    fi2, origen = subdivide_vertical(fi, G, caras, ext)

    rep = representantes(origen, fi2.tipos, fi)
    for original, nueva in rep.items():
        assert fi.tipos[original] == Kind.CREASE
        assert fi2.tipos[nueva] == Kind.CREASE
        assert origen[nueva] == original


def test_representantes_sin_subdividir_es_la_identidad():
    """En un patrón convexo cada pliegue es su propio representante."""
    from flatfold.subdivision import representantes
    from flatfold.io_fold import Kind

    fi = read_fold(TOY / "waterbomb4.fold")
    G, caras, ext = topologia(fi)
    fi2, origen = subdivide_vertical(fi, G, caras, ext)

    rep = representantes(origen, fi2.tipos, fi)
    assert all(k == v for k, v in rep.items())
    assert set(rep) == {i for i, t in enumerate(fi.tipos)
                        if t == Kind.CREASE}


def test_los_cortes_no_son_representantes():
    """Las aristas planas que introduce la subdivisión no vienen de
    ningún pliegue del fichero."""
    from flatfold.subdivision import representantes

    fi = read_fold(TOY / "dos_eles.fold")
    G, caras, ext = topologia(fi)
    fi2, origen = subdivide_vertical(fi, G, caras, ext)

    rep = representantes(origen, fi2.tipos, fi)
    cortes = {i for i, o in origen.items() if o is None}
    assert cortes and not (cortes & set(rep.values()))
    
    
def test_la_subdivision_activa_tortilla_tortilla():
    """
    Los cortes de la subdivisión son aristas planas, y dos de ellas
    colineales generan la familia tortilla-tortilla, que ninguna
    instancia del conjunto público ejercita.
    """
    a = analiza_completo(TOY / "dos_eles.fold")
    assert a.subdividido
    assert len(a.oo) == 1

# !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
