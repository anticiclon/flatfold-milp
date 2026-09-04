# -*- coding: utf-8 -*-
import pytest
from pathlib import Path
from flatfold.io_fold import read_fold, Kind
from flatfold.pattern import (crea_grafo_dirigido, vecinos_ordenados,
                              caras_grafo, caras_definitivas,
                              area_cara, aristas_caras, plegar, grafo_dual,
                              cara_plegada, verifica_consistencia)

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"

RUTAS = {
    "single": TOY / "single_fold.fold",
    "square": TOY / "square.fold",
    "map13":  TOY / "map_1x3.fold",
    "bomb":   TOY / "waterbomb4.fold",
    "bomb8":   TOY / "waterbomb8.fold",
}


def monta(nombre):
    """Cadena completa hasta el plegado."""
    fi, G, caras = construye(nombre)
    ext = indice_exterior(G, caras)
    dual = grafo_dual(G, caras, ext)
    ref = min(dual.nodes)
    iso = plegar(G, dual, ref)
    return fi, G, caras, dual, iso


def construye(nombre):
    """Cadena completa hasta las caras."""
    fi = read_fold(RUTAS[nombre])
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    return fi, G, caras


def indice_exterior(G, caras):
    _, exterior = caras_definitivas(G, caras)
    return caras.index(exterior)


# --- grafo -----------------------------------------------------------------

@pytest.mark.parametrize("nombre,n_vert,n_aristas", [
    ("single", 6,  7),
    ("square", 4,  4),
    ("map13",  8, 10),
])
def test_grafo_tiene_ambos_sentidos(nombre, n_vert, n_aristas):
    _, G, _ = construye(nombre)
    assert G.number_of_nodes() == n_vert
    assert G.number_of_edges() == 2 * n_aristas


def test_orden_angular_vertice_1():
    _, G, _ = construye("single")
    assert vecinos_ordenados(G)[1] == [2, 4, 0]   # Este, Norte, Oeste


# --- caras -----------------------------------------------------------------

@pytest.mark.parametrize("nombre,n_ciclos,n_interiores,area_papel", [
    ("single", 3, 2, 16),
    ("square", 2, 1,  8),
    ("map13",  4, 3, 24),
])
def test_caras(nombre, n_ciclos, n_interiores, area_papel):
    _, G, caras = construye(nombre)
    assert len(caras) == n_ciclos
    interiores, exterior = caras_definitivas(G, caras)
    assert len(interiores) == n_interiores
    assert sum(area_cara(G, c) for c in interiores) == area_papel
    assert area_cara(G, exterior) == -area_papel


@pytest.mark.parametrize("nombre", ["single", "square", "map13"])
def test_ciclos_particionan_las_aristas(nombre):
    _, G, caras = construye(nombre)
    todas = [h for cara in caras for h in cara]
    assert len(todas) == G.number_of_edges()
    assert len(set(todas)) == len(todas)


@pytest.mark.parametrize("nombre", ["single", "square", "map13"])
def test_euler(nombre):
    _, G, caras = construye(nombre)
    assert G.number_of_nodes() - G.number_of_edges() // 2 + len(caras) == 2


def test_caras_interiores_iguales_map13():
    _, G, caras = construye("map13")
    interiores, _ = caras_definitivas(G, caras)
    assert sorted(area_cara(G, c) for c in interiores) == [8, 8, 8]


# --- adyacencias -----------------------------------------------------------

@pytest.mark.parametrize("nombre,n_aristas,n_internas", [
    ("single",  7, 1),
    ("square",  4, 0),
    ("map13",  10, 2),
])
def test_aristas_interiores(nombre, n_aristas, n_internas):
    _, G, caras = construye(nombre)
    ext = indice_exterior(G, caras)
    dic = aristas_caras(G, caras)
    internas = [a for a, (f, g) in dic.items() if f != ext and g != ext]
    assert len(dic) == n_aristas
    assert len(internas) == n_internas


def test_arista_interior_es_el_pliegue():
    _, G, caras = construye("single")
    ext = indice_exterior(G, caras)
    dic = aristas_caras(G, caras)
    internas = [a for a, (f, g) in dic.items() if f != ext and g != ext]
    assert set(internas[0]) == {1, 4}


@pytest.mark.parametrize("nombre", ["single", "square", "map13"])
def test_borde_geometrico_coincide_con_el_fichero(nombre):
    fi, G, caras = construye(nombre)
    ext = indice_exterior(G, caras)
    dic = aristas_caras(G, caras)
    geom = {frozenset(a) for a, (f, g) in dic.items() if ext in (f, g)}
    fich = {frozenset(a) for a, k in zip(fi.aristas, fi.tipos)
            if k == Kind.BOUNDARY}
    assert geom == fich
    
    
def test_tipo_en_ambos_sentidos():
    _, G, _ = construye("single")
    assert G[1][4]["tipo"] == Kind.CREASE
    assert G[4][1]["tipo"] == Kind.CREASE
    assert G[0][1]["tipo"] == Kind.BOUNDARY
    
    
    
    
    
# --- dual ------------------------------------------------------------------

@pytest.mark.parametrize("nombre,n_nodos,n_aristas", [
    ("single", 2, 1),
    ("square", 0, 0),      # sin pliegues: el dual queda vacío
    ("map13",  3, 2),
    ("bomb",   4, 4),
])
def test_dual(nombre, n_nodos, n_aristas):
    _, G, caras = construye(nombre)
    dual = grafo_dual(G, caras, indice_exterior(G, caras))
    assert dual.number_of_nodes() == n_nodos
    assert dual.number_of_edges() == n_aristas


def test_dual_bomb_tiene_ciclo():
    """Cuatro caras y cuatro aristas: hay una arista fuera del árbol."""
    _, G, caras = construye("bomb")
    dual = grafo_dual(G, caras, indice_exterior(G, caras))
    assert dual.number_of_edges() > dual.number_of_nodes() - 1


# --- plegado ---------------------------------------------------------------

@pytest.mark.parametrize("nombre", ["single", "map13", "bomb"])
def test_todas_las_caras_reciben_isometria(nombre):
    _, _, _, dual, iso = monta(nombre)
    assert set(iso) == set(dual.nodes)


def test_plegado_pliegue_unico():
    _, G, caras, dual, iso = monta("single")
    otras = [f for f in iso if iso[f] != ((1,0,0,1), (0,0))]
    assert len(otras) == 1
    assert iso[otras[0]] == ((-1, 0, 0, 1), (4, 0))   # reflexión en x = 2


def test_plegado_bomb_isometrias():
    _, G, caras, dual, iso = monta("bomb")
    esperadas = {((1, 0, 0, 1), (0, 0)),        # identidad
                 ((-1, 0, 0, 1), (4, 0)),       # reflexión en x = 2
                 ((1, 0, 0, -1), (0, 4)),       # reflexión en y = 2
                 ((-1, 0, 0, -1), (4, 4))}      # giro de 180 en (2,2)
    assert set(iso.values()) == esperadas


def test_orientaciones_alternan_bomb():
    _, G, caras, dual, iso = monta("bomb")
    dets = sorted((a*d - b*c) for (a,b,c,d), _ in iso.values())
    assert dets == [-1, -1, 1, 1]


# --- caras plegadas --------------------------------------------------------

def test_caras_plegadas_coinciden_bomb():
    """Los cuatro triángulos caen exactamente uno sobre otro."""
    _, G, caras, dual, iso = monta("bomb")
    conjuntos = {frozenset(cara_plegada(G, caras[f], iso[f])) for f in iso}
    assert len(conjuntos) == 1


def test_plegado_conserva_area():
    _, G, caras, dual, iso = monta("bomb")
    for f in iso:
        original = abs(area_cara(G, caras[f]))
        plegada = cara_plegada(G, caras[f], iso[f])
        w = plegada[1:] + plegada[:1]
        doble = abs(sum(p[0]*q[1] - q[0]*p[1] for p, q in zip(plegada, w)))
        assert doble == original


# --- consistencia ----------------------------------------------------------

@pytest.mark.parametrize("nombre", ["single", "map13", "bomb"])
def test_consistencia(nombre):
    _, G, _, dual, iso = monta(nombre)
    verifica_consistencia(G, dual, iso)


def test_referencia_no_afecta_a_la_consistencia():
    """Plegar desde cualquier cara produce un plegado consistente."""
    _, G, caras, dual, _ = monta("bomb")
    for r in sorted(dual.nodes):
        verifica_consistencia(G, dual, plegar(G, dual, r))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
def test_dual_waterbomb8():
    _, G, caras = construye("bomb8")
    dual = grafo_dual(G, caras, indice_exterior(G, caras))
    assert dual.number_of_nodes() == 8
    assert dual.number_of_edges() == 8


def test_consistencia_con_diagonales():
    """La instancia donde el orden de composición importa de verdad."""
    _, G, _, dual, iso = monta("bomb8")
    verifica_consistencia(G, dual, iso)


def test_orientaciones_alternan_bomb8():
    _, G, caras, dual, iso = monta("bomb8")
    dets = sorted((a*d - b*c) for (a, b, c, d), _ in iso.values())
    assert dets == [-1]*4 + [1]*4



def test_caras_plegadas_coinciden_bomb8():
    """Los ocho triángulos caen exactamente uno sobre otro."""
    _, G, caras, dual, iso = monta("bomb8")
    conjuntos = {frozenset(cara_plegada(G, caras[f], iso[f])) for f in iso}
    assert len(conjuntos) == 1
    
    
    
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
