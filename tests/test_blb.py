# -*- coding: utf-8 -*-
"""
Tests de la condición big-little-big.

BLB es una restricción IMPLICADA (Corolario 11.4): añadirla no puede
eliminar ningún estado válido. Ese es el criterio que valida los signos,
tanto en el modelo directo como a través del presolve.
"""
# import pytest
from pathlib import Path
from flatfold.pattern import vecinos_ordenados
from flatfold.validate import (vertices_interiores, sectores, minimos_locales)
from flatfold.objective import etiquetas_mv, restricciones_blb
from flatfold.model import construye_modelo, añade_restricciones
from flatfold.presolve import (literales_a_ecuaciones, reduce_completo,
                               componentes, reparte_por_componentes)
from flatfold.model import resuelve_por_componentes, analiza_completo

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
GRIDS = RAIZ / "instances" / "flatfolder" / "grids"

# instancia real pequeña con mínimos locales: 104 variables, 8 BLB
PEQUEÑA = GRIDS / "004_brown_3x3_A3.fold"
ESTADOS_PEQUEÑA = 13122


# def monta(ruta):
#     """Cadena completa, devolviendo todo lo que hace falta para BLB."""
#     fi = read_fold(ruta)
#     G = crea_grafo_dirigido(fi)
#     caras = caras_grafo(G)
#     _, ext = caras_definitivas(G, caras)
#     idx = caras.index(ext)
#     dual = grafo_dual(G, caras, idx)
#     dic = aristas_caras(G, caras)
#     iso = plegar(G, dual, min(dual.nodes))
#     plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
#     pares = pares_comparables(plegadas)
#     trip = tripletas(plegadas, pares)
#     cr = cruces(G, caras, plegadas, iso, dic, idx)
#     tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)
#     etiq = etiquetas_mv(G, fi, iso, dic, idx)
#     blb = restricciones_blb(G, fi, iso, dic, idx, etiq)
#     return dict(fi=fi, G=G, iso=iso, dic=dic, idx=idx, pares=pares,
#                 trip=trip, cr=cr, tt=tt, tto=tto, oo=oo, etiq=etiq, blb=blb)




def monta(ruta):
    a = analiza_completo(ruta)
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    blb = restricciones_blb(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext, etiq)
    # return dict(fi=a.fi, G=a.G, iso=a.iso, dic=a.dic_aristas, idx=a.idx_ext,
    #             pares=a.pares, trip=a.trip, cr=a.cr, tt=a.tt, tto=a.tto,
    #             oo=a.oo, etiq=etiq, blb=blb)
    return dict(fi=a.fi, G=a.G, iso=a.iso, dic=a.dic_aristas, idx=a.idx_ext,
                pares=a.pares, trip=a.trip, cr=a.cr + a.cr_oo,
                tt=a.tt, tto=a.tto, oo=a.oo, etiq=etiq, blb=blb)







def cuenta_directo(d, con_blb, limite=50000):
    m, x = construye_modelo(d["pares"])
    indice = {par: i for i, par in enumerate(d["pares"])}
    añade_restricciones(m, x, indice, d["trip"], d["cr"], d["tt"],
                        d["tto"], d["oo"], blb=d["blb"], con_blb=con_blb)
    m.Params.OutputFlag = 0
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.optimize()
    assert m.SolCount < limite, "el pool se llenó"
    return m.SolCount


def cuenta_presolve(d, con_blb, limite=50000):
    blb = d["blb"] if con_blb else None
    ig, par, fij = literales_a_ecuaciones(
        d["pares"], d["cr"], d["tto"], d["oo"], d["tt"], blb=blb)
    mapeo, restantes, vivas, _ = reduce_completo(
        len(d["pares"]), ig, par, fij, d["trip"], d["pares"])
    comps = componentes(mapeo, restantes, d["trip"], d["pares"])
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    n, _, _ = resuelve_por_componentes(datos, modo="recuento",
                                       limite_pool=limite)
    return n


# ===== sectores ============================================================

def test_angulos_suman_una_vuelta():
    """Los sectores de un vértice cubren 360 grados: 8 octantes."""
    d = monta(PEQUEÑA)
    vecinos = vecinos_ordenados(d["G"])
    for v in vertices_interiores(d["G"]):
        secs = sectores(d["G"], vecinos, v)
        if secs:
            assert sum(a for a, _, _ in secs) == 8


def test_bomba_de_agua_sin_minimos():
    """Ocho sectores iguales de 45 grados: ningún mínimo local."""
    d = monta(TOY / "waterbomb8.fold")
    vecinos = vecinos_ordenados(d["G"])
    for v in vertices_interiores(d["G"]):
        secs = sectores(d["G"], vecinos, v)
        assert len(secs) == 8
        assert all(a == 1 for a, _, _ in secs)
        assert minimos_locales(secs) == []
    assert d["blb"] == []


def test_sectores_de_vertice_de_borde():
    """vertices_interiores ya los excluye; aquí solo que no revienta."""
    d = monta(TOY / "map_1x3.fold")
    assert vertices_interiores(d["G"]) == set()
    assert d["blb"] == []


# ===== minimos_locales, aislada ===========================================

def _s(*angulos):
    return [(a, None, None) for a in angulos]


def test_minimo_en_el_medio():
    assert minimos_locales(_s(3, 1, 4)) == [1]


def test_minimo_con_cierre_ciclico():
    """El primero se compara con el último."""
    assert minimos_locales(_s(1, 4, 3)) == [0]


def test_empates_no_cuentan():
    """Estrictamente menor: con <= se generarían restricciones inválidas."""
    assert minimos_locales(_s(2, 2, 4)) == []
    assert minimos_locales(_s(3, 2, 2, 3)) == []


def test_todos_iguales():
    assert minimos_locales(_s(2, 2, 2, 2)) == []


def test_menos_de_tres_sectores():
    assert minimos_locales(_s(4, 4)) == []
    assert minimos_locales(_s(8)) == []
    assert minimos_locales([]) == []


def test_dos_minimos():
    assert minimos_locales(_s(1, 3, 1, 3)) == [0, 2]


# ===== restricciones_blb ==================================================

def test_numero_de_restricciones():
    d = monta(PEQUEÑA)
    assert len(d["blb"]) == 8


def test_formato():
    """Cada elemento son dos entradas de etiquetas: ((par), signo)."""
    d = monta(PEQUEÑA)
    for (p1, s1), (p2, s2) in d["blb"]:
        assert len(p1) == 2 and len(p2) == 2
        assert s1 in (0, 1) and s2 in (0, 1)
        assert p1 != p2


def test_coincide_con_los_minimos():
    d = monta(PEQUEÑA)
    vecinos = vecinos_ordenados(d["G"])
    n = sum(len(minimos_locales(sectores(d["G"], vecinos, v)))
            for v in vertices_interiores(d["G"]))
    assert len(d["blb"]) == n


# ===== lo que de verdad valida los signos =================================

def test_recuento_base():
    d = monta(PEQUEÑA)
    assert cuenta_directo(d, con_blb=False) == ESTADOS_PEQUEÑA


def test_blb_no_cambia_el_recuento():
    """
    BLB es implicada: no puede eliminar estados. Si el recuento baja,
    hay un error de signo en restricciones_blb o en añade_blb.
    """
    d = monta(PEQUEÑA)
    assert cuenta_directo(d, con_blb=True) == ESTADOS_PEQUEÑA


def test_blb_no_cambia_el_recuento_con_presolve():
    """
    Lo mismo por el camino reducido. Valida blb_a_igualdades, donde el
    complemento del signo puede ir mal puesto.
    """
    d = monta(PEQUEÑA)
    assert cuenta_presolve(d, con_blb=False) == ESTADOS_PEQUEÑA
    assert cuenta_presolve(d, con_blb=True) == ESTADOS_PEQUEÑA


def test_blb_reduce_o_iguala_variables():
    """BLB alimenta al presolve: no puede aumentar las supervivientes."""
    d = monta(PEQUEÑA)

    def supervivientes(blb):
        ig, par, fij = literales_a_ecuaciones(
            d["pares"], d["cr"], d["tto"], d["oo"], d["tt"], blb=blb)
        mapeo, _, _, _ = reduce_completo(
            len(d["pares"]), ig, par, fij, d["trip"], d["pares"])
        return len({e[1] for e in mapeo.values() if e[0] == "var"})

    assert supervivientes(d["blb"]) <= supervivientes(None)
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

    
    
    
    
    
    
    
    
    
    