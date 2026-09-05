# -*- coding: utf-8 -*-
# import pytest
from pathlib import Path
import gurobipy as gb

from flatfold.io_fold import Kind
from flatfold.model import (construye_modelo, 
                            añade_restricciones, 
                            analiza_completo)
from flatfold.objective import (etiquetas_mv, fija_etiquetas, energia,
                                etiquetas_por_pliegue_original)
from flatfold.exact import sigma


TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"


# def prepara(nombre):
#     """Modelo completo con restricciones, más las etiquetas M/V."""
#     fi = read_fold(TOY / f"{nombre}.fold")
#     G = crea_grafo_dirigido(fi)
#     caras = caras_grafo(G)
#     _, exterior = caras_definitivas(G, caras)
#     idx = caras.index(exterior)
#     dual = grafo_dual(G, caras, idx)
#     dic = aristas_caras(G, caras)
#     iso = plegar(G, dual, min(dual.nodes))
#     plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}

#     pares = pares_comparables(plegadas)
#     trip = tripletas(plegadas, pares)
#     cr = cruces(G, caras, plegadas, iso, dic, idx)
#     tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)

#     m, x = construye_modelo(pares)
#     # añade_restricciones(m, x, trip, cr, tt, tto, oo)
#     indice_pares = {par: i for i, par in enumerate(pares)}
#     añade_restricciones(m, x, indice_pares, trip, cr, tt, tto, oo)
#     m.Params.OutputFlag = 0

#     etiq = etiquetas_mv(G, fi, iso, dic, idx)
#     return fi, G, iso, dic, idx, m, x, etiq

def prepara(nombre):
    """Modelo completo con restricciones, más las etiquetas M/V."""
    a = analiza_completo(TOY / f"{nombre}.fold")
    m, x = construye_modelo(a.pares)
    indice_pares = {par: i for i, par in enumerate(a.pares)}
    añade_restricciones(m, x, indice_pares, a.trip, a.cr, a.tt, a.tto, a.oo)
    m.Params.OutputFlag = 0
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    return a.fi, a.G, a.iso, a.dic_aristas, a.idx_ext, m, x, etiq


# --- etiquetas -------------------------------------------------------------

def test_una_etiqueta_por_pliegue_interior():
    fi, G, iso, dic, idx, m, x, etiq = prepara("map_1x3")
    interiores = [a for a, (f, g) in dic.items() if idx not in (f, g)]
    creases = [a for a in interiores if G[a[0]][a[1]]["tipo"] == Kind.CREASE]
    assert len(etiq) == len(creases)


def test_etiqueta_usa_par_ordenado():
    _, _, _, _, _, _, _, etiq = prepara("waterbomb4")
    for (u, v), s in etiq.values():
        assert u < v
        assert s in (0, 1)


def test_etiqueta_no_depende_de_la_cara_elegida():
    """
    m_e = X_fg si sigma(f)=+1, y 1-X_fg si sigma(f)=-1. Como las dos caras
    de un pliegue tienen sigma opuesto y X_gf = 1-X_fg, el resultado no
    depende de cuál se tome como f.
    """
    fi, G, iso, dic, idx, m, x, etiq = prepara("waterbomb4")
    for arista, ((u, v), s) in etiq.items():
        f, g = dic[arista] if arista in dic else (None, None)
        # las dos caras del pliegue tienen orientación opuesta
        assert sigma(iso[u]) == -sigma(iso[v])
        # el signo guardado corresponde a la cara de índice menor
        assert s == (0 if sigma(iso[u]) == 1 else 1)


# --- asignación prescrita --------------------------------------------------

def test_asignacion_del_fichero_es_factible():
    """La M/V que trae el fichero debe ser realizable."""
    for nombre in ["map_1x3", "waterbomb4", "003_brown_3x3_A2"]:
        fi, _, _, _, _, m, x, etiq = prepara(nombre)
        n = fija_etiquetas(m, x, etiq, fi.dicc_montana_valle)
        assert n == len(fi.dicc_montana_valle)
        m.optimize()
        assert m.SolCount > 0, nombre
        
        
        
def test_asignacion_irrealizable_es_infactible():
    """
    waterbomb8 es geométricamente válido —cumple Kawasaki y se pliega—
    pero su asignación M,V,M,V,M,V,M,V viola Maekawa: en un vértice de
    grado 8 la diferencia entre montes y valles debe ser 2, no 0.

    El modelo no impone Maekawa; la detecta como consecuencia
    (Corolario 11.4). Este test comprueba justamente eso.
    """
    fi, _, _, _, _, m, x, etiq = prepara("waterbomb8")
    n = fija_etiquetas(m, x, etiq, fi.dicc_montana_valle)
    assert n == 8
    m.optimize()
    assert m.SolCount == 0
    assert m.Status == gb.GRB.INFEASIBLE


def test_waterbomb8_es_factible_sin_prescribir():
    """Sin fijar etiquetas sí existe orden de capas: el patrón se pliega,
    lo que falla es esa asignación concreta."""
    _, _, _, _, _, m, x, _ = prepara("waterbomb8")
    m.optimize()
    assert m.SolCount > 0


def test_reparacion_de_waterbomb8():
    """Cuántos pliegues hay que invertir para hacerla realizable."""
    fi, _, _, _, _, m, x, etiq = prepara("waterbomb8")
    m.setObjective(energia(x, etiq, fi.dicc_montana_valle), gb.GRB.MINIMIZE)
    m.optimize()
    assert m.SolCount > 0
    assert int(round(m.ObjVal)) >= 1
        
        
        
        
        
        
        
        
        


def test_la_solucion_respeta_las_etiquetas_fijadas():
    fi, _, _, _, _, m, x, etiq = prepara("map_1x3")
    fija_etiquetas(m, x, etiq, fi.dicc_montana_valle)
    m.optimize()
    assert m.SolCount > 0
    for arista, pedido in fi.dicc_montana_valle.items():
        (u, v), s = etiq[arista]
        xv = int(round(x[u, v].X))
        assert (xv if s == 0 else 1 - xv) == pedido


# --- energía ---------------------------------------------------------------

def test_energia_nula_si_la_asignacion_es_realizable():
    for nombre in ["map_1x3", "waterbomb4"]:
        fi, _, _, _, _, m, x, etiq = prepara(nombre)
        m.setObjective(energia(x, etiq, fi.dicc_montana_valle),
                       gb.GRB.MINIMIZE)
        m.optimize()
        assert m.SolCount > 0
        assert int(round(m.ObjVal)) == 0, nombre


def test_energia_cuenta_las_desviaciones():
    """Con la asignación invertida, el óptimo es el número de pliegues
    que no se pueden satisfacer a la vez; nunca negativo ni mayor que
    el número de pliegues."""
    fi, _, _, _, _, m, x, etiq = prepara("waterbomb4")
    invertida = {a: 1 - v for a, v in fi.dicc_montana_valle.items()}
    m.setObjective(energia(x, etiq, invertida), gb.GRB.MINIMIZE)
    m.optimize()
    assert m.SolCount > 0
    assert 0 <= int(round(m.ObjVal)) <= len(invertida)


def test_pesos_nulos_no_penalizan():
    fi, _, _, _, _, m, x, etiq = prepara("waterbomb4")
    invertida = {a: 1 - v for a, v in fi.dicc_montana_valle.items()}
    ceros = {a: 0 for a in invertida}
    m.setObjective(energia(x, etiq, invertida, ceros), gb.GRB.MINIMIZE)
    m.optimize()
    assert int(round(m.ObjVal)) == 0
    
# def test_energia_sobre_patron_subdividido():
#     a = analiza_completo(TOY / "dos_eles.fold")
#     etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
#     etiq_orig = etiquetas_por_pliegue_original(etiq, a.representantes)

#     assert set(etiq_orig) == set(a.fi_original.dicc_montana_valle)  # ojo
    
def test_etiquetas_por_pliegue_original():
    """
    Tras subdividir, las etiquetas se indexan por arista del patrón
    nuevo; la energía las necesita por arista del fichero.
    """
    a = analiza_completo(TOY / "dos_eles.fold")
    assert a.subdividido

    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)

    assert set(orig) == set(a.fi_original.dicc_montana_valle)


def test_sin_subdividir_no_cambia_nada():
    """En un patrón convexo la traducción es la identidad."""
    from flatfold.objective import etiquetas_por_pliegue_original
    a = analiza_completo(TOY / "waterbomb4.fold")
    assert not a.subdividido

    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    assert orig == etiq
    
    
def test_energia_sobre_patron_subdividido():
    """
    La cadena completa con subdivisión: el objetivo se expresa sobre los
    pliegues del fichero aunque el patrón se haya partido. Con un solo
    pliegue y su asignación realizable, el óptimo es 0.
    """
    from flatfold.objective import etiquetas_por_pliegue_original
    a = analiza_completo(TOY / "dos_eles.fold")
    assert a.subdividido

    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr, a.tt, a.tto, a.oo)
    m.Params.OutputFlag = 0

    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    m.setObjective(energia(x, orig, a.fi_original.dicc_montana_valle),
                   gb.GRB.MINIMIZE)
    m.optimize()

    assert m.SolCount > 0
    assert int(round(m.ObjVal)) == 0


def test_asignacion_prescrita_sobre_patron_subdividido():
    """Fijar la etiqueta del pliegue original debe seguir siendo factible."""
    from flatfold.objective import etiquetas_por_pliegue_original
    a = analiza_completo(TOY / "dos_eles.fold")

    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr, a.tt, a.tto, a.oo)
    m.Params.OutputFlag = 0

    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    n = fija_etiquetas(m, x, orig, a.fi_original.dicc_montana_valle)
    assert n == 1
    m.optimize()
    assert m.SolCount > 0


def test_invertir_el_pliegue_tambien_es_realizable():
    """
    Con un solo pliegue, las dos etiquetas son posibles: el papel se
    dobla hacia un lado o hacia el otro. Así que la energía con la
    asignación invertida también es 0.
    """
    from flatfold.objective import etiquetas_por_pliegue_original
    a = analiza_completo(TOY / "dos_eles.fold")

    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr, a.tt, a.tto, a.oo)
    m.Params.OutputFlag = 0

    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    invertida = {k: 1 - v for k, v in a.fi_original.dicc_montana_valle.items()}
    m.setObjective(energia(x, orig, invertida), gb.GRB.MINIMIZE)
    m.optimize()

    assert m.SolCount > 0
    assert int(round(m.ObjVal)) == 0
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
