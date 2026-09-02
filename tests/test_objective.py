# -*- coding: utf-8 -*-
import pytest
from pathlib import Path
import gurobipy as gb

from flatfold.io_fold import read_fold, Kind
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.model import construye_modelo, añade_restricciones, literal
from flatfold.objective import etiquetas_mv, fija_etiquetas, energia
from flatfold.exact import sigma

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"


def prepara(nombre):
    """Modelo completo con restricciones, más las etiquetas M/V."""
    fi = read_fold(TOY / f"{nombre}.fold")
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, exterior = caras_definitivas(G, caras)
    idx = caras.index(exterior)
    dual = grafo_dual(G, caras, idx)
    dic = aristas_caras(G, caras)
    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}

    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    cr = cruces(G, caras, plegadas, iso, dic, idx)
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)

    m, x = construye_modelo(pares)
    añade_restricciones(m, x, trip, cr, tt, tto, oo)
    m.Params.OutputFlag = 0

    etiq = etiquetas_mv(G, fi, iso, dic, idx)
    return fi, G, iso, dic, idx, m, x, etiq


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
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"
