# -*- coding: utf-8 -*-
"""
Tests de la conexión del presolve con el modelo. El criterio de fondo es
que el sistema reducido tiene EXACTAMENTE las mismas soluciones que el
completo: mismo recuento de estados y mismo óptimo.
"""
import pytest
from pathlib import Path
import gurobipy as gb

from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.objective import etiquetas_mv, energia
from flatfold.model import (construye_modelo, añade_restricciones,
                            construye_reducido, resuelve_por_componentes)
from flatfold.presolve import (literales_a_ecuaciones, reduce_completo,
                               componentes, energia_reducida, expande,
                               reparte_por_componentes)

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"
GRIDS = (Path(__file__).resolve().parent.parent
         / "instances" / "flatfolder" / "grids")


def monta(ruta):
    fi = read_fold(ruta)
    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, ext = caras_definitivas(G, caras)
    idx = caras.index(ext)
    dual = grafo_dual(G, caras, idx)
    dic = aristas_caras(G, caras)
    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    cr = cruces(G, caras, plegadas, iso, dic, idx)
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic, idx)
    etiq = etiquetas_mv(G, fi, iso, dic, idx)
    return fi, pares, trip, cr, tt, tto, oo, etiq


def reduce_todo(pares, trip, cr, tt, tto, oo):
    ig, par, fij = literales_a_ecuaciones(pares, cr, tto, oo, tt)
    mapeo, restantes, vivas, _ = reduce_completo(
        len(pares), ig, par, fij, trip, pares)
    comps = componentes(mapeo, restantes, trip, pares)
    return mapeo, restantes, vivas, comps


def cuenta_completo(pares, trip, cr, tt, tto, oo, limite=100000):
    m, x = construye_modelo(pares)
    # añade_restricciones(m, x, trip, cr, tt, tto, oo)
    indice_pares = {par: i for i, par in enumerate(pares)}
    añade_restricciones(m, x, indice_pares, trip, cr, tt, tto, oo)
    m.Params.OutputFlag = 0
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.optimize()
    assert m.SolCount < limite, "el pool se llenó"
    return m.SolCount


# ===== reparto =============================================================

@pytest.mark.parametrize("nombre", ["waterbomb4", "waterbomb8", "map_1x5"])
def test_reparto_no_cruza_componentes(nombre):
    """Toda restricción cae entera en una componente; si no, los assert
    de reparte_por_componentes saltan."""
    _, pares, trip, cr, tt, tto, oo, _ = monta(TOY / f"{nombre}.fold")
    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    assert sum(len(d["tripletas"]) for d in datos) == len(vivas)
    assert sum(len(d["paridades"]) for d in datos) == len(restantes)


def test_reparto_instancia_real():
    _, pares, trip, cr, tt, tto, oo, _ = monta(GRIDS / "003_brown_3x3_A2.fold")
    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    raices = {r for d in datos for r in d["raices"]}
    supervivientes = {e[1] for e in mapeo.values() if e[0] == "var"}
    assert raices == supervivientes


# ===== el recuento no cambia ===============================================

@pytest.mark.parametrize("nombre", ["map_1x3", "map_1x5", "waterbomb4"])
def test_recuento_coincide_juguete(nombre):
    """El presolve no puede alterar el número de estados."""
    _, pares, trip, cr, tt, tto, oo, _ = monta(TOY / f"{nombre}.fold")
    esperado = cuenta_completo(pares, trip, cr, tt, tto, oo)

    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    obtenido, _, _ = resuelve_por_componentes(datos, modo="recuento")
    assert obtenido == esperado


def test_recuento_coincide_instancia_real():
    """Aquí el presolve sí reduce (de 12 variables a 8)."""
    _, pares, trip, cr, tt, tto, oo, _ = monta(GRIDS / "003_brown_3x3_A2.fold")
    esperado = cuenta_completo(pares, trip, cr, tt, tto, oo)

    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    obtenido, _, _ = resuelve_por_componentes(datos, modo="recuento")
    assert obtenido == esperado


@pytest.mark.parametrize("n,esperados", [(3, 6), (4, 16), (5, 50)])
def test_stamp_folding_con_presolve(n, esperados):
    """La secuencia clásica debe salir igual con el sistema reducido."""
    _, pares, trip, cr, tt, tto, oo, _ = monta(TOY / f"map_1x{n}.fold")
    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    obtenido, _, _ = resuelve_por_componentes(datos, modo="recuento")
    assert obtenido == esperados


# ===== factibilidad ========================================================

@pytest.mark.parametrize("nombre", ["map_1x3", "waterbomb4", "waterbomb8"])
def test_factibilidad_coincide(nombre):
    _, pares, trip, cr, tt, tto, oo, _ = monta(TOY / f"{nombre}.fold")
    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    ok, valores, _ = resuelve_por_componentes(datos, modo="factibilidad")
    assert ok
    assert set(valores) == {r for d in datos for r in d["raices"]}


def test_solucion_reducida_es_valida():
    """Expandir la solución del sistema reducido debe dar una solución
    válida del completo."""
    _, pares, trip, cr, tt, tto, oo, _ = monta(GRIDS / "003_brown_3x3_A2.fold")
    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    ok, beta, _ = resuelve_por_componentes(datos, modo="factibilidad")
    assert ok

    vals = expande(mapeo, pares, beta)

    idx = {p: i for i, p in enumerate(pares)}

    def X(u, v):
        return vals[(u, v)] if u < v else 1 - vals[(v, u)]

    for f, g, h in trip:
        assert 1 <= X(f, g) + X(g, h) + X(h, f) <= 2
    for f, g, t in cr:
        assert X(f, t) == X(g, t)
    for cuatro in tt:
        assert sum(X(u, v) for u, v in cuatro) % 2 == 0


# ===== energía =============================================================

def test_energia_reducida_coincide_con_la_completa():
    """
    Resuelve el completo, traduce su solución a las raíces, y comprueba
    que la energía reducida da el mismo valor.
    """
    fi, pares, trip, cr, tt, tto, oo, etiq = monta(TOY / "waterbomb4.fold")

    m, x = construye_modelo(pares)
    # añade_restricciones(m, x, trip, cr, tt, tto, oo)
    indice_pares = {par: i for i, par in enumerate(pares)}
    añade_restricciones(m, x, indice_pares, trip, cr, tt, tto, oo)
    m.Params.OutputFlag = 0
    invertida = {a: 1 - v for a, v in fi.dicc_montana_valle.items()}
    m.setObjective(energia(x, etiq, invertida), gb.GRB.MINIMIZE)
    m.optimize()
    assert m.SolCount > 0
    optimo_completo = int(round(m.ObjVal))

    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    terminos, const = energia_reducida(mapeo, pares, etiq, invertida)
    datos = reparte_por_componentes(comps, vivas, restantes, terminos)
    valor, _, _ = resuelve_por_componentes(datos, const, modo="optimizacion")
    assert int(round(valor)) == optimo_completo


def test_optimo_cero_si_la_asignacion_es_realizable():
    fi, pares, trip, cr, tt, tto, oo, etiq = monta(TOY / "waterbomb4.fold")
    mapeo, restantes, vivas, comps = reduce_todo(pares, trip, cr, tt, tto, oo)
    terminos, const = energia_reducida(mapeo, pares, etiq,
                                       fi.dicc_montana_valle)
    datos = reparte_por_componentes(comps, vivas, restantes, terminos)
    valor, _, _ = resuelve_por_componentes(datos, const, modo="optimizacion")
    assert int(round(valor)) == 0
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

    
    
    
    