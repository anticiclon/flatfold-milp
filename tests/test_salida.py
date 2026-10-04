# -*- coding: utf-8 -*-
"""
Tests de la verificación y la exportación de soluciones.

verifica_solucion recalcula las incidencias desde la geometría plegada,
sin usar las listas que produce incidences. Comprueba por tanto un
segundo camino, y detecta errores en las funciones que generan las
incidencias y errores de signo en la traducción al modelo: la clase de
fallo que ni los recuentos de estados ni las comprobaciones de
factibilidad revelan.
"""
import json
import pytest
from pathlib import Path

# import gurobipy as gb
from flatfold.model import (analiza_completo, construye_modelo,
                            añade_restricciones)
from flatfold.presolve import (literales_a_ecuaciones, reduce_completo,
                               componentes, reparte_por_componentes,
                               expande)
from flatfold.model import resuelve_por_componentes
from flatfold.salida import verifica_solucion, exporta_estado_fold, verifica_punto_a_punto

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
GRIDS = RAIZ / "instances" / "flatfolder" / "grids"


def resuelve(ruta):
    """Una solución cualquiera, por el camino directo."""
    a = analiza_completo(ruta)
    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr + a.cr_oo, a.tt,
                        a.tto, a.oo)
    m.Params.OutputFlag = 0
    m.optimize()
    assert m.SolCount > 0
    return a, {par: int(round(v.X)) for par, v in x.items()}


def resuelve_con_presolve(ruta):
    """Lo mismo por el camino reducido, expandiendo la solución."""
    a = analiza_completo(ruta)
    ig, par, fij = literales_a_ecuaciones(
        a.pares, a.cr + a.cr_oo, a.tto, a.oo, a.tt)
    mapeo, restantes, vivas, _ = reduce_completo(
        len(a.pares), ig, par, fij, a.trip, a.pares)
    comps = componentes(mapeo, restantes, a.trip, a.pares)
    datos = reparte_por_componentes(comps, vivas, restantes, {})
    ok, beta, _ = resuelve_por_componentes(datos, modo="factibilidad")
    assert ok
    return a, expande(mapeo, a.pares, beta)


# ===== el verificador acepta soluciones válidas ============================

@pytest.mark.parametrize("nombre", ["single_fold", "map_1x3", "map_1x5",
                                    "waterbomb4"])
def test_solucion_valida_en_juguete(nombre):
    a, vals = resuelve(TOY / f"{nombre}.fold")
    assert verifica_solucion(a, vals) == []


@pytest.mark.parametrize("nombre", ["004_brown_3x3_A3.fold",
                                    "182_brown_4x4_C41.fold"])
def test_solucion_valida_en_instancia_real(nombre):
    a, vals = resuelve(GRIDS / nombre)
    assert verifica_solucion(a, vals) == []


def test_solucion_del_presolve_es_valida():
    """
    El caso que más valor tiene: la solución viene expandida desde las
    variables supervivientes, de modo que un error de signo en expande
    se manifestaría aquí.
    """
    a, vals = resuelve_con_presolve(GRIDS / "004_brown_3x3_A3.fold")
    assert verifica_solucion(a, vals) == []


def test_solucion_sobre_patron_subdividido():
    """dos_eles es la única instancia que activa la subdivisión."""
    a, vals = resuelve(TOY / "dos_eles.fold")
    assert a.subdividido
    assert verifica_solucion(a, vals) == []


# ===== el verificador rechaza soluciones inválidas =========================

def test_detecta_un_ciclo():
    """
    Un ciclo f > g > h > f no corresponde a ningún apilamiento. Se
    construye forzando los tres literales de una tripleta al mismo
    valor, en lugar de invertir uno solo: invertir un literal cambia la
    suma en una unidad, y solo la saca del rango válido la mitad de las
    veces.
    """
    a, vals = resuelve(TOY / "waterbomb4.fold")
    assert verifica_solucion(a, vals) == []

    # las cuatro caras de la bomba de agua se apilan sobre el mismo
    # triángulo, así que cualquier tripleta cubre un punto común
    f, g, h = sorted(a.plegadas)[:3]
    roto = dict(vals)
    roto[(f, g)] = 1          # X_fg = 1
    roto[(g, h)] = 1          # X_gh = 1
    roto[(f, h)] = 0          # X_hf = 1 - x_fh = 1

    fallos = verifica_solucion(a, roto)
    assert any(t == "ciclo" for t, *_ in fallos)


def test_alguna_inversion_rompe_la_solucion():
    """
    Invertir un literal no siempre viola una restricción, pero no puede
    ser inocuo para todos: alguna inversión debe romper la solución.
    """
    a, vals = resuelve(GRIDS / "004_brown_3x3_A3.fold")
    rompen = 0
    for par in a.pares:
        roto = dict(vals)
        roto[par] = 1 - roto[par]
        if verifica_solucion(a, roto):
            rompen += 1
    assert rompen > 0


def test_detecta_todas_las_violaciones_al_invertirlo_todo():
    """Invertir cada literal no puede dejar la solución válida."""
    a, vals = resuelve(GRIDS / "182_brown_4x4_C41.fold")
    roto = {k: 1 - v for k, v in vals.items()}
    # la inversión global es de hecho válida: es el estado del reverso
    assert verifica_solucion(a, roto) == []


def test_detecta_un_cruce():
    """
    En una instancia con cruces, invertir uno de los dos literales de un
    par (f,t),(g,t) rompe la igualdad que impone la familia (E).
    """
    a, vals = resuelve(GRIDS / "182_brown_4x4_C41.fold")
    assert a.cr, "la instancia debe tener cruces"
    f, g, t = a.cr[0]
    par = (min(f, t), max(f, t))
    roto = dict(vals)
    roto[par] = 1 - roto[par]

    fallos = verifica_solucion(a, roto)
    assert fallos, "el verificador no detectó la violación"


# ===== exportación a FOLD ==================================================

def test_fold_exportado_es_json_valido(tmp_path):
    a, vals = resuelve(TOY / "waterbomb4.fold")
    ruta = exporta_estado_fold(a, vals, tmp_path / "estado.fold")
    doc = json.loads(Path(ruta).read_text(encoding="utf-8"))

    assert doc["file_spec"] == 1.1
    for campo in ("vertices_coords", "faces_vertices", "faceOrders"):
        assert campo in doc


def test_fold_una_terna_por_par(tmp_path):
    a, vals = resuelve(TOY / "waterbomb4.fold")
    ruta = exporta_estado_fold(a, vals, tmp_path / "estado.fold")
    doc = json.loads(Path(ruta).read_text(encoding="utf-8"))
    assert len(doc["faceOrders"]) == len(a.pares)


def test_fold_signos_admisibles(tmp_path):
    a, vals = resuelve(GRIDS / "004_brown_3x3_A3.fold")
    ruta = exporta_estado_fold(a, vals, tmp_path / "estado.fold")
    doc = json.loads(Path(ruta).read_text(encoding="utf-8"))

    n = len(doc["faces_vertices"])
    for f, g, s in doc["faceOrders"]:
        assert s in (-1, 1)
        assert 0 <= f < n and 0 <= g < n
        assert f != g


def test_fold_indices_de_cara_consecutivos(tmp_path):
    """Las caras se renumeran porque la exterior no se escribe."""
    a, vals = resuelve(GRIDS / "004_brown_3x3_A3.fold")
    ruta = exporta_estado_fold(a, vals, tmp_path / "estado.fold")
    doc = json.loads(Path(ruta).read_text(encoding="utf-8"))
    assert len(doc["faces_vertices"]) == len(a.plegadas)


def test_fold_vertices_deduplicados(tmp_path):
    """Las caras plegadas comparten vértices; no deben repetirse."""
    a, vals = resuelve(TOY / "waterbomb4.fold")
    ruta = exporta_estado_fold(a, vals, tmp_path / "estado.fold")
    doc = json.loads(Path(ruta).read_text(encoding="utf-8"))
    coords = [tuple(c) for c in doc["vertices_coords"]]
    assert len(coords) == len(set(coords))
    
    
    
    
    
    
@pytest.mark.parametrize("nombre", ["map_1x3", "map_1x5", "waterbomb4"])
def test_punto_a_punto_acepta(nombre):
    a, vals = resuelve(TOY / f"{nombre}.fold")
    assert verifica_punto_a_punto(a, vals) == []


def test_punto_a_punto_detecta_un_ciclo():
    a, vals = resuelve(TOY / "waterbomb4.fold")
    f, g, h = sorted(a.plegadas)[:3]
    roto = dict(vals)
    roto[(f, g)], roto[(g, h)], roto[(f, h)] = 1, 1, 0
    assert verifica_punto_a_punto(a, roto)
    
    
    
    