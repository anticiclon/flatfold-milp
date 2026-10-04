# -*- coding: utf-8 -*-
import pytest
from pathlib import Path

from flatfold.model import (construye_modelo, literal, resuelve,
                            orden_de_capas, añade_restricciones, analiza, 
                            separa_tripletas)

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"

RUTAS = {
    "single": TOY / "single_fold.fold",
    "square": TOY / "square.fold",
    "map13":  TOY / "map_1x3.fold",
    "bomb":   TOY / "waterbomb4.fold",
    "bomb8":  TOY / "waterbomb8.fold",
}


def modelo_completo(datos, **banderas):
    """Modelo con las restricciones activadas según las banderas."""
    pares, trip, cr, tt, tto, oo = datos
    m, x = construye_modelo(pares)
    # añade_restricciones(m, x, trip, cr, tt, tto, oo, **banderas)
    indice_pares = {par: i for i, par in enumerate(pares)}
    añade_restricciones(m, x, indice_pares, trip, cr, tt, tto, oo, **banderas)
    return m, x


def cuenta(datos, limite=5000, **banderas):
    m, _ = modelo_completo(datos, **banderas)
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.Params.OutputFlag = 0
    m.optimize()
    n = m.SolCount
    assert n < limite, "el pool se ha llenado: sube el límite"
    return n


# === LA VALIDACIÓN EXTERNA =================================================

@pytest.mark.parametrize("n,esperados", [
    (2, 2), (3, 6), (4, 16), (5, 50), (6, 144), (7, 462),
])
def test_stamp_folding(n, esperados):
    """
    Número de plegados planos de una tira de n celdas. Secuencia clásica
    del stamp folding, calculada por otros autores y por otro método: es
    la validación externa de toda la tubería.
    """
    assert cuenta(analiza(TOY / f"map_1x{n}.fold")) == esperados


# === recuentos conocidos ===================================================

def test_pliegue_unico_dos_estados():
    """Monte y valle."""
    assert cuenta(analiza(RUTAS["single"])) == 2


def test_cuadrado_un_estado():
    assert cuenta(analiza(RUTAS["square"])) == 1


def test_waterbomb4():
    assert cuenta(analiza(RUTAS["bomb"])) == 8


def test_cada_familia_restringe():
    """La caída 64 -> 24 -> 8 en la bomba de agua de grado 4."""
    d = analiza(RUTAS["bomb"])
    sin = dict(con_transitividad=False, con_cruces=False, con_colineales=False)
    solo_trans = dict(con_cruces=False, con_colineales=False)
    assert cuenta(d, **sin) == 64            # 2^6
    assert cuenta(d, **solo_trans) == 24     # 4! ordenes totales
    assert cuenta(d) == 8


def test_transitividad_da_los_ordenes_totales():
    """Con n caras apiladas y solo transitividad quedan n! estados."""
    for n, factorial in [(3, 6), (4, 24), (5, 120)]:
        d = analiza(TOY / f"map_1x{n}.fold")
        assert cuenta(d, con_cruces=False, con_colineales=False) == factorial


def test_numero_de_estados_es_par():
    """El conjunto de estados es cerrado bajo la inversión global."""
    for nombre in ["single", "map13", "bomb"]:
        assert cuenta(analiza(RUTAS[nombre])) % 2 == 0


# === tamaños del modelo ====================================================

@pytest.mark.parametrize("nombre,n_vars", [
    ("single", 1), ("square", 0), ("map13", 3), ("bomb", 6), ("bomb8", 28),
])
def test_numero_de_variables(nombre, n_vars):
    pares = analiza(RUTAS[nombre])[0]
    m, x = construye_modelo(pares)
    assert len(x) == n_vars
    assert m.NumVars == n_vars


@pytest.mark.parametrize("n", [3, 4, 5, 6])
def test_tamanos_crecen_como_toca(n):
    """n caras apiladas: C(n,2) variables y C(n,3) tripletas."""
    from math import comb
    pares, trip, *_ = analiza(TOY / f"map_1x{n}.fold")
    assert len(pares) == comb(n, 2)
    assert len(trip) == comb(n, 3)


# === validaciones enganchadas a la cadena ==================================

def test_pendiente_ilegal_para_la_cadena():
    with pytest.raises(ValueError):
        analiza(TOY / "bad_slope.fold")


def test_kawasaki_ilegal_para_la_cadena():
    with pytest.raises(ValueError):
        analiza(TOY / "bad_kawasaki.fold")


# === literal ===============================================================

def test_literal_devuelve_la_variable():
    m, x = construye_modelo([(0, 2)])
    assert literal(x, 0, 2) is x[0, 2]


def test_literal_antisimetrico():
    m, x = construye_modelo([(0, 2)])
    suma = literal(x, 0, 2) + literal(x, 2, 0)
    assert suma.getConstant() == 1
    assert suma.size() == 2
    assert suma.getCoeff(0) + suma.getCoeff(1) == 0


def test_literal_sobre_la_misma_cara_falla():
    m, x = construye_modelo([(0, 2)])
    with pytest.raises(ValueError):
        literal(x, 1, 1)


def test_literal_de_par_inexistente_falla():
    m, x = construye_modelo([(0, 2)])
    with pytest.raises(KeyError):
        literal(x, 0, 5)


# === restricciones sueltas =================================================

def test_transitividad_excluye_los_ciclos():
    """De las 8 asignaciones de tres comparaciones, sobreviven 6."""
    from flatfold.model import añade_transitividad
    m, x = construye_modelo([(0, 1), (0, 2), (1, 2)])
    añade_transitividad(m, x, [(0, 1, 2)])
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = 20
    m.Params.OutputFlag = 0
    m.optimize()
    assert m.SolCount == 6


def test_igualdad_reduce_a_la_mitad():
    from flatfold.model import añade_igualdades
    m, x = construye_modelo([(0, 1), (2, 3)])
    añade_igualdades(m, x, [((0, 1), (2, 3))])
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = 20
    m.Params.OutputFlag = 0
    m.optimize()
    assert m.SolCount == 2          # de 4 a 2


def test_taco_taco_deja_paridad_par():
    """De las 16 asignaciones de cuatro literales, las 8 de paridad par."""
    from flatfold.model import añade_taco_taco
    pares = [(0, 1), (0, 2), (1, 3), (2, 3)]
    m, x = construye_modelo(pares)
    añade_taco_taco(m, x, [[(0, 1), (0, 2), (1, 3), (2, 3)]])
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = 50
    m.Params.OutputFlag = 0
    m.optimize()
    assert m.SolCount == 8


# === orden de capas ========================================================

def test_orden_de_capas_dos_caras():
    assert orden_de_capas([(0, 2)], {(0, 2): 1}) == [[0, 2]]
    assert orden_de_capas([(0, 2)], {(0, 2): 0}) == [[2, 0]]


def test_orden_de_capas_tres_caras():
    pares = [(0, 1), (0, 2), (1, 2)]
    assert orden_de_capas(pares, {(0,1): 1, (0,2): 1, (1,2): 1}) == [[0, 1, 2]]
    assert orden_de_capas(pares, {(0,1): 0, (0,2): 0, (1,2): 0}) == [[2, 1, 0]]


def test_orden_de_capas_dos_grupos():
    assert orden_de_capas([(0, 1), (2, 3)],
                          {(0, 1): 1, (2, 3): 0}) == [[0, 1], [3, 2]]


def test_modelo_sin_variables_es_factible():
    m, x = construye_modelo([])
    estado, valores = resuelve(m, x)
    assert valores == {}
    assert orden_de_capas([], valores) == []


def test_la_cadena_produce_solucion():
    for nombre in RUTAS:
        datos = analiza(RUTAS[nombre])
        m, x = modelo_completo(datos)
        estado, valores = resuelve(m, x)
        assert valores is not None
        
def test_separa_detecta_el_ciclo():
    """f encima de g, g encima de h, h encima de f: imposible."""
    valores = {(0, 1): 1, (1, 2): 1, (0, 2): 0}
    # X(0,1)=1, X(1,2)=1, X(2,0)=1-X(0,2)=1  -> suma 3
    assert separa_tripletas(valores, [(0, 1, 2)]) == [(0, 1, 2)]


def test_separa_acepta_orden_total():
    valores = {(0, 1): 1, (1, 2): 1, (0, 2): 1}
    # X(0,1)=1, X(1,2)=1, X(2,0)=0  -> suma 2
    assert separa_tripletas(valores, [(0, 1, 2)]) == []


def test_separa_tolera_flotantes():
    valores = {(0, 1): 0.9999999, (1, 2): 1.0, (0, 2): 1e-9}
    assert separa_tripletas(valores, [(0, 1, 2)]) == [(0, 1, 2)]
    
def test_separa_detecta_el_ciclo_inverso():
    valores = {(0, 1): 0, (1, 2): 0, (0, 2): 1}
    assert separa_tripletas(valores, [(0, 1, 2)]) == [(0, 1, 2)]


def test_separa_acepta_los_ordenes_totales():
    """De las 8 asignaciones, 6 son órdenes totales y no se separan."""
    from itertools import product
    n_violadas = 0
    for a, b, c in product((0, 1), repeat=3):
        valores = {(0, 1): a, (1, 2): b, (0, 2): c}
        if separa_tripletas(valores, [(0, 1, 2)]):
            n_violadas += 1
    assert n_violadas == 2          # los dos ciclos

def test_separa_lista_vacia():
    assert separa_tripletas({}, []) == []
        
        
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

        