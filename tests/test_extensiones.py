# -*- coding: utf-8 -*-
"""
Tests de los objetivos y restricciones adicionales: profundidad de una
cara, exposición en superficie, cota de profundidad de una región,
caras intermedias entre dos caras, recuento de montes, distancia de
Hamming y enumeración de asignaciones.

Todas estas cantidades son lineales en las variables de capa y no alteran
la estructura del programa, de modo que lo que hay que comprobar es que
calculan lo que dicen y que las restricciones que generan no cambian el
conjunto de estados válidos más de lo que deben.
"""
import pytest
from pathlib import Path

import gurobipy as gb
from flatfold.model import (analiza_completo, construye_modelo,
                            añade_restricciones, profundidad,
                            añade_en_superficie, añade_profundidad_maxima,
                            literal, añade_caras_intermedias)
from flatfold.objective import (etiquetas_mv, etiquetas_por_pliegue_original,
                                energia, cuenta_montes, distancia_hamming,
                                enumera_asignaciones)

RAIZ = Path(__file__).resolve().parent.parent
TOY = RAIZ / "instances" / "toy"
GRIDS = RAIZ / "instances" / "flatfolder" / "grids"


def modelo(ruta):
    """Analisis y modelo con las cinco familias, sin objetivo."""
    a = analiza_completo(ruta)
    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr + a.cr_oo, a.tt,
                        a.tto, a.oo)
    m.Params.OutputFlag = 0
    return a, m, x


def cuenta_estados(a, m, x, limite=5000):
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.optimize()
    assert m.SolCount < limite, "el pool se llenó"
    return m.SolCount


# ===== profundidad =========================================================

def test_profundidad_suma_las_caras_encima():
    """Sobre una solución concreta, el valor de la expresión debe ser el
    número de caras que están por encima."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    m.optimize()
    assert m.SolCount > 0

    vals = {par: int(round(v.X)) for par, v in x.items()}

    def X(u, v):
        return vals[(u, v)] if u < v else 1 - vals[(v, u)]

    for f in a.plegadas:
        encima = sum(1 for g in a.plegadas
                     if g != f and (min(f, g), max(f, g)) in vals
                     and X(g, f) == 1)
        assert round(profundidad(x, a.pares, f).getValue()) == encima


def test_profundidades_son_una_permutacion():
    """Las cuatro caras de la bomba de agua se apilan, así que sus
    profundidades son 0, 1, 2 y 3."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    m.optimize()
    ds = sorted(round(profundidad(x, a.pares, f).getValue())
                for f in a.plegadas)
    assert ds == [0, 1, 2, 3]


def test_profundidad_de_cara_aislada_es_cero():
    """Una cara que no solapa con ninguna está siempre en la superficie."""
    a, m, x = modelo(TOY / "map_1x3.fold")
    aisladas = [f for f in a.plegadas
                if not any(f in par for par in a.pares)]
    for f in aisladas:
        assert profundidad(x, a.pares, f).size() == 0


# ===== en superficie =======================================================

def test_en_superficie_deja_la_cara_arriba():
    a, m, x = modelo(TOY / "waterbomb4.fold")
    f = sorted(a.plegadas)[0]
    n = añade_en_superficie(m, x, a.pares, f)
    assert n == sum(1 for par in a.pares if f in par)

    m.optimize()
    assert m.SolCount > 0
    assert round(profundidad(x, a.pares, f).getValue()) == 0


@pytest.mark.parametrize("i", [0, 1, 2, 3])
def test_cualquier_cara_puede_estar_arriba(i):
    """En la bomba de agua las cuatro caras se apilan sobre el mismo
    triángulo, y cualquiera de ellas admite quedar en la superficie."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    f = sorted(a.plegadas)[i]
    añade_en_superficie(m, x, a.pares, f)
    m.optimize()
    assert m.SolCount > 0


def test_dos_caras_solapadas_no_pueden_estar_ambas_arriba():
    """Si f y g se solapan, exigir las dos en superficie es imposible."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    f, g = a.pares[0]
    añade_en_superficie(m, x, a.pares, f)
    añade_en_superficie(m, x, a.pares, g)
    m.optimize()
    assert m.SolCount == 0


def test_en_superficie_restringe_los_estados():
    """Fijar una cara arriba deja un subconjunto propio de los estados."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    total = cuenta_estados(a, m, x)

    a2, m2, x2 = modelo(TOY / "waterbomb4.fold")
    añade_en_superficie(m2, x2, a2.pares, sorted(a2.plegadas)[0])
    con = cuenta_estados(a2, m2, x2)

    assert 0 < con < total


# ===== profundidad máxima ==================================================

def test_cota_acota_todas_las_caras():
    a, m, x = modelo(TOY / "waterbomb4.fold")
    region = sorted(a.plegadas)[:2]
    t = añade_profundidad_maxima(m, x, a.pares, region)
    m.setObjective(t, gb.GRB.MINIMIZE)
    m.optimize()
    assert m.SolCount > 0

    valor = round(t.X)
    for f in region:
        assert round(profundidad(x, a.pares, f).getValue()) <= valor


def test_cota_minima_de_dos_caras_apiladas():
    """Dos caras que se solapan no pueden estar ambas en la superficie,
    así que la profundidad máxima del par es al menos 1."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    f, g = a.pares[0]
    t = añade_profundidad_maxima(m, x, a.pares, [f, g])
    m.setObjective(t, gb.GRB.MINIMIZE)
    m.optimize()
    assert round(t.X) >= 1


def test_cota_con_conjunto_vacio():
    a, m, x = modelo(TOY / "waterbomb4.fold")
    with pytest.raises(ValueError):
        añade_profundidad_maxima(m, x, a.pares, [])


# ===== recuento de montes ==================================================

def test_montes_mas_valles_es_el_total():
    a, m, x = modelo(TOY / "waterbomb4.fold")
    m.optimize()
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)

    n = round(cuenta_montes(x, orig).getValue())
    assert 0 <= n <= len(orig)


def test_montes_coincide_con_las_etiquetas_leidas():
    """El recuento debe igualar al que se obtiene evaluando cada m_e."""
    a, m, x = modelo(GRIDS / "004_brown_3x3_A3.fold")
    m.optimize()
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)

    a_mano = 0
    for (u, v), s in orig.values():
        xv = int(round(x[u, v].X))
        a_mano += xv if s == 0 else 1 - xv
    assert round(cuenta_montes(x, orig).getValue()) == a_mano


def test_prescribir_el_numero_de_montes():
    """Fijar el recuento al que tiene la asignación del fichero debe
    seguir siendo factible, puesto que esa asignación es realizable."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    objetivo = sum(a.fi_original.dicc_montana_valle.values())

    m.addConstr(cuenta_montes(x, orig) == objetivo)
    m.optimize()
    assert m.SolCount > 0


# ===== distancia de Hamming ================================================

def test_distancia_a_si_misma_es_cero():
    a, m, x = modelo(TOY / "waterbomb4.fold")
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    m0 = a.fi_original.dicc_montana_valle

    m.setObjective(distancia_hamming(x, orig, m0), gb.GRB.MINIMIZE)
    m.optimize()
    assert m.SolCount > 0
    assert round(m.ObjVal) == 0, "la asignación del fichero es realizable"


def test_distancia_coincide_con_la_energia_unitaria():
    """La distancia es la energía de pesos unitarios; comprobamos que
    ambas dan el mismo valor sobre la misma solución."""
    a, m, x = modelo(GRIDS / "004_brown_3x3_A3.fold")
    m.optimize()
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    m0 = {k: 1 - v for k, v in a.fi_original.dicc_montana_valle.items()}

    assert (round(distancia_hamming(x, orig, m0).getValue())
            == round(energia(x, orig, m0).getValue()))


def test_distancia_contada_a_mano():
    a, m, x = modelo(GRIDS / "004_brown_3x3_A3.fold")
    m.optimize()
    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    m0 = a.fi_original.dicc_montana_valle

    difieren = 0
    for arista, deseado in m0.items():
        if arista not in orig:
            continue
        (u, v), s = orig[arista]
        xv = int(round(x[u, v].X))
        me = xv if s == 0 else 1 - xv
        if me != deseado:
            difieren += 1
    assert round(distancia_hamming(x, orig, m0).getValue()) == difieren


# ===== enumeración de asignaciones =========================================

def test_asignaciones_son_distintas_entre_si():
    a = analiza_completo(TOY / "map_1x3.fold")
    asigs, agotado = enumera_asignaciones(a)
    assert agotado, "deberían agotarse en una instancia tan pequeña"
    claves = {tuple(sorted(d.items())) for d in asigs}
    assert len(claves) == len(asigs)


def test_hay_menos_asignaciones_que_estados():
    """Un mismo etiquetado puede realizarse por varios apilamientos, de
    modo que las asignaciones nunca superan a los estados."""
    a, m, x = modelo(TOY / "map_1x3.fold")
    n_estados = cuenta_estados(a, m, x)
    asigs, _ = enumera_asignaciones(analiza_completo(TOY / "map_1x3.fold"))
    assert 0 < len(asigs) <= n_estados


def test_cada_asignacion_es_realizable():
    """Fijar cada una de las asignaciones enumeradas debe dar un modelo
    factible."""
    from flatfold.objective import fija_etiquetas
    ruta = TOY / "map_1x3.fold"
    asigs, _ = enumera_asignaciones(analiza_completo(ruta))

    for m0 in asigs:
        a, m, x = modelo(ruta)
        etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
        orig = etiquetas_por_pliegue_original(etiq, a.representantes)
        fija_etiquetas(m, x, orig, m0)
        m.optimize()
        assert m.SolCount > 0, m0


def test_la_asignacion_del_fichero_esta_entre_las_enumeradas():
    ruta = TOY / "waterbomb4.fold"
    a = analiza_completo(ruta)
    asigs, agotado = enumera_asignaciones(a)
    assert agotado
    assert a.fi_original.dicc_montana_valle in asigs


def test_una_asignacion_irrealizable_no_aparece():
    """waterbomb8 es geométricamente válido pero su asignación viola
    Maekawa, de modo que no puede estar entre las enumeradas."""
    ruta = TOY / "waterbomb8.fold"
    a = analiza_completo(ruta)
    asigs, _ = enumera_asignaciones(a, limite=200)
    assert a.fi_original.dicc_montana_valle not in asigs


def test_el_limite_se_respeta():
    a = analiza_completo(TOY / "waterbomb8.fold")
    asigs, agotado = enumera_asignaciones(a, limite=3)
    assert len(asigs) == 3
    assert not agotado
    
# ===== caras intermedias ===================================================

def modelo_sintetico(pares):
    """Modelo con solo las variables de capa, sin ninguna familia.

    Las caras intermedias no dependen de la geometría, sino de qué pares
    existen y de qué valen sus literales, de modo que estos tests fijan
    los literales a mano sobre configuraciones inventadas. Así se recorre
    la tabla de verdad completa, que ninguna instancia real garantiza
    visitar.
    """
    m, x = construye_modelo(pares)
    m.Params.OutputFlag = 0
    return m, x


def pares_completos(n):
    """Todas las caras de 0..n-1 solapan entre sí."""
    return [(u, v) for u in range(n) for v in range(u + 1, n)]


@pytest.mark.parametrize("v_ft,v_gt,esperado", [(0, 0, 0), (1, 0, 1),
                                                (0, 1, 1), (1, 1, 0)])
def test_caras_intermedias_es_el_xor(v_ft, v_gt, esperado):
    """Las cuatro desigualdades fijan y_t al XOR en los cuatro casos.

    Se resuelve dos veces, minimizando y maximizando: si solo se
    minimizara, una y_t meramente acotada por abajo pasaría el test
    igual, y no se comprobaría que las dos cotas superiores actúan.
    Exigiendo que ambos sentidos coincidan se confirma que y_t queda
    determinada y no solo acotada.
    """
    caras = list(range(3))
    pares = pares_completos(3)
    f, g, t = 0, 1, 2

    for sentido in (gb.GRB.MINIMIZE, gb.GRB.MAXIMIZE):
        m, x = modelo_sintetico(pares)
        N, y = añade_caras_intermedias(m, x, pares, caras, None, f, g)
        m.addConstr(literal(x, f, t) == v_ft)
        m.addConstr(literal(x, g, t) == v_gt)
        m.setObjective(N, sentido)
        m.optimize()
        assert m.Status == gb.GRB.OPTIMAL
        assert round(y[t].X) == esperado
        assert round(N.getValue()) == esperado


def test_caras_intermedias_sin_candidatas_vale_cero():
    """Si ninguna cara solapa con f y con g a la vez, N_fg es la
    constante 0: no hay nada capaz de intercalarse."""
    caras = list(range(4))
    pares = [(0, 1), (0, 2), (1, 3)]      # 2 solo con 0; 3 solo con 1
    m, x = modelo_sintetico(pares)
    N, y = añade_caras_intermedias(m, x, pares, caras, None, 0, 1)
    assert y == {}
    assert N.getValue() == 0


def test_caras_intermedias_exige_que_f_y_g_se_solapen():
    """Sin región común no hay grosor que medir."""
    caras = list(range(3))
    pares = [(0, 2), (1, 2)]              # 0 y 1 no solapan entre sí
    m, x = modelo_sintetico(pares)
    with pytest.raises(ValueError):
        añade_caras_intermedias(m, x, pares, caras, None, 0, 1)


def test_caras_intermedias_excluye_la_cara_exterior():
    """La cara exterior está en caras pero no es una capa del plegado."""
    caras = list(range(4))
    pares = pares_completos(4)
    m, x = modelo_sintetico(pares)
    _, y = añade_caras_intermedias(m, x, pares, caras, 3, 0, 1)
    assert set(y) == {2}


@pytest.mark.parametrize("v2,v3,esperado", [((0, 0), (0, 0), 0),
                                            ((1, 0), (0, 0), 1),
                                            ((1, 0), (1, 0), 2),
                                            ((1, 1), (0, 1), 1)])
def test_caras_intermedias_cuenta_varias(v2, v3, esperado):
    """Con dos candidatas, N_fg recorre 0, 1 y 2 según los literales."""
    caras = list(range(4))
    pares = pares_completos(4)
    f, g = 0, 1
    m, x = modelo_sintetico(pares)
    N, _ = añade_caras_intermedias(m, x, pares, caras, None, f, g)
    m.addConstr(literal(x, f, 2) == v2[0])
    m.addConstr(literal(x, g, 2) == v2[1])
    m.addConstr(literal(x, f, 3) == v3[0])
    m.addConstr(literal(x, g, 3) == v3[1])
    m.setObjective(N, gb.GRB.MINIMIZE)
    m.optimize()
    assert m.Status == gb.GRB.OPTIMAL
    assert round(N.getValue()) == esperado


def test_caras_intermedias_sobre_instancia_real():
    """Sobre un estado válido de la bomba de agua, N_fg debe coincidir
    con el recuento hecho leyendo las variables."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    f, g = a.pares[0]
    N, y = añade_caras_intermedias(m, x, a.pares, list(a.plegadas),
                                   a.idx_ext, f, g)
    m.optimize()
    assert m.SolCount > 0

    vals = {par: int(round(v.X)) for par, v in x.items()}

    def X(u, v):
        return vals[(u, v)] if u < v else 1 - vals[(v, u)]

    a_mano = sum(1 for t in y if X(f, t) != X(g, t))
    assert round(N.getValue()) == a_mano


def test_caras_intermedias_coincide_con_la_diferencia_de_profundidades():
    """En la bomba de agua las cuatro caras se apilan sobre el mismo
    triángulo, de modo que entre dos de ellas quedan exactamente las
    que su diferencia de profundidades indica."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    f, g = a.pares[0]
    N, _ = añade_caras_intermedias(m, x, a.pares, list(a.plegadas),
                                   a.idx_ext, f, g)
    m.optimize()
    df = round(profundidad(x, a.pares, f).getValue())
    dg = round(profundidad(x, a.pares, g).getValue())
    assert round(N.getValue()) == abs(df - dg) - 1


def test_caras_intermedias_acota_el_grosor():
    """Uso previsto: N_fg <= k como cota de grosor."""
    a, m, x = modelo(TOY / "waterbomb4.fold")
    f, g = a.pares[0]
    N, _ = añade_caras_intermedias(m, x, a.pares, list(a.plegadas),
                                   a.idx_ext, f, g)
    m.addConstr(N == 0)
    m.setObjective(N, gb.GRB.MAXIMIZE)
    m.optimize()
    assert m.SolCount > 0, "f y g deben poder quedar contiguas"
    assert round(N.getValue()) == 0
