# -*- coding: utf-8 -*-
import pytest
from pathlib import Path
from itertools import product

from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.incidences import (pares_comparables, tripletas, cruces,
                                 colineales)
from flatfold.objective import etiquetas_mv
from flatfold.presolve import (Infactible, UnionFindParidad, indice_literal,
                               literales_a_ecuaciones, reduce, expande,
                               componentes)

TOY = Path(__file__).resolve().parent.parent / "instances" / "toy"


# ===== UnionFindParidad ====================================================

def test_uf_transitividad():
    uf = UnionFindParidad(3)
    uf.union(0, 1, 0)
    uf.union(1, 2, 0)
    assert uf.find(0) == uf.find(2)


def test_uf_paridad_se_propaga():
    """x0 != x1 y x1 != x2  =>  x0 = x2."""
    uf = UnionFindParidad(3)
    uf.union(0, 1, 1)
    uf.union(1, 2, 1)
    r0, p0 = uf.find(0)
    r2, p2 = uf.find(2)
    assert r0 == r2 and p0 == p2


def test_uf_contradiccion():
    uf = UnionFindParidad(2)
    uf.union(0, 1, 0)
    with pytest.raises(Infactible):
        uf.union(0, 1, 1)


def test_uf_union_repetida_no_une():
    uf = UnionFindParidad(2)
    assert uf.union(0, 1, 0) is True
    assert uf.union(0, 1, 0) is False


def test_uf_fijacion_se_propaga():
    uf = UnionFindParidad(2)
    uf.fija(0, 1)
    uf.union(0, 1, 1)              # x1 != x0
    assert uf.determinado(0) == 1
    assert uf.determinado(1) == 0


def test_uf_no_determinado():
    uf = UnionFindParidad(2)
    uf.union(0, 1, 0)
    assert uf.determinado(0) is None


def test_uf_compresion_conserva_paridad():
    """Cadena larga: tras comprimir, las paridades deben seguir siendo
    correctas respecto de la raíz nueva."""
    n = 20
    uf = UnionFindParidad(n)
    for i in range(n - 1):
        uf.union(i, i + 1, 1)      # alterna
    r0, p0 = uf.find(0)
    for i in range(n):
        r, p = uf.find(i)
        assert r == r0
        assert p == (p0 ^ (i % 2))


# ===== indice_literal ======================================================

def test_indice_literal_signos():
    idx = {(0, 2): 0, (1, 3): 1}
    assert indice_literal(idx, 0, 2) == (0, 0)
    assert indice_literal(idx, 2, 0) == (0, 1)


def test_indice_literal_misma_cara():
    with pytest.raises(ValueError):
        indice_literal({(0, 2): 0}, 1, 1)


def test_indice_literal_par_inexistente():
    with pytest.raises(KeyError):
        indice_literal({(0, 2): 0}, 0, 5)


# ===== reduce ==============================================================

def test_reduce_igualdad_simple():
    """x0 = x1 deja una sola variable superviviente."""
    mapeo, restantes, _ = reduce(2, [(0, 0, 1, 0)], [], [])
    assert restantes == []
    assert mapeo[0][0] == "var" and mapeo[1][0] == "var"
    assert mapeo[0][1] == mapeo[1][1]        # misma raíz
    assert mapeo[0][2] == mapeo[1][2]        # mismo signo


def test_reduce_complementarios():
    """X_0 = 1 - X_1 deja signos opuestos."""
    mapeo, _, _ = reduce(2, [(0, 0, 1, 1)], [], [])
    assert mapeo[0][1] == mapeo[1][1]
    assert mapeo[0][2] != mapeo[1][2]


def test_reduce_fijacion():
    mapeo, _, _ = reduce(2, [], [], [(0, 0, 1)])
    assert mapeo[0] == ("const", 1)
    assert mapeo[1][0] == "var"


def test_reduce_paridad_colapsa_a_igualdad():
    """Si dos de los cuatro literales ya están ligados, el soporte baja a 2."""
    igualdades = [(0, 0, 1, 0)]                    # x0 = x1
    paridades = [[(0, 0), (1, 0), (2, 0), (3, 0)]]
    mapeo, restantes, _ = reduce(4, igualdades, paridades, [])
    assert restantes == []                          # colapsó
    assert mapeo[2][1] == mapeo[3][1]               # x2 y x3 ligados


def test_reduce_paridad_se_conserva():
    """Cuatro literales independientes: la paridad sobrevive."""
    _, restantes, _ = reduce(4, [], [[(i, 0) for i in range(4)]], [])
    assert len(restantes) == 1
    assert len(restantes[0]) == 4


def test_reduce_infactible():
    """x0 = x1 y a la vez X_0 = 1 - X_1."""
    with pytest.raises(Infactible):
        reduce(2, [(0, 0, 1, 0), (0, 0, 1, 1)], [], [])


def test_reduce_paridad_infactible():
    """Los cuatro fijados con suma impar."""
    fijaciones = [(0, 0, 1), (1, 0, 1), (2, 0, 1), (3, 0, 0)]
    paridades = [[(i, 0) for i in range(4)]]
    with pytest.raises(Infactible):
        reduce(4, [], paridades, fijaciones)


# ===== expande =============================================================

def test_expande_ida_y_vuelta():
    pares = [(0, 1), (0, 2), (1, 2)]
    mapeo, _, _ = reduce(3, [(0, 0, 1, 0)], [], [(2, 0, 1)])
    raiz = mapeo[0][1]
    vals = expande(mapeo, pares, {raiz: 1})
    assert set(vals) == set(pares)
    assert vals[pares[0]] == vals[pares[1]]     # estaban ligados
    assert vals[pares[2]] == 1                  # fijado


def test_expande_falta_una_raiz():
    pares = [(0, 1), (0, 2)]
    mapeo, _, _ = reduce(2, [], [], [])
    with pytest.raises(KeyError):
        expande(mapeo, pares, {})


# ===== la prueba de fondo: el presolve preserva las soluciones =============

def monta(nombre):
    fi = read_fold(TOY / f"{nombre}.fold")
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
    return fi, pares, trip, cr, tt, tto, oo


def satisface(valores, pares, trip, cr, tt, tto, oo):
    """Comprueba las cinco familias sobre una asignación completa."""
    idx = {p: i for i, p in enumerate(pares)}

    def X(u, v):
        i, s = indice_literal(idx, u, v)
        return valores[i] ^ s

    for f, g, h in trip:
        s = X(f, g) + X(g, h) + X(h, f)
        if not (1 <= s <= 2):
            return False
    for f, g, t in cr:
        if X(f, t) != X(g, t):
            return False
    for familia in (tto, oo):
        for (u1, v1), (u2, v2) in familia:
            if X(u1, v1) != X(u2, v2):
                return False
    for cuatro in tt:
        if sum(X(u, v) for u, v in cuatro) % 2 != 0:
            return False
    return True


@pytest.mark.parametrize("nombre", ["map_1x3", "waterbomb4"])
def test_presolve_preserva_las_soluciones(nombre):
    """
    Enumera exhaustivamente el sistema original y el reducido, expande
    este último, y comprueba que los conjuntos de soluciones coinciden.
    Es la única forma segura de validar el presolve.
    """
    fi, pares, trip, cr, tt, tto, oo = monta(nombre)
    n = len(pares)
    assert n <= 12, "demasiadas variables para enumerar"

    # soluciones del sistema original
    originales = {v for v in product((0, 1), repeat=n)
                  if satisface(list(v), pares, trip, cr, tt, tto, oo)}

    # sistema reducido
    ig, par, fij = literales_a_ecuaciones(pares, cr, tto, oo, tt)
    mapeo, restantes, _ = reduce(n, ig, par, fij)
    raices = sorted({e[1] for e in mapeo.values() if e[0] == "var"})

    # enumerar sobre las supervivientes y expandir
    reconstruidas = set()
    for combo in product((0, 1), repeat=len(raices)):
        beta = dict(zip(raices, combo))
        # las paridades restantes hay que comprobarlas aparte
        ok = all(sum(beta[r] ^ p for r, p in lits) % 2 == 0
                 for lits in restantes)
        if not ok:
            continue
        vals = expande(mapeo, pares, beta)
        completo = [vals[p] for p in pares]
        if satisface(completo, pares, trip, cr, tt, tto, oo):
            reconstruidas.add(tuple(completo))

    assert reconstruidas == originales


# ===== componentes =========================================================

def test_una_componente_en_waterbomb4():
    """Las cuatro caras están apiladas: todo queda conectado."""
    fi, pares, trip, cr, tt, tto, oo = monta("waterbomb4")
    ig, par, fij = literales_a_ecuaciones(pares, cr, tto, oo, tt)
    mapeo, restantes, _ = reduce(len(pares), ig, par, fij)
    comps = componentes(mapeo, restantes, trip, pares)
    assert len(comps) == 1


def test_componentes_particionan():
    fi, pares, trip, cr, tt, tto, oo = monta("waterbomb8")
    ig, par, fij = literales_a_ecuaciones(pares, cr, tto, oo, tt)
    mapeo, restantes, _ = reduce(len(pares), ig, par, fij)
    comps = componentes(mapeo, restantes, trip, pares)
    todas = [r for c in comps for r in c]
    assert len(todas) == len(set(todas))
    supervivientes = {e[1] for e in mapeo.values() if e[0] == "var"}
    assert set(todas) == supervivientes
    
    
    
    
    
    
# Ejecutar: !pytest "C:\Users\Trabajador\Desktop\origami\flatfold-milp"

    