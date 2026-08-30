# -*- coding: utf-8 -*-
import gurobipy as gb
from collections import defaultdict
from flatfold.pattern import aristas_caras
from flatfold.incidences import tripletas, cruces, colineales



def construye_modelo(pares, nombre="flatfold"):
    """
    Input:  pares — lista de pares (f, g) con f < g que se solapan.
    Output: (modelo, x) donde x es {(f,g): variable binaria}.

    x[f,g] = 1 significa que la cara f está por encima de la g.
    Solo hay variable para los pares que se solapan con área positiva;
    los demás nunca compiten por el mismo punto del plano.

    El modelo sale sin restricciones y sin objetivo: es factibilidad pura.
    """
    m = gb.Model(nombre)
    x = m.addVars(pares, vtype=gb.GRB.BINARY, name="x")
    m.update()
    return m, x



def literal(x, u, v):
    """
    Input:  x — {(f,g): variable} con f < g;  u, v — dos caras.
    Output: expresión lineal que vale 1 si u está por encima de v.

    Los pares se guardan una sola vez, con el índice menor primero. Esta
    función traduce cualquier consulta a esa convención:
        u < v  ->  x[u,v]
        v < u  ->  1 - x[v,u]

    La antisimetría queda garantizada por construcción, sin restricción
    que la imponga.
    """
    if u == v:
        raise ValueError(f"literal sobre la misma cara: {u}")
    if u < v:
        return x[u, v]
    return 1 - x[v, u]


def añade_transitividad(m, x, tripletas):
    """Tres caras sobre un mismo punto: su orden debe ser total, no cíclico."""
    for f, g, h in tripletas:
        s = literal(x, f, g) + literal(x, g, h) + literal(x, h, f)
        m.addConstr(s >= 1)
        m.addConstr(s <= 2)


def añade_cruces(m, x, cruces):
    """Una cara que cruza la línea de una arista no puede colarse entre
    los dos trozos unidos por ella."""
    for f, g, t in cruces:
        m.addConstr(literal(x, f, t) == literal(x, g, t))


def añade_taco_taco(m, x, taco_taco):
    """Dos pliegues colineales que abren al mismo lado no pueden
    entrelazarse: un número par de los cuatro literales vale 1."""
    for cuatro in taco_taco:
        k = m.addVar(vtype=gb.GRB.INTEGER, lb=0, ub=2)
        m.addConstr(gb.quicksum(literal(x, u, v) for u, v in cuatro) == 2 * k)


def añade_igualdades(m, x, pares_de_pares):
    """Familias taco-tortilla y tortilla-tortilla: igualdad de dos literales."""
    for (u1, v1), (u2, v2) in pares_de_pares:
        m.addConstr(literal(x, u1, v1) == literal(x, u2, v2))




def añade_restricciones(m, x, tripletas, cruces, taco_taco,
                        taco_tortilla, tortilla_tortilla,
                        con_transitividad=True, con_cruces=True,
                        con_colineales=True):
    """
    Añade las cinco familias al modelo. Las banderas permiten activarlas
    por separado, que es lo que hace falta para medir el efecto de cada una.
    """
    if con_transitividad:
        añade_transitividad(m, x, tripletas)
    if con_cruces:
        añade_cruces(m, x, cruces)
    if con_colineales:
        añade_taco_taco(m, x, taco_taco)
        añade_igualdades(m, x, taco_tortilla)
        añade_igualdades(m, x, tortilla_tortilla)
    
def resuelve(m, x, tiempo_limite=None, semilla=0, silencioso=True):
    """
    Input:  el modelo y sus variables.
    Output: (estado, valores) donde estado es el Model.Status de Gurobi
            y valores es {(f,g): 0 o 1}, o None si no hay solución.

    Sin objetivo, esto resuelve factibilidad pura: devuelve un orden de
    capas válido si existe.
    """
    m.Params.IntegralityFocus = 1     # el modelo es puramente combinatorio
    m.Params.IntFeasTol = 1e-9
    m.Params.MIPGap = 0               # óptimo exacto
    m.Params.Seed = semilla           # reproducibilidad
    if silencioso:
        m.Params.LogToConsole = 0
    if tiempo_limite is not None:
        m.Params.TimeLimit = tiempo_limite

    m.optimize()

    if m.SolCount == 0:
        return m.Status, None

    valores = {k: int(round(v.X)) for k, v in x.items()}
    return m.Status, valores

def cuenta_soluciones(m, limite=1000):
    """Número de órdenes de capas válidos. Solo para instancias pequeñas."""
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = limite
    m.Params.LogToConsole = 0
    m.optimize()
    return m.SolCount
    
    


def orden_de_capas(pares, valores):
    """
    Input:  pares — lista de (f,g) con f < g
            valores — {(f,g): 0 o 1} que devuelve el solver
    Output: lista de grupos, cada uno con sus caras ordenadas de arriba abajo.

    Las caras solo son comparables dentro de un mismo grupo de solapamiento.
    Entre grupos distintos no hay relación de orden, así que se devuelven
    por separado.
    """
    # Grupos: componentes conexas del grafo de solapamiento
    vecinos = defaultdict(set)
    for f, g in pares:
        vecinos[f].add(g)
        vecinos[g].add(f)

    grupos, visitadas = [], set()
    for c in sorted(vecinos):
        if c in visitadas:
            continue
        grupo, pila = [], [c]
        visitadas.add(c)
        while pila:
            v = pila.pop()
            grupo.append(v)
            for w in vecinos[v]:
                if w not in visitadas:
                    visitadas.add(w)
                    pila.append(w)
        grupos.append(grupo)

    # Dentro de cada grupo: cuántas caras hay por encima de cada una
    def encima(u, v):
        return valores[(u, v)] if u < v else 1 - valores[(v, u)]

    resultado = []
    for grupo in grupos:
        clave = {c: sum(1 for o in grupo if o != c and encima(o, c) == 1)
                 for c in grupo}
        resultado.append(sorted(grupo, key=clave.get))
    return resultado




def main():
    from pathlib import Path
    import gurobipy as gb
    from flatfold.io_fold import read_fold
    from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                                  caras_definitivas, grafo_dual, plegar,
                                  cara_plegada)
    from flatfold.incidences import pares_comparables

    TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

    for nombre in ["single_fold", "square", "map_1x3",
                   "waterbomb4", "waterbomb8"]:
        print("=" * 60)
        print(nombre)

        # --- etapas 1 a 5 ---
        fi = read_fold(TOY / f"{nombre}.fold")
        G = crea_grafo_dirigido(fi)
        caras = caras_grafo(G)
        _, exterior = caras_definitivas(G, caras)
        idx_ext = caras.index(exterior)
        dual = grafo_dual(G, caras, idx_ext)
        dic_aristas = aristas_caras(G, caras)

        if len(dual) == 0:
            print("  sin caras interiores")
            continue

        iso = plegar(G, dual, min(dual.nodes))
        plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
        pares = pares_comparables(plegadas)
        trip = tripletas(plegadas, pares)
        cr = cruces(G, caras, plegadas, iso, dic_aristas, idx_ext)
        tt, tto, oo = colineales(G, caras, plegadas, iso, dic_aristas, idx_ext)

        print(f"  variables: {len(pares)}   tripletas: {len(trip)}   "
              f"cruces: {len(cr)}   taco-taco: {len(tt)}")

        # --- recuento con cada combinación de familias ---
        for etiqueta, banderas in [
                ("sin restricciones", dict(con_transitividad=False,
                                           con_cruces=False,
                                           con_colineales=False)),
                ("+ transitividad",   dict(con_cruces=False,
                                           con_colineales=False)),
                ("+ cruces",          dict(con_colineales=False)),
                ("todas",             dict())]:
            m, x = construye_modelo(pares)
            añade_restricciones(m, x, trip, cr, tt, tto, oo, **banderas)
            print(f"    {etiqueta:20s}: {cuenta_soluciones(m)} estados")

if __name__ == "__main__":
    main()
