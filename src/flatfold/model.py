# -*- coding: utf-8 -*-
import gurobipy as gb
from collections import defaultdict



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


def añade_restricciones(m, x, restricciones):
    """
    Input:  el modelo, las variables, y las listas de la etapa (5).
    Output: nada; modifica el modelo.

    PENDIENTE: la etapa (5) todavía no genera restricciones.
    Cuando las genere, aquí van las cinco familias, cada una en su
    propia función para poder activarlas por separado:
        - transitividad sobre tripletas que cubren un punto común
        - cruces: una cara que atraviesa la imagen de una arista
        - tres familias de solapamientos colineales
    """
    pass
    
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
        dual = grafo_dual(G, caras, caras.index(exterior))

        if len(dual) == 0:
            plegadas, pares = {}, []
        else:
            iso = plegar(G, dual, min(dual.nodes))
            plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
            pares = pares_comparables(plegadas)

        print(f"  caras interiores: {len(plegadas)}   variables: {len(pares)}")

        # --- etapa 6 ---
        m, x = construye_modelo(pares)
        estado, valores = resuelve(m, x)

        if valores is None:
            print(f"  sin solución (estado {estado})")
            continue

        print(f"  solución: {valores}")
        print(f"  capas (de arriba abajo): {orden_de_capas(pares, valores)}")

        # --- recuento de estados ---
        m2, x2 = construye_modelo(pares)
        m2.Params.PoolSearchMode = 2
        m2.Params.PoolSolutions = 1000
        m2.Params.LogToConsole = 0
        m2.optimize()
        print(f"  órdenes de capas admitidos por el modelo: {m2.SolCount}")


if __name__ == "__main__":
    main()
