# -*- coding: utf-8 -*-
import gurobipy as gb
from collections import defaultdict
from flatfold.pattern import aristas_caras
from flatfold.incidences import cruces, tripletas, colineales



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




def añade_restricciones(m, x, indice_pares, tripletas, cruces, taco_taco,
                        taco_tortilla, tortilla_tortilla, blb=None,
                        con_transitividad=True, con_cruces=True,
                        con_colineales=True, con_blb=False):
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
    if con_blb and blb:
        añade_blb(m, x, indice_pares, blb)
        
        
    
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




def analiza(ruta):
    """
    Tubería completa: de fichero FOLD a datos del modelo.
    Devuelve (pares, tripletas, cruces, taco_taco, taco_tortilla,
              tortilla_tortilla).
    """
    from flatfold.io_fold import read_fold
    from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                                  caras_definitivas, grafo_dual, plegar,
                                  cara_plegada, verifica_consistencia)
    from flatfold.validate import valida_pendientes, valida_kawasaki
    from flatfold.geometry import es_convexa
    from flatfold.incidences import pares_comparables

    fi = read_fold(ruta)
    valida_pendientes(fi)

    G = crea_grafo_dirigido(fi)
    valida_kawasaki(G)

    caras = caras_grafo(G)
    _, exterior = caras_definitivas(G, caras)
    idx_ext = caras.index(exterior)
    dual = grafo_dual(G, caras, idx_ext)
    dic_aristas = aristas_caras(G, caras)

    if len(dual) == 0:
        return [], [], [], [], [], []

    iso = plegar(G, dual, min(dual.nodes))
    verifica_consistencia(G, dual, iso)

    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
    for f, p in plegadas.items():
        if not es_convexa(p):
            raise ValueError(f"la cara {f} no es convexa: {p}")

    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    cr = cruces(G, caras, plegadas, iso, dic_aristas, idx_ext, set(pares))
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic_aristas, idx_ext)
    return pares, trip, cr, tt, tto, oo



def construye_reducido(datos, nombre="componente", minimizar=False):
    """
    Input:  datos — un diccionario de reparte_por_componentes:
                    {"raices", "tripletas", "paridades", "terminos"}
            minimizar — si hay objetivo, lo pone; si no, factibilidad pura
    Output: (modelo, beta) con beta indexado por raíz.

    Una binaria por raíz de la componente. Las tripletas imponen que los
    tres literales no valgan lo mismo; las paridades, que la suma de sus
    literales sea par (o impar, según la constante que arrastren).
    """
    m = gb.Model(nombre)
    raices = sorted(datos["raices"])
    beta = m.addVars(raices, vtype=gb.GRB.BINARY, name="b")

    def lit(r, s):
        return beta[r] if s == 0 else 1 - beta[r]

    for lits in datos["tripletas"]:
        s = gb.quicksum(lit(r, sg) for r, sg in lits)
        m.addConstr(s >= 1)
        m.addConstr(s <= 2)

    for lits in datos["paridades"]:
        n = len(lits)
        k = m.addVar(vtype=gb.GRB.INTEGER, lb=0, ub=n // 2)
        m.addConstr(gb.quicksum(lit(r, s) for r, s in lits) == 2 * k)

    if minimizar and datos["terminos"]:
        m.setObjective(
            gb.quicksum(c * beta[r] for r, c in datos["terminos"].items()),
            gb.GRB.MINIMIZE)

    m.update()
    return m, beta




# def resuelve_por_componentes(datos_comps, constante=0, modo="factibilidad",
#                              limite_pool=10**7, tiempo=None, silencioso=True):
def resuelve_por_componentes(datos_comps, constante=0, modo="factibilidad",
                             limite_pool=10**7, tiempo=None, silencioso=True,
                             presolve_gurobi=-1):
    """
    Input:  datos_comps — lista de reparte_por_componentes()
            constante — el término independiente de energia_reducida()
            modo — "factibilidad" | "recuento" | "optimizacion"
    Output: según el modo:
            factibilidad  -> (bool, valores)
            recuento      -> (n_estados, None)
            optimizacion  -> (valor_optimo, valores)
            donde valores es {raiz: 0 o 1}, o None si no procede.

    Cada componente se resuelve por separado. Los recuentos se MULTIPLICAN
    y los óptimos se SUMAN. Una componente sin restricciones es una
    variable libre: factible, dos estados, óptimo cero.
    """
    total_estados = 1
    total_valor = constante
    valores = {}
    total_nodes = 0

    for k, datos in enumerate(datos_comps):
        m, beta = construye_reducido(datos, f"comp{k}",
                                     minimizar=(modo == "optimizacion"))
        if silencioso:
            m.Params.OutputFlag = 0
        if tiempo is not None:
            m.Params.TimeLimit = tiempo
        m.Params.IntegralityFocus = 1
        m.Params.MIPGap = 0
        m.Params.Presolve = presolve_gurobi      # <- nueva

        if modo == "recuento":
            m.Params.PoolSearchMode = 2
            m.Params.PoolSolutions = limite_pool
            m.optimize()
            total_nodes += m.NodeCount
            if m.SolCount == 0:
                return 0, None
            if m.SolCount >= limite_pool:
                raise RuntimeError(
                    f"componente {k}: el pool se llenó ({limite_pool})")
            total_estados *= m.SolCount
            continue

        if modo == "factibilidad":
            m.Params.SolutionLimit = 1

        m.optimize()
        # if m.SolCount == 0:
        #     return (False, None) if modo == "factibilidad" else (None, None)
        total_nodes += m.NodeCount
        
        if m.SolCount == 0:
            return (False, None, total_nodes) if modo == "factibilidad" \
                   else (None, None, total_nodes)

        valores.update({r: int(round(v.X)) for r, v in beta.items()})
        if modo == "optimizacion":
            total_valor += m.ObjVal

    if modo == "recuento":
        # return total_estados, None
        return total_estados, None, total_nodes
    if modo == "factibilidad":
        # return True, valores
        return True, valores, total_nodes
    # return total_valor, valores
    return total_valor, valores, total_nodes



# def añade_blb(m, x, indice_pares, blb):
#     """Los dos pliegues que delimitan un mínimo local tienen etiquetas
#     opuestas: m_e1 != m_e2."""
#     for ((u1, v1), s1), ((u2, v2), s2) in blb:
#         l1 = literal(x, indice_pares, u1, v1)
#         l2 = literal(x, indice_pares, u2, v2)
#         if s1:
#             l1 = 1 - l1
#         if s2:
#             l2 = 1 - l2
#         m.addConstr(l1 == 1 - l2)

def añade_blb(m, x, indice_pares, blb):
    """Los dos pliegues que delimitan un mínimo local tienen etiquetas
    opuestas: m_e1 != m_e2."""
    for ((u1, v1), s1), ((u2, v2), s2) in blb:
        l1 = literal(x, u1, v1)
        l2 = literal(x, u2, v2)
        if s1:
            l1 = 1 - l1
        if s2:
            l2 = 1 - l2
        m.addConstr(l1 == 1 - l2)


# def main():
#     from pathlib import Path
#     import gurobipy as gb
#     from flatfold.io_fold import read_fold
#     from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
#                                   caras_definitivas, grafo_dual, plegar,
#                                   cara_plegada)
#     from flatfold.incidences import pares_comparables

#     TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

#     for nombre in ["single_fold", "square",
#                    "map_1x2", "map_1x3", "map_1x4", "map_1x5",
#                    "map_1x6", "map_1x7", "map_1x8",
#                    "waterbomb4", "waterbomb8"]:
#         print("=" * 60)
#         print(nombre)

#         # --- etapas 1 a 5 ---
#         fi = read_fold(TOY / f"{nombre}.fold")
#         G = crea_grafo_dirigido(fi)
#         caras = caras_grafo(G)
#         _, exterior = caras_definitivas(G, caras)
#         idx_ext = caras.index(exterior)
#         dual = grafo_dual(G, caras, idx_ext)
#         dic_aristas = aristas_caras(G, caras)

#         if len(dual) == 0:
#             print("  sin caras interiores")
#             continue

#         iso = plegar(G, dual, min(dual.nodes))
#         plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
#         pares = pares_comparables(plegadas)
#         trip = tripletas(plegadas, pares)
#         cr = cruces(G, caras, plegadas, iso, dic_aristas, idx_ext)
#         tt, tto, oo = colineales(G, caras, plegadas, iso, dic_aristas, idx_ext)

#         print(f"  variables: {len(pares)}   tripletas: {len(trip)}   "
#               f"cruces: {len(cr)}   taco-taco: {len(tt)}")

#         # --- recuento con cada combinación de familias ---
#         for etiqueta, banderas in [
#                 ("sin restricciones", dict(con_transitividad=False,
#                                            con_cruces=False,
#                                            con_colineales=False)),
#                 ("+ transitividad",   dict(con_cruces=False,
#                                            con_colineales=False)),
#                 ("+ cruces",          dict(con_colineales=False)),
#                 ("todas",             dict())]:
#             m, x = construye_modelo(pares)
#             añade_restricciones(m, x, trip, cr, tt, tto, oo, **banderas)
#             print(f"    {etiqueta:20s}: {cuenta_soluciones(m, 5000)} estados")
            
            
            
            
            
            
            
            

# if __name__ == "__main__":
#     main()
