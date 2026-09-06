# -*- coding: utf-8 -*-
import gurobipy as gb
from dataclasses import dataclass
from collections import defaultdict
from itertools import combinations
from flatfold.incidences import cruces, tripletas, colineales, pares_comparables

from flatfold.io_fold import read_fold
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                                  caras_definitivas, grafo_dual, plegar,
                                  cara_plegada, verifica_consistencia,
                                  aristas_caras)
from flatfold.validate import valida_pendientes, valida_kawasaki
from flatfold.subdivision import (necesita_subdivision,
                                      subdivide_vertical, representantes)
from flatfold.geometry import es_convexa




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


# def añade_taco_taco(m, x, taco_taco):
#     """Dos pliegues colineales que abren al mismo lado no pueden
#     entrelazarse: un número par de los cuatro literales vale 1."""
#     for cuatro in taco_taco:
#         k = m.addVar(vtype=gb.GRB.INTEGER, lb=0, ub=2)
#         m.addConstr(gb.quicksum(literal(x, u, v) for u, v in cuatro) == 2 * k)


def añade_taco_taco(m, x, taco_taco, facetas=False):
    """
    Dos pliegues colineales que abren al mismo lado no pueden
    entrelazarse: un número par de los cuatro literales vale 1.

    Dos formulaciones equivalentes en los puntos enteros:
      - con variable de acarreo: sum X_i = 2k, una restricción y una
        variable entera;
      - con las facetas del politopo de paridad: 2^(n-1) desigualdades
        sin variables nuevas, y relajación lineal estrictamente más
        fuerte.
    """
    for cuatro in taco_taco:
        if not facetas:
            k = m.addVar(vtype=gb.GRB.INTEGER, lb=0, ub=len(cuatro) // 2)
            m.addConstr(gb.quicksum(literal(x, u, v)
                                    for u, v in cuatro) == 2 * k)
            continue

        for S in facetas_paridad(len(cuatro)):
            m.addConstr(
                gb.quicksum(
                    (1 - literal(x, u, v)) if i in S else literal(x, u, v)
                    for i, (u, v) in enumerate(cuatro)) >= 1)


def añade_igualdades(m, x, pares_de_pares):
    """Familias taco-tortilla y tortilla-tortilla: igualdad de dos literales."""
    for (u1, v1), (u2, v2) in pares_de_pares:
        m.addConstr(literal(x, u1, v1) == literal(x, u2, v2))




def añade_restricciones(m, x, indice_pares, tripletas, cruces, taco_taco,
                        taco_tortilla, tortilla_tortilla, blb=None,
                        maekawa=None, con_transitividad=True,
                        con_cruces=True, con_colineales=True,
                        con_blb=False, con_maekawa=False,
                        facetas_paridad_flag=False,
                        transitividad_lazy=False):
    """
    Añade las cinco familias al modelo. Las banderas permiten activarlas
    por separado, que es lo que hace falta para medir el efecto de cada una.

    Con transitividad_lazy, las tripletas NO se escriben: se cuelgan del
    modelo para que callback_lazy las separe cuando haga falta. Quien
    llame debe entonces poner m.Params.LazyConstraints = 1 y optimizar
    con m.optimize(callback_lazy).
    """
    if con_transitividad:
        if transitividad_lazy:
            m._x = x
            m._tripletas = tripletas
            m._n_lazy = 0
        else:
            añade_transitividad(m, x, tripletas)
    if con_cruces:
        añade_cruces(m, x, cruces)
        # añade_cruces(m, x, cruces_oo)
    if con_colineales:
        añade_taco_taco(m, x, taco_taco, facetas=facetas_paridad_flag)
        añade_igualdades(m, x, taco_tortilla)
        añade_igualdades(m, x, tortilla_tortilla)
    if con_blb and blb:
        añade_blb(m, x, indice_pares, blb)
    if con_maekawa and maekawa:
        añade_maekawa(m, x, maekawa)

        
    
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




# def analiza(ruta):
#     """
#     Tubería completa: de fichero FOLD a datos del modelo.
#     Devuelve (pares, tripletas, cruces, taco_taco, taco_tortilla,
#               tortilla_tortilla).
#     """
#     # from flatfold.io_fold import read_fold
#     # from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
#     #                               caras_definitivas, grafo_dual, plegar,
#     #                               cara_plegada, verifica_consistencia)
#     # from flatfold.validate import valida_pendientes, valida_kawasaki
#     # from flatfold.geometry import es_convexa
#     # from flatfold.incidences import pares_comparables

#     fi = read_fold(ruta)
#     valida_pendientes(fi)

#     G = crea_grafo_dirigido(fi)
#     valida_kawasaki(G)

#     caras = caras_grafo(G)
#     _, exterior = caras_definitivas(G, caras)
#     idx_ext = caras.index(exterior)
#     dual = grafo_dual(G, caras, idx_ext)
#     dic_aristas = aristas_caras(G, caras)

#     if len(dual) == 0:
#         return [], [], [], [], [], []

#     iso = plegar(G, dual, min(dual.nodes))
#     verifica_consistencia(G, dual, iso)

#     plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
#     for f, p in plegadas.items():
#         if not es_convexa(p):
#             raise ValueError(f"la cara {f} no es convexa: {p}")

#     pares = pares_comparables(plegadas)
#     trip = tripletas(plegadas, pares)
#     cr = cruces(G, caras, plegadas, iso, dic_aristas, idx_ext, set(pares))
#     tt, tto, oo = colineales(G, caras, plegadas, iso, dic_aristas, idx_ext)
#     return pares, trip, cr, tt, tto, oo

def analiza(ruta, subdividir=True):
    """Compatibilidad: solo las seis listas."""
    a = analiza_completo(ruta, subdividir)
    return a.pares, a.trip, a.cr, a.tt, a.tto, a.oo



def construye_reducido(datos, nombre="componente", minimizar=False,
                       facetas=False):
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

    # for lits in datos["paridades"]:
    #     n = len(lits)
    #     k = m.addVar(vtype=gb.GRB.INTEGER, lb=0, ub=n // 2)
    #     m.addConstr(gb.quicksum(lit(r, s) for r, s in lits) == 2 * k)
    for lits in datos["paridades"]:
        if facetas:
            for S in facetas_paridad(len(lits)):
                m.addConstr(
                    gb.quicksum((1 - lit(r, s)) if i in S else lit(r, s)
                                for i, (r, s) in enumerate(lits)) >= 1)
        else:
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
                             presolve_gurobi=-1, facetas=False):
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
                                     minimizar=(modo == "optimizacion"),
                                     facetas=facetas)
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
        
        
        
        


def facetas_paridad(n):
    """
    Input:  n — número de literales de la restricción de paridad.
    Output: lista de tuplas con los índices de cada subconjunto de
            cardinal impar.

    La envoltura convexa de los puntos de {0,1}^n con paridad par se
    describe exactamente con una desigualdad por cada subconjunto S de
    cardinal impar:

        sum_{i en S} (1 - X_i)  +  sum_{i fuera de S} X_i  >=  1

    Cada una excluye el punto en que los de S valen 1 y el resto 0, que
    son justamente los 2^(n-1) puntos de paridad impar.

    Son facetas del politopo, así que la relajación lineal es la
    envoltura convexa exacta: más fuerte que la formulación con variable
    de acarreo, y sin variables añadidas.
    """
    return [s for k in range(1, n + 1, 2)
            for s in combinations(range(n), k)]
        


def separa_tripletas(valores, tripletas, tol=1e-6):
    """
    Input:  valores — {(f,g): valor} de la solución candidata
            tripletas — [(f,g,h)]
    Output: la lista de tripletas violadas.

    La condición es 1 <= X_fg + X_gh + X_hf <= 2: los tres literales no
    pueden valer lo mismo. Se viola cuando la suma es 0 o 3, es decir
    cuando las tres comparaciones forman un ciclo.

    Es lógica pura, sin Gurobi: el callback solo la envuelve, y así se
    puede probar con valores puestos a mano.
    """
    def X(u, v):
        return valores[(u, v)] if u < v else 1 - valores[(v, u)]

    violadas = []
    for f, g, h in tripletas:
        s = X(f, g) + X(g, h) + X(h, f)
        if s < 1 - tol or s > 2 + tol:
            violadas.append((f, g, h))
    return violadas
        

def callback_lazy(m, where):
    """
    Callback de separación para la transitividad.

    Se invoca con (modelo, punto del algoritmo). En MIPSOL Gurobi ha
    encontrado una solución entera candidata: se obtienen sus valores con
    cbGetSolution, se buscan las tripletas violadas y se añaden con
    cbLazy.

    Los datos que necesita se cuelgan del modelo antes de optimizar:
        m._x          el diccionario de variables
        m._tripletas  la lista de (f,g,h)
        m._n_lazy     contador de restricciones añadidas

    Requiere m.Params.LazyConstraints = 1.
    """
    if where != gb.GRB.Callback.MIPSOL:
        return

    valores = {par: v for par, v in
               zip(m._x.keys(), m.cbGetSolution(list(m._x.values())))}

    for f, g, h in separa_tripletas(valores, m._tripletas):
        s = (literal(m._x, f, g) + literal(m._x, g, h) + literal(m._x, h, f))
        m.cbLazy(s >= 1)
        m.cbLazy(s <= 2)
        m._n_lazy += 2
        
        
def resuelve_con_lazy(m):
    """Optimiza con separación perezosa, poniendo el parámetro obligatorio."""
    assert hasattr(m, "_tripletas"), "el modelo no tiene tripletas colgadas"
    m.Params.LazyConstraints = 1
    m.optimize(callback_lazy)
    return m

def añade_maekawa(m, x, maekawa):
    """
    Teorema de Maekawa: en cada vértice interior, |M - V| = 2.

    Con m_e = 1 para monte y 0 para valle, la suma de (2*m_e - 1) sobre
    los pliegues del vértice vale M - V. Y 4*s - 2 vale -2 o +2 según la
    binaria auxiliar s, de modo que la igualdad

        sum_e (2*m_e - 1) = 4*s_v - 2

    dice exactamente que la diferencia es 2 en valor absoluto.

    Es la única de las cuatro variantes que añade variables: una binaria
    por vértice. No es afín sobre F2, así que no puede alimentar el
    presolve.
    """
    for lits in maekawa:
        s = m.addVar(vtype=gb.GRB.BINARY)
        suma = gb.quicksum(
            2 * (literal(x, u, v) if sg == 0 else 1 - literal(x, u, v)) - 1
            for (u, v), sg in lits)
        m.addConstr(suma == 4 * s - 2)
        
        
# @dataclass(frozen=True)
# class Analisis:
#     """Todo lo que produce la cadena, de fichero FOLD a datos del modelo."""
#     fi: object                  # el FoldInput, subdividido si hizo falta
#     G: object
#     caras: list
#     idx_ext: int
#     dic_aristas: dict
#     iso: dict
#     plegadas: dict
#     pares: list
#     trip: list
#     cr: list
#     tt: list
#     tto: list
#     oo: list
#     representantes: dict        # {arista_del_fichero: arista_actual}
#     subdividido: bool
@dataclass(frozen=True)
class Analisis:
    fi: object                  # el FoldInput, subdividido si hizo falta
    fi_original: object         # el del fichero, antes de subdividir
    G: object
    caras: list
    idx_ext: int
    dic_aristas: dict
    iso: dict
    plegadas: dict
    pares: list
    trip: list
    cr: list          # cruces con taco: las dos caras de la arista se solapan
    cr_oo: list       # cruces sin taco: tortilla-tortilla
    tt: list
    tto: list
    oo: list
    origen: dict                # {arista_actual: arista_del_fichero}
    representantes: dict        # {arista_del_fichero: arista_actual}
    subdividido: bool


def analiza_completo(ruta, subdividir=True):
    """
    Como analiza(), pero devuelve toda la cadena en lugar de solo las
    seis listas de restricciones. Es lo que necesita quien vaya a
    calcular etiquetas o energía, porque tras subdividir el FoldInput
    cambia y los índices de arista con él.
    """

    fi0 = read_fold(ruta)
    valida_pendientes(fi0)
    G = crea_grafo_dirigido(fi0)
    valida_kawasaki(G)
    caras = caras_grafo(G)
    _, exterior = caras_definitivas(G, caras)
    idx_ext = caras.index(exterior)

    fi = fi0
    subdividido = False
    origen = {i: i for i in range(len(fi0.aristas))}

    if necesita_subdivision(G, caras, idx_ext):
        if not subdividir:
            raise ValueError("hay caras no convexas y subdividir=False")
        fi, origen = subdivide_vertical(fi0, G, caras, idx_ext)
        subdividido = True
        G = crea_grafo_dirigido(fi)
        valida_kawasaki(G)
        caras = caras_grafo(G)
        _, exterior = caras_definitivas(G, caras)
        idx_ext = caras.index(exterior)
        assert not necesita_subdivision(G, caras, idx_ext)

    rep = representantes(origen, fi.tipos, fi0)
    dic_aristas = aristas_caras(G, caras)
    dual = grafo_dual(G, caras, idx_ext)

    if len(dual) == 0:
        return Analisis(fi, fi0, G, caras, idx_ext, dic_aristas, {}, {},
                        [], [], [], [], [], [], [], origen, rep, subdividido)

    iso = plegar(G, dual, min(dual.nodes))
    verifica_consistencia(G, dual, iso)
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}

    pares = pares_comparables(plegadas)
    trip = tripletas(plegadas, pares)
    # cr = cruces(G, caras, plegadas, iso, dic_aristas, idx_ext, set(pares))
    cr, cr_oo = cruces(G, caras, plegadas, iso, dic_aristas, idx_ext,
                   set(pares))
    tt, tto, oo = colineales(G, caras, plegadas, iso, dic_aristas, idx_ext)

    return Analisis(fi, fi0, G, caras, idx_ext, dic_aristas, iso, plegadas,
                    pares, trip, cr, cr_oo, tt, tto, oo, origen, rep,
                    subdividido)

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
