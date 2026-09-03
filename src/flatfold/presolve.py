# -*- coding: utf-8 -*-
from collections import defaultdict


def indice_literal(indice_pares, u, v):
    """
    Input:  indice_pares — {(f,g): i} con f < g;  u, v — dos caras.
    Output: (i, s) donde i es el índice de la variable y s el signo:
            s = 0 significa el literal x_i, s = 1 significa 1 - x_i.

    Es la versión algebraica de literal(): devuelve el índice en lugar
    de la expresión de Gurobi, para poder operar antes de construir
    ningún modelo.
    """
    if u == v:
        raise ValueError(f"literal sobre la misma cara: {u}")
    if u < v:
        return indice_pares[(u, v)], 0
    return indice_pares[(v, u)], 1





class Infactible(Exception):
    """El sistema sobre F2 no tiene solución."""


class UnionFindParidad:
    """
    Union-find donde cada nodo guarda, además del padre, un bit de paridad
    relativa a su raíz:  x_v = x_raiz XOR paridad(v).

    Permite imponer relaciones de la forma  x_u XOR x_v = p, que es la
    forma que toman sobre F2 todas las restricciones del modelo salvo la
    transitividad: una igualdad entre literales del mismo signo es p = 0,
    y entre literales complementarios es p = 1.

    El nodo n (uno más que las variables) es virtual y representa la
    constante 0. Fijar x_v = c es unir v con ese nodo con paridad c, de
    modo que igualdades y fijaciones usan el mismo mecanismo.
    """

    def __init__(self, n):
        self.n = n
        self.BOT = n                       # nodo virtual: la constante 0
        self.padre = list(range(n + 1))
        self.par = [0] * (n + 1)

    def find(self, v):
        """Devuelve (raíz, paridad) con x_v = x_raiz XOR paridad."""
        raiz, p = v, 0
        while self.padre[raiz] != raiz:
            p ^= self.par[raiz]
            raiz = self.padre[raiz]

        # compresión de caminos, recalculando la paridad hasta la raíz
        actual, p_actual = v, p
        while self.padre[actual] != actual:
            siguiente = self.padre[actual]
            p_siguiente = p_actual ^ self.par[actual]
            self.padre[actual] = raiz
            self.par[actual] = p_actual
            actual, p_actual = siguiente, p_siguiente

        return raiz, p

    def union(self, u, v, p):
        """
        Impone x_u XOR x_v = p.
        Devuelve True si se unieron dos clases, False si ya estaban
        ligados de forma coherente. Lanza Infactible si contradice.
        """
        ru, pu = self.find(u)
        rv, pv = self.find(v)
        if ru == rv:
            if pu ^ pv != p:
                raise Infactible(
                    f"contradicción: x_{u} XOR x_{v} = {p} "
                    f"pero ya estaban ligados con {pu ^ pv}")
            return False
        # cuelga rv de ru;  x_rv = x_ru XOR (pu XOR pv XOR p)
        self.padre[rv] = ru
        self.par[rv] = pu ^ pv ^ p
        return True

    def fija(self, v, c):
        """Impone x_v = c."""
        return self.union(v, self.BOT, c)

    def determinado(self, v):
        """Si v está ligado a la constante, devuelve su valor; si no, None."""
        rv, pv = self.find(v)
        rb, pb = self.find(self.BOT)
        return (pv ^ pb) if rv == rb else None




def literales_a_ecuaciones(pares, cruces, tto, oo, tt,
                           etiquetas=None, prescritas=None):
    """
    Input:  pares — lista de (f,g) con f<g, que define los índices
            cruces, tto, oo, tt — las familias generadas por incidences
            etiquetas, prescritas — opcionales, para asignación prescrita
    Output: (igualdades, paridades, fijaciones)
            igualdades: [(i1, s1, i2, s2)]     significa X_1 = X_2
            paridades:  [[(i, s), ...]]        suma par de los literales
            fijaciones: [(i, s, valor)]        significa X = valor

    Traduce cada restricción a su forma algebraica sobre F2. Sobre F2 todo
    el sistema salvo la transitividad es lineal: cruces y colineales con
    tortillas son igualdades, taco-taco es una paridad, y una etiqueta
    prescrita es una fijación.
    """
    idx = {par: i for i, par in enumerate(pares)}

    igualdades = []
    paridades = []
    fijaciones = []

    # (E) cruces: la cara t no puede colarse entre f y g  ->  X_ft = X_gt
    for f, g, t in cruces:
        i1, s1 = indice_literal(idx, f, t)
        i2, s2 = indice_literal(idx, g, t)
        igualdades.append((i1, s1, i2, s2))

    # (S2) y (S3): igualdad entre dos literales
    for familia in (tto, oo):
        for (u1, v1), (u2, v2) in familia:
            i1, s1 = indice_literal(idx, u1, v1)
            i2, s2 = indice_literal(idx, u2, v2)
            igualdades.append((i1, s1, i2, s2))

    # (S1) taco-taco: un número PAR de los cuatro literales vale 1
    for cuatro in tt:
        paridades.append([indice_literal(idx, u, v) for u, v in cuatro])

    # etiquetas prescritas
    if etiquetas and prescritas:
        for arista, valor in prescritas.items():
            if arista not in etiquetas:
                continue
            (u, v), s_orient = etiquetas[arista]
            i, s = indice_literal(idx, u, v)
            fijaciones.append((i, s ^ s_orient, valor))

    return igualdades, paridades, fijaciones



def reduce(n_vars, igualdades, paridades, fijaciones):
    """
    Input:  n_vars — número de variables x_fg
            igualdades — [(i1, s1, i2, s2)], significa X_1 = X_2
            paridades  — [[(i, s), ...]], suma par de los literales
            fijaciones — [(i, s, valor)], significa X = valor
    Output: (mapeo, paridades_restantes, uf)
            mapeo[i] es ("const", v) o ("var", raiz, signo), donde
            signo = 0 significa x_i = x_raiz y 1 que es su complemento.

    Punto fijo: alimenta igualdades y fijaciones al union-find, reescribe
    las paridades en términos de raíces, y las que colapsan producen
    nuevas igualdades o fijaciones que a su vez colapsan otras.

    NO se reduce más allá de esto. Las tripletas se escriben como
    1 <= X + X + X <= 2 sobre literales simples; si una variable pasara
    a ser una suma larga de paridad habría que introducir auxiliares y
    la relajación lineal empeoraría.
    """
    uf = UnionFindParidad(n_vars)

    # X_1 = X_2  <=>  x_i1 XOR x_i2 = s1 XOR s2
    for i1, s1, i2, s2 in igualdades:
        uf.union(i1, i2, s1 ^ s2)

    # X = valor  <=>  x_i = valor XOR s
    for i, s, valor in fijaciones:
        uf.fija(i, valor ^ s)

    pendientes = list(paridades)

    while True:
        nuevas_igualdades = []
        nuevas_fijaciones = []
        siguientes = []
        cambio = False

        for lits in pendientes:
            # reescribir cada literal en términos de su raíz
            soporte = {}          # raiz -> paridad acumulada
            const = 0
            raiz_bot, p_bot = uf.find(uf.BOT)

            for i, s in lits:
                raiz, p = uf.find(i)
                if raiz == raiz_bot:
                    const ^= p ^ p_bot ^ s
                    continue
                # raices repetidas se cancelan:  beta XOR beta = 0
                if raiz in soporte:
                    const ^= soporte.pop(raiz) ^ p ^ s
                else:
                    soporte[raiz] = p ^ s

            n = len(soporte)
            if n != len(lits):
                cambio = True

            if n == 0:
                if const != 0:
                    raise Infactible("paridad sin solución tras reducir")
            elif n == 1:
                raiz, p = soporte.popitem()
                nuevas_fijaciones.append((raiz, p, const))
            elif n == 2:
                (r1, p1), (r2, p2) = soporte.items()
                nuevas_igualdades.append((r1, p1, r2, p2 ^ const))
            else:
                siguientes.append([(r, p) for r, p in soporte.items()])
                # la constante se absorbe en el signo del primer literal
                if const:
                    r, p = siguientes[-1][0]
                    siguientes[-1][0] = (r, p ^ 1)

        for i1, s1, i2, s2 in nuevas_igualdades:
            if uf.union(i1, i2, s1 ^ s2):
                cambio = True
        for i, s, valor in nuevas_fijaciones:
            if uf.fija(i, valor ^ s):
                cambio = True

        pendientes = siguientes
        if not cambio:
            break

    # mapeo final
    raiz_bot, p_bot = uf.find(uf.BOT)
    mapeo = {}
    for i in range(n_vars):
        raiz, p = uf.find(i)
        mapeo[i] = ("const", p ^ p_bot) if raiz == raiz_bot else ("var", raiz, p)

    return mapeo, pendientes, uf



def expande(mapeo, pares, valores_beta):
    """
    Input:  mapeo — el de reduce(): {i: ("const", v) | ("var", raiz, s)}
            pares — lista de (f,g), para traducir índice a par de caras
            valores_beta — {raiz: 0 o 1} con raiz índice de variable
    Output: {(f,g): 0 o 1} completo, como lo consume orden_de_capas.

    Reconstruye todas las variables originales a partir de las
    supervivientes. Las determinadas toman su valor constante; las demás
    valen su raíz o su complemento según el signo.
    """
    valores = {}
    for i, entrada in mapeo.items():
        if entrada[0] == "const":
            valores[pares[i]] = entrada[1]
        else:
            _, raiz, s = entrada
            if raiz not in valores_beta:
                raise KeyError(f"falta el valor de la raíz {raiz} "
                               f"(par {pares[raiz]})")
            valores[pares[i]] = valores_beta[raiz] ^ s
    return valores


def componentes(mapeo, paridades_restantes, tripletas, pares):
    """
    Input:  mapeo — el de reduce()
            paridades_restantes — las que sobrevivieron
            tripletas — [(f,g,h)] en caras, de incidences
            pares — lista de (f,g), para traducir
    Output: lista de conjuntos de raíces (variables supervivientes).

    Grafo cuyos nodos son las variables supervivientes y cuyas aristas
    ligan las que aparecen juntas en alguna restricción. Cada componente
    se resuelve por separado y el número de estados total es el producto
    de los de cada una.

    Las aristas de las tripletas hacen falta AUNQUE la transitividad se
    genere de forma perezosa: ligan variables aunque las restricciones
    no lleguen a escribirse. Omitirlas da una descomposición incorrecta.
    """
    idx = {par: i for i, par in enumerate(pares)}

    def raiz_de(i):
        entrada = mapeo[i]
        return None if entrada[0] == "const" else entrada[1]

    vecinos = defaultdict(set)
    supervivientes = {e[1] for e in mapeo.values() if e[0] == "var"}
    for r in supervivientes:
        vecinos[r]                      # asegura que existe, aunque aislada

    def conecta(raices):
        rs = [r for r in raices if r is not None]
        for a in range(len(rs)):
            for b in range(a + 1, len(rs)):
                if rs[a] != rs[b]:
                    vecinos[rs[a]].add(rs[b])
                    vecinos[rs[b]].add(rs[a])

    for lits in paridades_restantes:
        conecta([r for r, _ in lits])

    for f, g, h in tripletas:
        raices = []
        for u, v in ((f, g), (g, h), (f, h)):
            i, _ = indice_literal(idx, u, v)
            raices.append(raiz_de(i))
        conecta(raices)

    # componentes conexas
    comps, visitadas = [], set()
    for c in sorted(vecinos):
        if c in visitadas:
            continue
        comp, pila = set(), [c]
        visitadas.add(c)
        while pila:
            v = pila.pop()
            comp.add(v)
            for w in vecinos[v]:
                if w not in visitadas:
                    visitadas.add(w)
                    pila.append(w)
        comps.append(comp)
    return comps



