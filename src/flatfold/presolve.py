# -*- coding: utf-8 -*-
from collections import defaultdict, Counter


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
                           etiquetas=None, prescritas=None, blb=None):
    """
    Input:  pares — lista de (f,g) con f<g, que define los índices
            cruces, tto, oo, tt — las familias generadas por incidences
            etiquetas, prescritas — opcionales, para asignación prescrita,
            blb — opcionales
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
    
    if blb:
        igualdades.extend(blb_a_igualdades(blb, pares))

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









def diagnostica_tripletas(mapeo, tripletas, pares):
    """
    Cuenta cómo quedan las tripletas al sustituir por raíces, sin decidir
    qué hacer con las degeneradas. Sirve para saber qué casos aparecen
    de verdad antes de escribir el tratamiento.
    """
    idx = {par: i for i, par in enumerate(pares)}
    conteo = Counter()

    for f, g, h in tripletas:
        lits = []
        for u, v in ((f, g), (g, h), (h, f)):
            i, s = indice_literal(idx, u, v)
            e = mapeo[i]
            lits.append(("c", e[1] ^ s) if e[0] == "const"
                        else ("v", e[1], e[2] ^ s))

        n_const = sum(1 for l in lits if l[0] == "c")
        raices = [l[1] for l in lits if l[0] == "v"]
        n_repes = len(raices) - len(set(raices))

        conteo[(n_const, n_repes)] += 1

    return conteo




def reescribe_tripletas(mapeo, tripletas, pares, traza = None):
    """
    Input:  mapeo — el de reduce(), tripletas — [(f,g,h)], pares.
    Output: (vivas, igualdades_nuevas)
            vivas — [[(raiz, signo) x3]] con tres raíces distintas
            igualdades_nuevas — [(r1, s1, r2, s2)] a realimentar en reduce

    La condición 1 <= X1+X2+X3 <= 2 dice que los tres literales no valen
    lo mismo. Al sustituir por raíces puede degenerar:

      - tres raíces distintas   -> sobrevive
      - dos raíces coinciden    -> si esos dos literales son iguales, el
                                   tercero debe diferir (nueva igualdad);
                                   si son complementarios, ya se cumple
      - las tres coinciden      -> trivial, salvo que los tres literales
                                   tengan el mismo signo, que sería
                                   infactible

    Medido sobre 30 millones de tripletas de grids e instagram, nunca
    aparece un literal determinado por el presolve, de modo que el caso
    "no ambos iguales a c" —que no es afín sobre F2— no se da.
    """
    idx = {par: i for i, par in enumerate(pares)}
    vivas, igualdades = [], []
    if traza is None:
        traza = Counter()


    for f, g, h in tripletas:
        lits = []
        for u, v in ((f, g), (g, h), (h, f)):
            i, s = indice_literal(idx, u, v)
            entrada = mapeo[i]
            if entrada[0] == "const":
                raise Infactible(
                    "tripleta con literal determinado: caso no previsto "
                    f"({f}, {g}, {h})")
            lits.append((entrada[1], entrada[2] ^ s))

        raices = [r for r, _ in lits]
        distintas = len(set(raices))

        if distintas == 3:
            vivas.append(lits)
            traza["viva"] += 1

        elif distintas == 2:
            # localizar el par repetido y dejar el tercero aparte
            if raices[0] == raices[1]:
                (ra, sa), (rb, sb), (rc, sc) = lits
            elif raices[0] == raices[2]:
                (ra, sa), (rc, sc), (rb, sb) = lits
            else:
                (rc, sc), (ra, sa), (rb, sb) = lits

            if sa == sb:
                # los dos son el mismo literal: el tercero debe diferir
                # si sa != sb ya difieren entre sí: trivial
                igualdades.append((ra, sa, rc, sc ^ 1))
                traza["dos_iguales -> igualdad"] += 1
            else:
                traza["dos_complementarios -> trivial"] += 1

        else:                       # las tres sobre la misma raíz
            signos = {s for _, s in lits}
            if len(signos) == 1:
                raise Infactible(
                    f"tripleta con tres literales idénticos: ({f}, {g}, {h})")
            # hay al menos uno distinto: trivial
            traza["tres_misma_raiz -> trivial"] += 1

    return vivas, igualdades





def reduce_completo(n_vars, igualdades, paridades, fijaciones,
                    tripletas, pares, max_iter=50):
    """
    Input:  lo mismo que reduce(), más las tripletas en caras y los pares.
    Output: (mapeo, paridades_restantes, tripletas_vivas, n_iter)

    Punto fijo exterior. reduce() colapsa el sistema afín; al reescribir
    las tripletas sobre las raíces resultantes, las que degeneran producen
    igualdades nuevas, que a su vez pueden colapsar más paridades y más
    tripletas. Se itera hasta que no aparecen igualdades nuevas.

    Termina siempre: cada iteración une al menos dos clases del union-find,
    y el número de clases es finito.
    """
    ig = list(igualdades)

    for n_iter in range(1, max_iter + 1):
        mapeo, restantes, _ = reduce(n_vars, ig, paridades, fijaciones)
        vivas, nuevas = reescribe_tripletas(mapeo, tripletas, pares)

        if not nuevas:
            return mapeo, restantes, vivas, n_iter

        ig.extend(nuevas)

    raise RuntimeError(f"el punto fijo no cerró en {max_iter} iteraciones")



def energia_reducida(mapeo, pares, etiquetas, objetivo, pesos=None):
    """
    Input:  mapeo — el de reduce()
            pares — lista de (f,g), para traducir a índices
            etiquetas — {arista: ((u,v), s)} de etiquetas_mv
            objetivo — {arista: 0 o 1}, la asignación deseada
            pesos — {arista: w}, o None para peso 1
    Output: (terminos, constante) con terminos = {raiz: coeficiente}

    Traduce la energía a las variables supervivientes:

        sum_{e: m^_e = 0} w_e * m_e  +  sum_{e: m^_e = 1} w_e * (1 - m_e)

    Cada m_e es un literal sobre una variable x_i. Tras el presolve, x_i
    vale beta_r o su complemento, o es una constante. En el primer caso
    aporta un término lineal en beta_r; en el segundo, una constante.

    Los coeficientes se ACUMULAN: dos pliegues distintos pueden acabar
    sobre la misma raíz.
    """
    idx = {par: i for i, par in enumerate(pares)}
    terminos = defaultdict(int)
    constante = 0

    for arista, deseado in objetivo.items():
        if arista not in etiquetas:
            continue
        w = 1 if pesos is None else pesos.get(arista, 1)
        if w == 0:
            continue

        (u, v), s_orient = etiquetas[arista]
        i, s_par = indice_literal(idx, u, v)
        s = s_par ^ s_orient          # signo total del literal m_e

        entrada = mapeo[i]
        if entrada[0] == "const":
            m_e = entrada[1] ^ s
            # penaliza si m_e difiere de lo deseado
            constante += w * (m_e if deseado == 0 else 1 - m_e)
            continue

        _, raiz, s_raiz = entrada
        # m_e = beta_raiz XOR (s ^ s_raiz)
        t = s ^ s_raiz

        # penalización = w * m_e            si deseado = 0
        #              = w * (1 - m_e)      si deseado = 1
        # y m_e = beta si t = 0, o 1 - beta si t = 1
        if (deseado == 0) == (t == 0):
            terminos[raiz] += w           # w * beta
        else:
            terminos[raiz] -= w           # w * (1 - beta)
            constante += w

    return dict(terminos), constante


def reparte_por_componentes(comps, tripletas_vivas, paridades, terminos):
    """
    Input:  comps — lista de conjuntos de raíces
            tripletas_vivas — [[(raiz, signo) x3]]
            paridades — [[(raiz, signo), ...]] con soporte 3 o 4
            terminos — {raiz: coeficiente} del objetivo
    Output: lista de diccionarios, uno por componente:
            {"raices", "tripletas", "paridades", "terminos"}

    Cada componente se resuelve por separado. Toda restricción debe caer
    entera dentro de una componente: si tuviera raíces de dos componentes
    distintas, la descomposición sería incorrecta y los resultados,
    también. De ahí los assert.
    """
    de_quien = {}
    for k, comp in enumerate(comps):
        for r in comp:
            de_quien[r] = k

    datos = [{"raices": set(comp), "tripletas": [], "paridades": [],
              "terminos": {}} for comp in comps]

    def componente_de(lits, que):
        ks = {de_quien[r] for r, _ in lits}
        assert len(ks) == 1, (
            f"{que} cruza {len(ks)} componentes: {lits}")
        return ks.pop()

    for lits in tripletas_vivas:
        datos[componente_de(lits, "tripleta")]["tripletas"].append(lits)

    for lits in paridades:
        datos[componente_de(lits, "paridad")]["paridades"].append(lits)

    for raiz, coef in terminos.items():
        assert raiz in de_quien, f"término sobre una raíz inexistente: {raiz}"
        datos[de_quien[raiz]]["terminos"][raiz] = coef

    return datos


def blb_a_igualdades(blb, pares):
    """
    Input:  blb — [((par1, s1), (par2, s2))] de restricciones_blb()
            pares — lista de (f,g), para traducir a índices
    Output: [(i1, s1, i2, s2)] con el mismo formato que usa reduce()

    BLB dice m_e1 != m_e2, y el formato de igualdades significa X_1 = X_2.
    Para expresar la desigualdad se COMPLEMENTA uno de los dos signos:
        m_e1 != m_e2   <=>   X_1 = 1 - X_2
    """
    idx = {par: i for i, par in enumerate(pares)}
    salida = []
    for ((u1, v1), s_o1), ((u2, v2), s_o2) in blb:
        i1, s_p1 = indice_literal(idx, u1, v1)
        i2, s_p2 = indice_literal(idx, u2, v2)
        salida.append((i1, s_p1 ^ s_o1, i2, (s_p2 ^ s_o2) ^ 1))
    return salida

