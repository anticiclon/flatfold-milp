# -*- coding: utf-8 -*-
from flatfold.geometry import es_convexa
from flatfold.io_fold import Kind
from flatfold.io_fold import FoldInput


def necesita_subdivision(G, caras, idx_exterior):
    """
    Input:  el grafo, la lista de caras y el índice de la exterior.
    Output: la lista de índices de caras interiores no convexas.
            Vacía si el patrón se puede procesar sin subdividir.

    Toda la maquinaria posterior asume convexidad: la representación por
    caja octogonal solo describe exactamente un polígono si es convexo.
    Con una cara en L, la caja sería su envolvente y los predicados
    darían resultados incorrectos en silencio.

    La convexidad se comprueba en el patrón SIN plegar. Da lo mismo: una
    isometría no puede volver convexa una cara que no lo era.
    """
    no_convexas = []
    for i, cara in enumerate(caras):
        if i == idx_exterior:
            continue
        poligono = [G.nodes[a[0]]["coor"] for a in cara]
        if not es_convexa(poligono):
            no_convexas.append(i)
    return no_convexas




def abscisas_de_corte(fi):
    """
    Input:  el FoldInput.
    Output: lista ordenada de las abscisas distintas de los vértices.

    Define por dónde pasan los cortes verticales del Lema 6.4: una recta
    por cada abscisa de vértice. Es la construcción garantizada, y deja
    todas las caras convertidas en trapecios de lados verticales.

    No es la más económica: cortar solo por los vértices reflejos daría
    menos piezas. Tenerlo aparte permite probar esa variante sin tocar
    el resto.
    """
    return sorted({x for x, _ in fi.vertices})




def franjas(poligono):
    """
    Input:  poligono — lista de vértices de una cara, en coordenadas
            del patrón sin plegar.
    Output: lista de pares (a, b) de abscisas consecutivas, las franjas
            verticales en que queda dividida la cara.

    Son los intervalos entre abscisas consecutivas de vértices, que es
    donde el Lema 6.4 garantiza que la estructura de la sección vertical
    es constante: ninguna arista empieza ni acaba dentro de la franja, y
    las que la cruzan son gráficas de funciones afines disjuntas.
    """
    xs = sorted({x for x, _ in poligono})
    return list(zip(xs, xs[1:]))





def aristas_en_franja(poligono, a, b):
    """
    Input:  poligono — vértices de la cara, en orden cíclico
            (a, b) — una franja de franjas()
    Output: lista de aristas que cruzan la franja, ordenadas de abajo
            arriba, cada una como (x1, y1, x2, y2) con x1 = a y x2 = b.

    Por el Lema 6.4, en el interior de la franja no hay ningún vértice,
    así que toda arista que la toca la cruza de lado a lado como gráfica
    de una función afín, y estas gráficas son disjuntas dos a dos. Por
    tanto están totalmente ordenadas, y el orden no depende de la
    abscisa.

    Las aristas verticales quedan sobre las rectas frontera de la franja
    y no la cruzan: se descartan.

    El orden se establece por el valor en a, desempatando por pendiente:
    si dos aristas comparten el extremo en a, dentro de la franja va
    debajo la de menor pendiente. Así no hace falta evaluar en el punto
    medio, que podría caer fuera de la retícula.
    """
    n = len(poligono)
    cruzan = []

    for i in range(n):
        (x1, y1), (x2, y2) = poligono[i], poligono[(i + 1) % n]
        if x1 == x2:
            continue                        # vertical: no cruza la franja
        if x1 > x2:
            (x1, y1), (x2, y2) = (x2, y2), (x1, y1)
        if x1 > a or x2 < b:
            continue                        # no cubre la franja entera

        # y de la recta en a y en b; la pendiente es 0 o +-1
        pend_num, pend_den = (y2 - y1), (x2 - x1)
        ya = y1 + (a - x1) * pend_num // pend_den
        yb = y1 + (b - x1) * pend_num // pend_den
        cruzan.append((ya, yb))

    cruzan.sort()                           # por ya, y luego por yb
    return [(a, ya, b, yb) for ya, yb in cruzan]


def trapecios_de_cara(poligono):
    """
    Input:  poligono — vértices de una cara, en orden cíclico.
    Output: lista de polígonos convexos, las piezas de la subdivisión.

    Aplica el Lema 6.4: en cada franja vertical, las aristas que la
    cruzan están totalmente ordenadas, y las componentes conexas de la
    cara son las regiones entre pares consecutivos. Cada una es un
    trapecio de lados verticales, degenerado en triángulo cuando las dos
    rectas se encuentran en un extremo de la franja.

    Los vértices salen en orden antihorario.
    """
    piezas = []

    for a, b in franjas(poligono):
        cruzan = aristas_en_franja(poligono, a, b)
        assert len(cruzan) % 2 == 0, (
            f"franja ({a},{b}): {len(cruzan)} aristas, debería ser par")

        for j in range(0, len(cruzan), 2):
            (_, ya_baja, _, yb_baja) = cruzan[j]
            (_, ya_alta, _, yb_alta) = cruzan[j + 1]

            # recorrido antihorario: base de izquierda a derecha, luego
            # el lado derecho hacia arriba, la tapa hacia la izquierda y
            # el lado izquierdo hacia abajo
            vertices = [(a, ya_baja), (b, yb_baja),
                        (b, yb_alta), (a, ya_alta)]

            # eliminar los repetidos: el trapecio degenera en triángulo
            # cuando las dos rectas se tocan en a o en b
            limpio = [v for i, v in enumerate(vertices)
                      if v != vertices[i - 1]]
            if len(limpio) >= 3:
                piezas.append(limpio)

    return piezas



def cortes_de_vertical(vertices, aristas, x):
    """
    Input:  vertices — lista de Point del patrón
            aristas — lista de (u, v), índices en vertices
            x — la abscisa de la recta vertical
    Output: {indice_de_arista: punto_de_corte} para las aristas que la
            recta atraviesa por su interior.

    Solo cuenta el cruce por el interior: si la vertical pasa por un
    extremo de la arista, ese punto ya es un vértice y no hay nada que
    insertar. Las aristas verticales sobre la propia recta tampoco
    cuentan.

    El punto de corte cae en la retícula: una arista de pendiente 0 o
    +-1 con extremos de coordenadas pares corta a una vertical de
    abscisa par en un punto de coordenadas pares.
    """
    cortes = {}
    for i, (u, v) in enumerate(aristas):
        (x1, y1), (x2, y2) = vertices[u], vertices[v]
        if x1 == x2:
            continue                       # vertical: no la atraviesa
        if not (min(x1, x2) < x < max(x1, x2)):
            continue                       # fuera, o justo en un extremo

        num, den = (y2 - y1), (x2 - x1)
        assert num % den == 0 or abs(num) == abs(den), \
            f"pendiente no admisible en la arista {i}"
        y = y1 + (x - x1) * num // den
        cortes[i] = (x, y)
    return cortes





def inserta_cortes(vertices, aristas, tipos, cortes, origen=None):
    """
    Input:  vertices, aristas, tipos — las tres listas del patrón
            cortes — {indice_de_arista: punto}, de cortes_de_vertical
            origen — {indice_de_arista: arista_original}, opcional, para
                     mantener el vínculo con los pliegues del fichero
    Output: (vertices, aristas, tipos, origen) nuevos.

    Cada arista con corte se sustituye por dos, ambas con su mismo tipo.
    Los puntos nuevos se deduplican contra los que ya existen: si el
    mismo punto se generase dos veces con índices distintos, el grafo se
    partiría y la DCEL saldría mal.

    El mapa origen sigue la pista de qué pliegue del fichero produjo cada
    trozo, que es lo que necesita la energía para elegir representantes.
    """
    verts = list(vertices)
    indice = {p: i for i, p in enumerate(verts)}

    def idx_de(p):
        if p not in indice:
            indice[p] = len(verts)
            verts.append(p)
        return indice[p]

    if origen is None:
        origen = {i: i for i in range(len(aristas))}

    nuevas, nuevos_tipos, nuevo_origen = [], [], {}
    for i, (u, v) in enumerate(aristas):
        if i not in cortes:
            nuevo_origen[len(nuevas)] = origen[i]
            nuevas.append((u, v))
            nuevos_tipos.append(tipos[i])
            continue

        w = idx_de(cortes[i])
        for par in ((u, w), (w, v)):
            nuevo_origen[len(nuevas)] = origen[i]
            nuevas.append(par)
            nuevos_tipos.append(tipos[i])

    return verts, nuevas, nuevos_tipos, nuevo_origen




def _parte_por_vertices(p, q, vertices):
    """
    Parte un segmento vertical por los vértices del patrón que caigan en
    su interior. Devuelve la lista de tramos.

    Un lado de pieza puede solaparse parcialmente con una arista ya
    existente en lugar de coincidir con ella: en el hexágono en L, el
    lado derecho de la pieza izquierda va de (2,0) a (2,4), pero solo el
    tramo de (2,0) a (2,2) es nuevo. Partir por los vértices resuelve
    esa comparación.
    """
    x = p[0]
    lo, hi = (p[1], q[1]) if p[1] <= q[1] else (q[1], p[1])
    ys = sorted({y for vx, y in vertices if vx == x and lo <= y <= hi}
                | {lo, hi})
    return [((x, ys[i]), (x, ys[i + 1])) for i in range(len(ys) - 1)]


def segmentos_verticales_nuevos(piezas, vertices, aristas):
    """
    Input:  piezas — los trapecios de trapecios_de_cara()
            vertices, aristas — el patrón
    Output: lista de pares de puntos: los tramos verticales que hay que
            añadir como aristas planas.

    Los lados verticales de los trapecios son los cortes de la
    subdivisión. Cada uno se parte por los vértices que contenga, y se
    conservan los tramos que no coinciden con ninguna arista existente.
    """
    def norm(p, q):
        return (p, q) if p <= q else (q, p)

    existentes = {norm(vertices[u], vertices[v]) for u, v in aristas}

    nuevos = set()
    for pieza in piezas:
        n = len(pieza)
        for i in range(n):
            p, q = pieza[i], pieza[(i + 1) % n]
            if p[0] != q[0]:
                continue
            for s in _parte_por_vertices(p, q, vertices):
                if s not in existentes:
                    nuevos.add(s)

    return sorted(nuevos)



def subdivide_vertical(fi, G, caras, idx_ext):
    """
    Input:  el FoldInput y la topología ya calculada.
    Output: (FoldInput nuevo, origen) donde origen es
            {arista_nueva: arista_del_fichero}.

    Aplica el Lema 6.4 a cada cara no convexa: traza una vertical por
    cada abscisa de sus vértices, inserta los cortes en TODAS las aristas
    del patrón —no solo en las de esa cara, o el grafo dejaría de ser
    planar bien formado— y añade como aristas planas los tramos
    verticales que no existían.

    Si no hay caras no convexas devuelve el patrón intacto, que es lo que
    ocurre en las 401 instancias reales del conjunto de pruebas.
    """
    no_convexas = necesita_subdivision(G, caras, idx_ext)
    if not no_convexas:
        return fi, {i: i for i in range(len(fi.aristas))}

    verts = list(fi.vertices)
    aristas = list(fi.aristas)
    tipos = list(fi.tipos)
    origen = None

    for idx_cara in no_convexas:
        # el polígono de la cara, ANTES de insertar nada: si se
        # recalculasen las abscisas tras cada inserción aparecerían
        # nuevas y el proceso no terminaría
        poligono = [G.nodes[a[0]]["coor"] for a in caras[idx_cara]]

        for x in sorted({p[0] for p in poligono}):
            cortes = cortes_de_vertical(verts, aristas, x)
            verts, aristas, tipos, origen = inserta_cortes(
                verts, aristas, tipos, cortes, origen)

        piezas = trapecios_de_cara(poligono)
        indice = {p: i for i, p in enumerate(verts)}
        for p, q in segmentos_verticales_nuevos(piezas, verts, aristas):
            for punto in (p, q):
                if punto not in indice:
                    indice[punto] = len(verts)
                    verts.append(punto)
            origen[len(aristas)] = None          # el corte no viene del fichero
            aristas.append((indice[p], indice[q]))
            tipos.append(Kind.FLAT)

    nuevo = FoldInput(vertices=verts, aristas=aristas, tipos=tipos,
                      dicc_montana_valle=fi.dicc_montana_valle,
                      escala=fi.escala)
    return nuevo, origen




def representantes(origen, tipos, fi_original):
    """
    Input:  origen — {arista_nueva: arista_original}, de subdivide_vertical
            tipos — los tipos del patrón subdividido
            fi_original — el FoldInput antes de subdividir
    Output: {arista_original: arista_nueva}, un trozo por cada pliegue
            del fichero.

    La energía se define sobre los pliegues del fichero de entrada, donde
    viven la asignación objetivo y los pesos. Al subdividir, un pliegue
    puede quedar partido en varios trozos; se elige uno como
    representante y su literal es el que entra en el objetivo.

    Basta con uno porque la familia de cruces obliga a que todos los
    trozos de un mismo pliegue reciban la misma etiqueta. Es una
    consecuencia del modelo, no una hipótesis: conviene comprobarla sobre
    la solución en lugar de darla por buena.

    Se descartan las aristas planas nuevas (origen None) y las que no
    eran pliegues en el patrón de partida.
    """
    elegido = {}
    for nueva, original in sorted(origen.items()):
        if original is None:
            continue                          # corte de la subdivisión
        if fi_original.tipos[original] != Kind.CREASE:
            continue                          # borde o arista plana
        if tipos[nueva] != Kind.CREASE:
            continue                          # no debería pasar
        elegido.setdefault(original, nueva)
    return elegido













