# -*- coding: utf-8 -*-
from fractions import Fraction
import numpy as np
Point   = tuple[int, int]
LineKey = tuple[str, int]


def es_convexa(poligono):
    """
    Determina si un polígono simple es convexo.

    Args:
        poligono: lista de puntos [(x1, y1), (x2, y2), ..., (xn, yn)]

    Returns:
        True si es convexo, False en caso contrario.
    """
    
    if len(poligono) < 3:
        return False
    
    aristas = list(zip(poligono, poligono[1:] + poligono[:1]))
    pares_aristas = list(zip(aristas, aristas[1:] + aristas[:1]))
    
    # signo de referencia
    signo = 0
    
    for i in pares_aristas:
        # Vectores consecutivos
        v1 = (i[0][1][0] - i[0][0][0], i[0][1][1] - i[0][0][1])
        v2 = (i[1][1][0] - i[1][0][0], i[1][1][1] - i[1][0][1])
        # Producto vectorial (componente z)
        cruz = v1[0] * v2[1] - v1[1] * v2[0]
        
        
        if cruz != 0:
            # Guardar el primer signo encontrado
            if signo == 0:
                if cruz > 0:
                    signo = 1
                else:
                    signo = -1
            # Buscar cambios de signo
            elif (cruz > 0 and signo < 0) or (cruz < 0 and signo > 0):
                return False
            
    return True




def caja(poligono):
    """
    Input:
        poligono: lista de puntos [(x, y), ...]

    Output:
        (min_x, max_x,
         min_y, max_y,
         min_s, max_s,
         min_d, max_d)

    donde:
        s = x + y
        d = x - y
    """
    xs = [i[0] for i in poligono]
    ys = [i[1] for i in poligono]
    suma = [i[0] + i[1] for i in poligono]
    resta = [i[0] - i[1] for i in poligono]

    return (
        min(xs), max(xs),
        min(ys), max(ys),
        min(suma), max(suma),
        min(resta), max(resta)
    )





def interseca_cajas(cajas):
    """
    Input:
        cajas: lista de cajas, donde cada caja tiene el formato
        (min_x, max_x,
         min_y, max_y,
         min_s, max_s,
         min_d, max_d)

    Output:
        caja correspondiente a la intersección.
    """

    inter = list(cajas[0])

    for c in cajas[1:]:
        inter[0] = max(inter[0], c[0])  # min_x
        inter[1] = min(inter[1], c[1])  # max_x

        inter[2] = max(inter[2], c[2])  # min_y
        inter[3] = min(inter[3], c[3])  # max_y

        inter[4] = max(inter[4], c[4])  # min_s
        inter[5] = min(inter[5], c[5])  # max_s

        inter[6] = max(inter[6], c[6])  # min_d
        inter[7] = min(inter[7], c[7])  # max_d

    return tuple(inter)




def recorta(poligono, tipo, umbral, sentido):
    """
    Recorta un polígono convexo por un semiplano.

    Input:
        poligono: lista de puntos en orden antihorario
        tipo:     "x", "y", "s" (x+y) o "d" (x-y)
        umbral:   el valor de la cota
        sentido:  +1 para  forma <= umbral,  -1 para  forma >= umbral

    Output: el polígono recortado (lista de puntos, quizá vacía)
    """
    def forma(p):
        """
        Evalúa la expresión que toca según el tipo. Los cuatro tipos 
        corresponden a las cuatro direcciones: "x" para rectas verticales, "y" 
        para horizontales, "s" para antidiagonales y "d" para diagonales.
        """
        x, y = p
        return {"x": x, "y": y, "s": x + y, "d": x - y}[tipo]

    def dentro(p):
        """
        Compara con la cota. El sentido distingue si la desigualdad es ≤ o ≥, 
        porque de cada tipo hay dos: una cota inferior y una superior.
        Los puntos sobre la recta cuentan como dentro, porque las comparaciones
        no son estrictas
        """
        if sentido == 1:
            return (forma(p) <= umbral)
        else:
            return (forma(p) >= umbral)
            
    def corte(p, q):
        """
        Calcula dónde la arista cruza la recta. La expresión es la 
        interpolación lineal habitual: avanzas de p hacia q la fracción num/den
        que hace que la forma alcance el umbral.
        """
        fp, fq = forma(p), forma(q)
        num = umbral - fp
        den = fq - fp
        assert (q[0] - p[0]) * num % den == 0, f"corte no entero en x: {p}->{q}"
        assert (q[1] - p[1]) * num % den == 0, f"corte no entero en y: {p}->{q}"
        return (p[0] + (q[0] - p[0]) * num // den,
                p[1] + (q[1] - p[1]) * num // den)

    salida = []
    n = len(poligono)
    for i in range(n):
        p = poligono[i]
        q = poligono[(i + 1) % n]
        p_dentro, q_dentro = dentro(p), dentro(q)
        if p_dentro and q_dentro:
            salida.append(q)
        elif p_dentro and not q_dentro:
            salida.append(corte(p, q))
        elif not p_dentro and q_dentro:
            salida.append(corte(p, q))
            salida.append(q)
    return salida


def poligono_de_caja(caja_lista):
    """
    Input:  los ocho enteros (min_x, max_x, min_y, max_y,
                              min_s, max_s, min_d, max_d)
    Output: el polígono convexo que describen, como lista de vértices.
            Lista vacía o degenerada si la región no tiene interior.
    """
    #  Componentens
    xmin, xmax, ymin, ymax, smin, smax, dmin, dmax = caja_lista

    # Comprobar que tienen sentido
    if xmin > xmax or ymin > ymax:
        return []

    # Partimos del rectángulo caracterizado por sus vértices en orden 
    # antihorario. Nos lo da las desigualdades x−≤x≤x+,y−≤y≤y+
    p = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)]
    
    # Hay que recortar las siguientes desigualdades restastantes:
    # x+y≥s−, x+y≤s+, x−y≥d−, x−y≤d+
    for tipo, umbral, sentido in [("s", smax, +1), ("s", smin, -1),
                                  ("d", dmax, +1), ("d", dmin, -1)]:
        p = recorta(p, tipo, umbral, sentido)
        if not p:
            return []

    return p



def doble_area_poligono(vertices):
    """
    vertices = [(x0,y0), (x1,y1), ..., (xn-1,yn-1)]

    Los vértices deben estar ordenados alrededor
    del polígono (horario o antihorario).
    """

    n = len(vertices)

    suma = 0

    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]

        suma += x1 * y2 - y1 * x2

    return abs(suma)


def tiene_interior(caja):
    p = poligono_de_caja(caja)
    if len(p) < 3:
        return False
    return doble_area_poligono(p) != 0




def se_solapan(poligono_f, poligono_g):
    """
    Input: dos polígonos ya plegados.
    Output: booleano.
    La caja de cada uno, la intersección de las dos cajas, y tiene_interior.
    """
    caja_f = caja(poligono_f)
    caja_g = caja(poligono_g)
    inter = interseca_cajas([caja_f, caja_g])
    return tiene_interior(inter)




def atraviesa(segmento, poligono):
    """
    Input:  segmento — (p, q), dos puntos
            poligono — lista de vértices
    Output: booleano.

    Cierto si el interior del segmento entra en el interior del polígono.
    Rozar un borde o tocar un vértice no cuenta: las ocho desigualdades se
    aplican de forma ESTRICTA, así que un segmento que yace sobre una recta
    soporte del polígono no lo atraviesa.
    """
    p, q = segmento
    xmin, xmax, ymin, ymax, smin, smax, dmin, dmax = caja(poligono)

    def forma(punto, tipo):
        x, y = punto
        return {"x": x, "y": y, "s": x + y, "d": x - y}[tipo]

    # El punto del segmento en el parámetro t es p + t (q - p), con t en (0,1).
    lo, hi = Fraction(0), Fraction(1)

    for tipo, umbral, sentido in [("x", xmax, +1), ("x", xmin, -1),
                                  ("y", ymax, +1), ("y", ymin, -1),
                                  ("s", smax, +1), ("s", smin, -1),
                                  ("d", dmax, +1), ("d", dmin, -1)]:
        fp = forma(p, tipo)
        fq = forma(q, tipo)
        # valor a lo largo del segmento:  fp + t (fq - fp)
        # se pide  sentido * valor < sentido * umbral
        den = sentido * (fq - fp)
        num = sentido * (umbral - fp)

        if den == 0:
            # la forma es constante sobre el segmento
            if num <= 0:
                return False          # no se cumple en ningún punto
            continue                  # se cumple en todos: no recorta

        t = Fraction(num, den)
        if den > 0:
            hi = min(hi, t)           # se cumple para t < t*
        else:
            lo = max(lo, t)           # se cumple para t > t*

        if lo >= hi:
            return False

    return lo < hi

def recta_soporte(segmento):
    """
    Input:  (p, q), los dos extremos de un segmento.
    Output: (tipo, valor), la clave canónica de la recta que lo contiene:
                ("H", a)  para  y = a
                ("V", b)  para  x = b
                ("D", c)  para  y - x = c      (pendiente +1)
                ("A", d)  para  x + y = d      (pendiente -1)

    Dos segmentos son colineales si y solo si tienen la misma clave. La
    comparación es de tuplas, sin geometría de por medio.
    """
    (px, py), (qx, qy) = segmento
    dx, dy = qx - px, qy - py

    if dx == 0 and dy == 0:
        raise ValueError(f"segmento degenerado: {segmento}")
    if dy == 0:
        return ("H", py)
    if dx == 0:
        return ("V", px)
    if dy == dx:
        return ("D", py - px)
    if dy == -dx:
        return ("A", px + py)
    raise ValueError(f"pendiente no admisible: {segmento}")




def lado(key: LineKey, p: Point) -> int:
    """
    Una recta parte el plano en dos mitades. side_of te dice en cuál de las 
    dos está un punto, o si está justo encima.
    
    Cuál de los dos lados sale +1 es arbitrario. Lo único que importa es que 
    sea consistente: el mismo lado siempre da el mismo signo. Nunca se va usar 
    el valor absoluto ni a interpretar el signo geométricamente; solo se va a 
    comparar signos entre sí.
    """
    tipo_pendiente = key[0]
    valor = key[1]
    if tipo_pendiente == "H":
        # y - Y0
        return np.sign(p[1] - valor)
    elif tipo_pendiente == "V":
        # x - X0
        return np.sign(p[0] - valor)
    elif tipo_pendiente == "D":
        # (y - x) - Y0
        return np.sign(p[1] - p[0] - valor)
    elif tipo_pendiente == "A":
        # (y + x) - Y0
        return np.sign(p[1] + p[0] - valor)
    raise ValueError(f"LineKey inválida: {key}")



def parametro(key, punto):
    """
    Posición de un punto a lo largo de su recta, como un solo entero.

    Se usa x, salvo en las verticales donde x es constante y no distingue
    nada; ahí se usa y. El valor es monótono a lo largo de la recta, que es
    lo único que hace falta para ordenar y comparar intervalos.
    """
    return punto[1] if key[0] == "V" else punto[0]


def solapan_colinealmente(seg1, seg2):
    """
    Input:  dos segmentos.
    Output: el tramo común como par de parámetros (lo, hi), o None.

    Cierto solo si los dos están sobre la misma recta y comparten un tramo
    de longitud positiva. Compartir un extremo no cuenta.
    """
    key = recta_soporte(seg1)
    if key != recta_soporte(seg2):
        return None

    a1, b1 = sorted(parametro(key, p) for p in seg1)
    a2, b2 = sorted(parametro(key, p) for p in seg2)

    lo, hi = max(a1, a2), min(b1, b2)
    return (lo, hi) if lo < hi else None
