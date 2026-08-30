# -*- coding: utf-8 -*-


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












