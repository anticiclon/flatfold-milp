# -*- coding: utf-8 -*-

from math import lcm
from fractions import Fraction


Point   = tuple[int, int]


def escala_coordenadas(coords):
    """
    Multiplica todas las coordenadas por un mismo factor s, elegido como el 
    número más pequeño que consigue dos cosas a la vez: que todas las 
    coordenadas queden enteras, y que además queden pares
    Input:
        coords: lista de pares de Fraction. 
    Output:
        Devuelve (lista de Point, s)."""
    # denoms = [c.denominator for par in coords for c in par]
    denominadores = []
    for coor in coords:
        denominadores.append(coor[0].denominator)
        denominadores.append(coor[1].denominator)
    # 2*mínimo común múltiplo
    # operador * desempaqueta la lista en argumentos separados
    s = 2 * lcm(*denominadores)
    vertices = [(int(x * s), int(y * s)) for x, y in coords]
    return vertices, s
    




direcciones = {
    ( 1,  0): 0,   # Este
    ( 1,  1): 1,   # Noreste
    ( 0,  1): 2,   # Norte
    (-1,  1): 3,   # Noroeste
    (-1,  0): 4,   # Oeste
    (-1, -1): 5,   # Suroeste
    ( 0, -1): 6,   # Sur
    ( 1, -1): 7,   # Sureste
}

 

def signo(n):
    return (n > 0) - (n < 0)





def direccion(G, u, v):
    """Da la dirección u -> v."""
    xu, yu = G.nodes[u]["coor"]
    xv, yv = G.nodes[v]["coor"]
    return direcciones[(signo(xv - xu), signo(yv - yu))]



def pendiente_admisible(p, q):
    """
    Input: dos puntos. Output: booleano.
    Comprueba que el segmento es horizontal, vertical, diagonal o antidiagonal
    """
    dx, dy = q[0] - p[0], q[1] - p[1]
    if dx == 0 and dy == 0:
        return False
    return dx == 0 or dy == 0 or abs(dx) == abs(dy)
        
        

def reflexion(p, q):
    """
    Dado dos puntos que delimitan una recta, la función devuelve la 
    transformación (reflexión más traslación) a la que se somete cualquier 
    punto al reflejarse respecto a esa recta
    """
    vector = (q[0]-p[0], q[1]-p[1])
    direccion = (int(signo(vector[0])), int(signo(vector[1])))
    # Recta vertical
    if direccion in ((0, 1), (0, -1)):
        A = (-1, 0, 0, 1)
        t = (2*p[0], 0)
    # Recta horizontal
    elif direccion in ((1, 0), (-1, 0)):
        A = (1, 0, 0, -1)
        t = (0, 2*p[1])
    # Recta diagonal
    elif direccion in ((1, 1), (-1, -1)):
        A = (0, 1, 1, 0)
        t = (-(p[1]-p[0]), p[1]-p[0])
    # Recta antidiagonal
    elif direccion in ((1, -1), (-1, 1)):
        A = (0, -1, -1, 0)
        t = (p[0]+p[1], p[0]+p[1])
    else:
        raise ValueError("Pendiente no admitida")
    return A, t

def aplica_isometria(isometria, punto):
    """Aplica una isometría a un punto: A p + t."""
    (a, b, c, d), (t0, t1) = isometria
    x, y = punto
    return (a*x + b*y + t0, c*x + d*y + t1)


def compon_isometrias(f, g):
    """
    Input: dos isometrías. 
    Output: la isometría que resulta de aplicar primero g y después f.
    # Es la operación que encadena reflexiones a lo largo del giro alrededor de
    un vértice. El orden importa porque las reflexiones no conmutan.
    """
    (a_f, b_f, c_f, d_f), (t_f0, t_f1) = f
    (a_g, b_g, c_g, d_g), (t_g0, t_g1) = g
    
    A = (a_f*a_g + b_f*c_g, a_f*b_g + b_f*d_g, 
          c_f*a_g + d_f*c_g, c_f*b_g + d_f*d_g)
    
    t = (a_f*t_g0 + b_f*t_g1 + t_f0, c_f*t_g0 + d_f*t_g1 + t_f1)
    
    return A, t
    

def es_identidad(isometria):
    """
    Input: una isometría. 
    Output: booleano.
    Es el veredicto de Kawasaki: la composición completa alrededor del vértice 
    tiene que dar esto.
    """
    A, t = isometria
    return A == (1, 0, 0, 1) and t == (0, 0)



def normaliza_a_reticula(coords, tol=1e-9, max_den=10000):
    """
    Input:  coords — lista de pares de números (con ruido de coma flotante)
            tol    — distancia máxima admitida a un punto de la retícula
            pasos  — denominadores candidatos, de menos a más fino
    Output: (lista de pares de Fraction exactos, k) donde la retícula es 1/k.

    Los ficheros producidos por herramientas que trabajan en coma flotante
    traen coordenadas como 0.9999999999999998 en lugar de 1. Esta función
    recupera la retícula original: busca el 1/k más grueso que explica todas
    las coordenadas dentro de la tolerancia, y ajusta cada una a él.

    Es el ÚNICO redondeo del proyecto. Ocurre una sola vez, en la frontera
    del sistema, y se comprueba que cada punto está realmente cerca de la
    retícula: si no lo está, la instancia se rechaza en lugar de forzarse.
    """
    fracs = [[Fraction(float(c)).limit_denominator(max_den) for c in par]
             for par in coords]

    k = 1
    for par in fracs:
        for f in par:
            k = lcm(k, f.denominator)
    if k > max_den:
        raise ValueError(f"retícula demasiado fina: 1/{k}")

    for par, orig in zip(fracs, coords):
        for f, c in zip(par, orig):
            if abs(float(f) - float(c)) > tol:
                raise ValueError(
                    f"la coordenada {c!r} no cae en una retícula racional; "
                    f"la mejor aproximación es {f} y dista {abs(float(f)-float(c)):.2e}")

    return [tuple(par) for par in fracs], k






    
def sigma(isometria):
    """
    Input:  una isometría (A, t).
    Output: +1 si conserva la orientación, -1 si la invierte.

    Es el determinante de la parte lineal. Con A ortogonal solo puede
    valer +1 o -1: +1 son giros y traslaciones, -1 reflexiones.

    En el plegado equivale a la paridad del número de pliegues cruzados
    desde la cara de referencia, y dice si la cara quedó del derecho o
    del revés. De ahí sale la etiqueta monte/valle.
    """
    (a, b, c, d), _ = isometria
    det = a * d - b * c
    assert det in (1, -1), f"parte lineal no ortogonal: {isometria}"
    return det
    
    
    
    
    
    
    
    
def indice_literal(indice_pares, u, v):
    """
    Input:  indice_pares — {(f,g): i} con f < g;  u, v — dos caras.
    Output: (i, s); s = 0 significa x_i, s = 1 significa 1 - x_i.
    """
    if u == v:
        raise ValueError(f"literal sobre la misma cara: {u}")
    if u < v:
        return indice_pares[(u, v)], 0
    return indice_pares[(v, u)], 1
    
    
    







