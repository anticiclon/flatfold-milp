# -*- coding: utf-8 -*-


# from fractions import Fraction
import numpy as np

Point   = tuple[int, int]
Iso     = tuple[int, tuple[int, int]]
LineKey = tuple[str, int]

D4 = [ (1,0,0,1), (0,-1,1,0), (-1,0,0,-1), (0,1,-1,0),
       (1,0,0,-1), (-1,0,0,1), (0,1,1,0), (0,-1,-1,0) ]

IDENTITY: Iso = (0, (0, 0))



def line_key(p: Point, q: Point) -> LineKey:
    """
    Dado dos puntos da una dupla que dice si están en la horizontal,
    la vertical, diagonal o antidiagonal. En la segunda posición de
    la dupla, se da el valor en el eje de abcisas o coordenadas o
    la ordenada en el origen.
    """
    if p == q:
        raise ValueError("p == q")
    elif p[1] == q[1]:
        return ("H", p[1])
    elif p[0] == q[0]:
        return ("V", p[0])
    elif (p[1] - p[0]) == (q[1] - q[0]):
        return ("D", p[1] - p[0])
    elif (p[0] + p[1]) == (q[0] + q[1]):
        return ("A", p[0] + p[1])
    elif abs(p[1] - p[0]) != abs(q[1] - q[0]):
        raise ValueError("pendiente diferent de 0, inf, 1 o -1")



def side_pf(key: LineKey, p: Point) -> int:
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


























def test_line_key():
    casos_ok = [
        (((2, 6), (10, 6)), ("H", 6)),
        (((0, 0), (4, 0)), ("H", 0)),
        (((-4, -2), (8, -2)), ("H", -2)),
        (((4, 0), (4, 12)), ("V", 4)),
        (((-6, 2), (-6, 10)), ("V", -6)),
        (((0, 0), (6, 6)), ("D", 0)),
        (((2, 8), (10, 16)), ("D", 6)),
        (((4, 0), (10, 6)), ("D", -4)),
        (((0, 8), (8, 0)), ("A", 8)),
        (((2, 2), (6, -2)), ("A", 4)),
        (((-4, 4), (4, -4)), ("A", 0)),
        (((10, 6), (2, 6)), ("H", 6)),
        (((0, 6), (4, 6)), ("H", 6)),
        (((20, 6), (30, 6)), ("H", 6)),
        (((0, 4), (6, 10)), ("D", 4)),
    ]

    casos_error = [
        ((0, 0), (6, 2)),
        ((0, 0), (2, 6)),
        ((4, 6), (4, 6)),
    ]

    for (p, q), esperado in casos_ok:
        obtenido = line_key(p, q)
        assert obtenido == esperado, (
            f"ERROR: line_key({p}, {q}) -> "
            f"{obtenido}, esperado {esperado}"
        )

    for p, q in casos_error:
        try:
            line_key(p, q)
            raise AssertionError(
                f"ERROR: line_key({p}, {q}) debería lanzar ValueError"
            )
        except ValueError:
            pass

    print("Todos los tests han pasado.")
    
    
    
    
    
def test_side_pf():
    casos_ok = [
        ((("H", 6), (2, 10)),  1),
        ((("H", 6), (2,  6)),  0),
        ((("H", 6), (2,  0)), -1),

        ((("V", 4), (10, 2)),  1),
        ((("V", 4), ( 4, 2)),  0),
        ((("V", 4), ( 0, 2)), -1),

        ((("D", 0), (2, 8)),   1),
        ((("D", 0), (8, 2)),  -1),
        ((("D", 0), (4, 4)),   0),

        ((("A", 8), (6, 6)),   1),
        ((("A", 8), (0, 0)),  -1),
        ((("A", 8), (4, 4)),   0),

        ((("D", 6), (10, 16)),  0),
        ((("D", 6), (10, 20)),  1),
        ((("D", 6), (10, 12)), -1),
    ]

    casos_error = [
        (("X", 0), (1, 1)),
        (("", 0), (1, 1)),
    ]

    for (key, p), esperado in casos_ok:
        obtenido = side_pf(key, p)
        assert obtenido == esperado, (
            f"ERROR: side_pf({key}, {p}) -> "
            f"{obtenido}, esperado {esperado}"
        )

    for key, p in casos_error:
        try:
            side_pf(key, p)
            raise AssertionError(
                f"ERROR: side_pf({key}, {p}) debería lanzar ValueError"
            )
        except ValueError:
            pass

    print("Todos los tests de side_pf han pasado.")
    
    
    



def main():
    # test_l    ine_key()
    test_side_pf()
     
if __name__ == "__main__":
    main()
     





















