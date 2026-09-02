# -*- coding: utf-8 -*-

from dataclasses import dataclass
from fractions import Fraction
import json
from enum import Enum
from flatfold.exact import Point, escala_coordenadas, normaliza_a_reticula




# Un Enum (enumeración) es un conjunto cerrado de valores con nombre. 
# Se usa cuando una variable solo puede tomar unos pocos valores fijos 
# y conocidos de antemano.
class Kind(Enum):
    CREASE = "C"
    FLAT = "F"
    BOUNDARY = "B"


@dataclass(frozen=True)
class FoldInput:
    vertices:     list[Point]          # [(0,0), (2,0), (4,0), ...] enteros pares
    aristas:     list[tuple[int,int]] # [(0,1), (1,2), ..., (1,4)] índices
    tipos:     list[Kind]           # [BOUNDARY, ..., CREASE]
    dicc_montana_valle: dict[int, int]       # {6: 1}
    escala:     int                  # 2





def read_fold(path):
    with open(path, encoding="utf-8") as f:
        # Fraction hace que los números decimales del JSON se lean como 
        # fracciones exactas en lugar de como float. Obtienes un objeto 
        # Fraction. Puedes acceder a esos valores así:
        # print(x.numerator)
        # print(x.denominator) 
        data = json.load(f, parse_float = Fraction)
    # vertices , escala = escala_coordenadas(data["vertices_coords"])
    coords, k = normaliza_a_reticula(data["vertices_coords"])
    vertices, escala = escala_coordenadas(coords)
    aristas = [(i, j) for i,j in data["edges_vertices"]]
    #  data.get("clave") devuelve None si no está
    asignaciones = data.get("edges_assignment")
    if asignaciones is None:
        # "U" de unassigned
        asignaciones = ["U"] * len(aristas)
    tipos = []
    dicc_montana_valle = {}
    for idx, letra in enumerate(asignaciones):
        if letra in ("M", "V", "U"):
            tipos.append(Kind.CREASE)
            if letra == "M":
                dicc_montana_valle[idx] = 1
            elif letra == "V":
                dicc_montana_valle[idx] = 0
        elif letra == "F":
            tipos.append(Kind.FLAT)
        elif letra == "B":
            tipos.append(Kind.BOUNDARY)
        else:
            raise ValueError(f"asignación desconocida: {letra!r}")
    return FoldInput(
        vertices = vertices,
        aristas = aristas,
        tipos = tipos,
        dicc_montana_valle = dicc_montana_valle,
        escala = escala)

















    











