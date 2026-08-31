# -*- coding: utf-8 -*-
"""
Genera las instancias de juguete en instances/toy/.
Se ejecuta a mano cuando haga falta rehacerlas.
"""
import json
from pathlib import Path

TOY = Path(__file__).resolve().parent / "toy"


def mapa_1xn(n):
    """Tira de n celdas con n-1 pliegues verticales, alternando M y V."""
    abajo = [[i, 0] for i in range(n + 1)]
    arriba = [[i, 1] for i in range(n, -1, -1)]
    verts = abajo + arriba

    borde = [[i, i + 1] for i in range(2 * n + 1)] + [[2 * n + 1, 0]]
    pliegues = [[i, 2 * n + 1 - i] for i in range(1, n)]

    asign = (["B"] * len(borde)
             + ["M" if k % 2 == 0 else "V" for k in range(len(pliegues))])

    return {"file_spec": 1.1,
            "vertices_coords": verts,
            "edges_vertices": borde + pliegues,
            "edges_assignment": asign}


def escribe(nombre, patron):
    TOY.mkdir(parents=True, exist_ok=True)
    ruta = TOY / f"{nombre}.fold"
    with open(ruta, "w") as f:
        json.dump(patron, f, indent=2)
    print(f"  escrito {ruta.name}: "
          f"{len(patron['vertices_coords'])} vértices, "
          f"{len(patron['edges_vertices'])} aristas")


def main():
    print("Generando instancias de juguete...")
    for n in range(2, 9):
        escribe(f"map_1x{n}", mapa_1xn(n))



if __name__ == "__main__":
    main()









