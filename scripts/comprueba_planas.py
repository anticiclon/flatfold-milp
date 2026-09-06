# -*- coding: utf-8 -*-
"""
Inspecciona las cuatro instancias que difieren de la referencia solo en
las familias con tortillas.

La hipótesis es que fueron reetiquetadas —aristas planas convertidas en
pliegues o al revés— después de generar el CSV de referencia. Si el
fichero actual tiene pocas o ninguna arista FLAT mientras el CSV espera
restricciones tortilla, la explicación queda confirmada.
"""
import csv
from collections import Counter
from pathlib import Path

from flatfold.io_fold import read_fold
from flatfold.model import analiza_completo

RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"

CASOS = [
    ("grids", "121_brown_4x4_A99"),
    ("grids", "171_brown_4x4_C30"),
    ("instagram", "086_miwu_5x5_Flippable_Pixel_Grid_1"),
    ("instagram", "087_miwu_5x5_Flippable_Pixel_Grid_2"),
]

CAMPOS = ["vertices", "edges", "faces", "variables", "taco-taco",
          "taco-tortilla", "tortilla-tortilla", "transitivity"]


def referencia(carpeta):
    ruta = EJEMPLOS / f"{carpeta}_data.csv"
    with open(ruta, encoding="utf-8") as f:
        return {r["number"]: r for r in csv.DictReader(f)}


def main():
    for carpeta, nombre in CASOS:
        ruta = EJEMPLOS / carpeta / f"{nombre}.fold"
        if not ruta.exists():
            print(f"{nombre}: no está en {carpeta}")
            continue

        fi = read_fold(ruta)
        tipos = Counter(t.name for t in fi.tipos)
        a = analiza_completo(ruta)
        ref = referencia(carpeta)[nombre.split("_")[0]]

        mio = {
            "vertices": len(a.fi_original.vertices),
            "edges": len(a.fi_original.aristas),
            "faces": len(a.plegadas),
            "variables": len(a.pares),
            "taco-taco": len(a.tt),
            "taco-tortilla": len(a.cr) + len(a.tto),
            "tortilla-tortilla": len(a.cr_oo) + len(a.oo),
            "transitivity": len(a.trip),
        }

        print("=" * 62)
        print(f"{carpeta}/{nombre}")
        print(f"  tipos de arista: {dict(tipos)}")
        print(f"  subdividido: {a.subdividido}")
        print(f"  {'campo':20s} {'nuestro':>10s} {'referencia':>12s}")
        for c in CAMPOS:
            marca = "  <--" if mio[c] != int(ref[c]) else ""
            print(f"  {c:20s} {mio[c]:>10d} {int(ref[c]):>12d}{marca}")


if __name__ == "__main__":
    main()