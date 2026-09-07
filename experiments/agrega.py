# -*- coding: utf-8 -*-
"""
Concatena las filas sueltas de un experimento en un único CSV.

run_experimentos.py escribe un fichero por ejecución —esa granularidad es
lo que permite paralelizar sin que las tareas se pisen—. Este guion los
junta para el análisis.

Uso:
    python experiments/agrega.py E1
    python experiments/agrega.py E1 --patron "*grids*"
"""
import argparse
import csv
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
RESULTADOS = RAIZ / "experiments" / "results"


def agrega(experimento, patron="*.csv"):
    carpeta = RESULTADOS / experimento
    if not carpeta.is_dir():
        raise SystemExit(f"no existe {carpeta}")

    ficheros = sorted(f for f in carpeta.glob(patron)
                      if f.name != f"{experimento}.csv")
    if not ficheros:
        raise SystemExit(f"no hay ficheros que casen con {patron}")

    filas, campos = [], []
    for f in ficheros:
        with open(f, encoding="utf-8") as h:
            for fila in csv.DictReader(h):
                filas.append(fila)
                for c in fila:
                    if c not in campos:
                        campos.append(c)

    # las filas pueden tener columnas distintas: las rechazadas no llegan
    # a resolverse y no llevan los campos del solver
    destino = RESULTADOS / f"{experimento}.csv"
    with open(destino, "w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=campos, restval="")
        w.writeheader()
        w.writerows(filas)

    print(f"{len(filas)} filas de {len(ficheros)} ficheros -> {destino}")
    return destino


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("experimento")
    p.add_argument("--patron", default="*.csv")
    args = p.parse_args()
    agrega(args.experimento, args.patron)


if __name__ == "__main__":
    main()
