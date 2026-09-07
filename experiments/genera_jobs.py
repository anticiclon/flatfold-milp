# -*- coding: utf-8 -*-
"""
Genera jobs.txt: una línea de argumentos por ejecución.

Toda la complejidad de decidir qué se ejecuta queda aquí, en un fichero
de texto legible, en lugar de en el guion de envío.
"""
import argparse
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--experimento", required=True)
    p.add_argument("--carpeta", required=True,
                   help="subcarpeta de instances/flatfolder")
    p.add_argument("--configs", nargs="+", required=True,
                   help="nombres sin .json, en experiments/configs/")
    p.add_argument("--semillas", nargs="+", type=int, default=[0])
    p.add_argument("--salida", default="experiments/jobs.txt")
    args = p.parse_args()

    carpeta = RAIZ / "instances" / "flatfolder" / args.carpeta
    instancias = sorted(carpeta.glob("*.fold"))
    if not instancias:
        raise SystemExit(f"no hay .fold en {carpeta}")

    lineas = []
    for ruta in instancias:
        for cfg in args.configs:
            for s in args.semillas:
                lineas.append(
                    f"--experimento {args.experimento} "
                    f"--instancia instances/flatfolder/{args.carpeta}/{ruta.name} "
                    f"--config experiments/configs/{cfg}.json "
                    f"--salida experiments/results/{args.experimento} "
                    f"--semilla {s}")

    destino = RAIZ / args.salida
    destino.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    print(f"{len(lineas)} ejecuciones -> {destino}")


if __name__ == "__main__":
    main()