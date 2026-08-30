# -*- coding: utf-8 -*-


from itertools import combinations
from flatfold.geometry import se_solapan


def pares_comparables(caras_plegadas):
    """
    Input:  {cara: polígono plegado}
    Output: lista de pares (f, g) con f < g cuyas imágenes se solapan
            con área positiva.

    Cada par es una variable binaria del modelo: x_fg = 1 significa que la
    cara f está por encima de la g. Los pares que no se solapan no generan
    variable porque nunca compiten por el mismo punto del plano.
    """
    pares = []
    # Todas las posibles combinaciones de 2 caras
    for f, g in combinations(sorted(caras_plegadas), 2):
        if se_solapan(caras_plegadas[f], caras_plegadas[g]):
            pares.append((f, g))
    return pares





def main():
    from pathlib import Path
    from flatfold.io_fold import read_fold
    from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                                  caras_definitivas, grafo_dual, plegar,
                                  cara_plegada)

    TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

    for nombre in ["single_fold", "square", "map_1x3",
                   "waterbomb4", "waterbomb8"]:
        fi = read_fold(TOY / f"{nombre}.fold")
        G = crea_grafo_dirigido(fi)
        caras = caras_grafo(G)
        _, exterior = caras_definitivas(G, caras)
        dual = grafo_dual(G, caras, caras.index(exterior))

        if len(dual) == 0:                 # patrón sin pliegues
            print(f"{nombre}: sin caras interiores conectadas")
            continue

        iso = plegar(G, dual, min(dual.nodes))
        plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
        
        
        print(nombre, "caras:", len(caras), "interiores:", len(iso))
        for f, p in sorted(plegadas.items()):
            print("  ", f, p)


        pares = pares_comparables(plegadas)

        print(f"{nombre}: {len(plegadas)} caras, {len(pares)} pares -> {pares}")


if __name__ == "__main__":
    main()





