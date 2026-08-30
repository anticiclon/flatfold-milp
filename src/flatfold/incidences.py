# -*- coding: utf-8 -*-


from itertools import combinations
from flatfold.geometry import se_solapan
from flatfold.exact import aplica_isometria, reflexion, compon_isometrias


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
                                  cara_plegada, aristas_caras)

    TOY = Path(__file__).resolve().parent.parent.parent / "instances" / "toy"

    for nombre in ["single_fold", "square", "map_1x3",
                   "waterbomb4", "waterbomb8"]:
        fi = read_fold(TOY / f"{nombre}.fold")
        G = crea_grafo_dirigido(fi)
        caras = caras_grafo(G)
        _, exterior = caras_definitivas(G, caras)
        idx_exterior = caras.index(exterior)
        dual = grafo_dual(G, caras, idx_exterior)
        dic_aristas = aristas_caras(G, caras)

        if len(dual) == 0:
            print(f"{nombre}: sin caras interiores conectadas")
            continue

        iso = plegar(G, dual, min(dual.nodes))
        plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
        pares = pares_comparables(plegadas)

        print("=" * 60)
        print(f"{nombre}: {len(plegadas)} caras, {len(pares)} pares")

        # --- lo que consume cruces ---
        # G, caras, plegadas, iso, dic_aristas, idx_exterior
        internas = [(a, f, g) for a, (f, g) in dic_aristas.items()
                    if f != idx_exterior and g != idx_exterior]
        print(f"  aristas interiores: {len(internas)}")
        for (u, v), f, g in internas:
            pu_f = aplica_isometria(iso[f], G.nodes[u]["coor"])
            pu_g = aplica_isometria(iso[g], G.nodes[u]["coor"])
            pv_f = aplica_isometria(iso[f], G.nodes[v]["coor"])
            pv_g = aplica_isometria(iso[g], G.nodes[v]["coor"])
            ok = (pu_f == pu_g and pv_f == pv_g)
            print(f"    arista {(u,v)} caras {f},{g}  "
                  f"desde {f}: {pu_f},{pv_f}   desde {g}: {pu_g},{pv_g}   "
                  f"{'ok' if ok else '*** DISCREPA ***'}")
            
    if nombre == "map_1x3":
        r = reflexion(G.nodes[2]["coor"], G.nodes[5]["coor"])
        print("  iso[2]:", iso[2])
        print("  iso[3]:", iso[3])
        print("  compon(r, iso[2]):", compon_isometrias(r, iso[2]))
        print("  compon(r, iso[3]):", compon_isometrias(r, iso[3]))




if __name__ == "__main__":
    main()





