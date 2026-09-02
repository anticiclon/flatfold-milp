# -*- coding: utf-8 -*-
"""
Reimplementa ExE_fill_BT de Flat-Folder sobre nuestros datos y compara su
salida con la de colineales(). Diagnóstico de una instancia concreta.
"""
from pathlib import Path
from collections import defaultdict, Counter

from flatfold.io_fold import read_fold, Kind
from flatfold.pattern import (crea_grafo_dirigido, caras_grafo,
                              caras_definitivas, grafo_dual, plegar,
                              cara_plegada, aristas_caras)
from flatfold.geometry import (recta_soporte, solapan_colinealmente,
                               se_solapan, caja, interseca_cajas,
                               poligono_de_caja, doble_area_poligono)
from flatfold.exact import aplica_isometria
from flatfold.incidences import pares_comparables, colineales

RAIZ = Path(__file__).resolve().parent.parent
EJEMPLOS = RAIZ / "instances" / "flatfolder"

# --- configuración ---
CARPETA = "grids"
FICHERO = "121_brown_4x4_A99.fold"
ARISTA = (5, 10)          # arista concreta a examinar; None para omitir
# ---------------------

RUTA = EJEMPLOS / CARPETA / FICHERO


def exe_al_estilo_ku(aristas_por_recta, conjunto_pares):
    """Traducción literal de ExE_fill_BT (conversion.js)."""
    def solapan(a, b):
        return 1 if (min(a, b), max(a, b)) in conjunto_pares else 0

    tt, tto, oo, motivos = [], [], [], Counter()

    for key, aristas in aristas_por_recta.items():
        for i in range(len(aristas)):
            for j in range(i + 1, len(aristas)):
                s1, _, f1, f2, a1 = aristas[i]
                s2, _, f3, f4, a2 = aristas[j]
                if solapan_colinealmente(s1, s2) is None:
                    continue

                f1f2 = solapan(f1, f2)
                f1f3 = solapan(f1, f3)
                f1f4 = solapan(f1, f4)
                choice = (f1f2 << 2) | (f1f3 << 1) | f1f4
                motivos[choice] += 1

                if choice == 7:
                    tt.append((f1, f2, f3, f4, key, a1, a2))
                elif choice == 6:
                    tto.append((f1, f2, f3, key, a1, a2))
                elif choice == 5:
                    tto.append((f1, f2, f4, key, a1, a2))
                elif choice == 3:
                    tto.append((f3, f4, f1, key, a1, a2))
                elif choice == 0:
                    tto.append((f3, f4, f2, key, a1, a2))
                elif choice == 1:
                    oo.append((f1, f2, f4, f3, key, a1, a2))
                elif choice == 2:
                    oo.append((f1, f2, f3, f4, key, a1, a2))
                # choice == 4: nada

    return tt, tto, oo, motivos


def examina_arista(G, dic_aristas, plegadas, arista):
    """Detalle geométrico de una arista concreta."""
    u, v = arista
    clave = (u, v) if u < v else (v, u)
    if clave not in dic_aristas:
        print(f"\n  la arista {clave} no está en dic_aristas (¿es frontera?)")
        return
    f, g = dic_aristas[clave]
    print(f"\n  arista {clave}: tipo {G[u][v]['tipo'].name}, caras {f} y {g}")
    print(f"    sin plegar: {G.nodes[u]['coor']} -> {G.nodes[v]['coor']}")
    print(f"    cara {f} plegada: {plegadas[f]}")
    print(f"    cara {g} plegada: {plegadas[g]}")
    poli = poligono_de_caja(interseca_cajas([caja(plegadas[f]),
                                             caja(plegadas[g])]))
    area = doble_area_poligono(poli) if len(poli) >= 3 else 0
    print(f"    doble área de la intersección: {area}")
    print(f"    se_solapan: {se_solapan(plegadas[f], plegadas[g])}")


def main():
    print(FICHERO)

    fi = read_fold(RUTA)
    print("  tipos de arista:", Counter(k.name for k in fi.tipos))

    G = crea_grafo_dirigido(fi)
    caras = caras_grafo(G)
    _, exterior = caras_definitivas(G, caras)
    idx_ext = caras.index(exterior)
    dual = grafo_dual(G, caras, idx_ext)
    dic_aristas = aristas_caras(G, caras)
    iso = plegar(G, dual, min(dual.nodes))
    plegadas = {f: cara_plegada(G, caras[f], iso[f]) for f in iso}
    pares = set(pares_comparables(plegadas))

    # --- borde geométrico frente al del fichero ---
    geom = {frozenset(a) for a, (f, g) in dic_aristas.items()
            if idx_ext in (f, g)}
    fich = {frozenset(a) for a, k in zip(fi.aristas, fi.tipos)
            if k == Kind.BOUNDARY}
    print(f"  aristas: {len(dic_aristas)} totales, "
          f"{len(dic_aristas) - len(geom)} interiores, {len(geom)} de borde")
    print(f"  borde geométrico vs fichero: "
          f"{'coinciden' if geom == fich else 'DIFIEREN'}")

    # --- grupos por recta soporte (una sola vez) ---
    grupos = defaultdict(list)
    for (u, v), (f, g) in dic_aristas.items():
        if f == idx_ext or g == idx_ext:
            continue
        seg = (aplica_isometria(iso[f], G.nodes[u]["coor"]),
               aplica_isometria(iso[f], G.nodes[v]["coor"]))
        grupos[recta_soporte(seg)].append((seg, G[u][v]["tipo"], f, g, (u, v)))

    entradas = [t for aristas in grupos.values() for t in aristas]
    print(f"  en grupos: {len(entradas)} entradas, "
          f"{len({t[4] for t in entradas})} aristas distintas, "
          f"{len(grupos)} rectas")
    print("  tipos en los grupos:",
          Counter(t[1].name for t in entradas))

    # --- comparación de criterios ---
    ku_tt, ku_tto, ku_oo, motivos = exe_al_estilo_ku(grupos, pares)
    mi_tt, mi_tto, mi_oo = colineales(G, caras, plegadas, iso,
                                      dic_aristas, idx_ext)

    print(f"\n  Ku:  taco-taco {len(ku_tt)}, taco-tortilla {len(ku_tto)}, "
          f"tortilla-tortilla {len(ku_oo)}")
    print(f"  mío: taco-taco {len(mi_tt)}, taco-tortilla {len(mi_tto)}, "
          f"tortilla-tortilla {len(mi_oo)}")
    print("  reparto por choice (f1f2 f1f3 f1f4):")
    for c in sorted(motivos):
        print(f"    {c:03b} ({c}): {motivos[c]}")

    mias = {frozenset(sum(cuatro, ())) for cuatro in mi_tt}
    suyas = {frozenset((f1, f2, f3, f4)) for f1, f2, f3, f4, *_ in ku_tt}
    print(f"  cuádruplas solo en Ku: {len(suyas - mias)}")
    print(f"  cuádruplas solo mías:  {len(mias - suyas)}")

    if ARISTA:
        examina_arista(G, dic_aristas, plegadas, ARISTA)


if __name__ == "__main__":
    main()