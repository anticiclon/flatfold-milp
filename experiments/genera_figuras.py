# -*- coding: utf-8 -*-
"""
Consumidor de resultados: lee los CSV agregados y produce las cifras y
las figuras del artículo.

No resuelve nada ni toca Gurobi. Las tablas salen como CSV con las cifras
ya calculadas, para maquetarlas a mano; las figuras, en PDF vectorial.

Uso:
    python experiments/genera_figuras.py
"""
import csv
import statistics
from collections import defaultdict
from pathlib import Path

RAIZ = Path(__file__).resolve().parent.parent
RESULTADOS = RAIZ / "experiments" / "results"
SALIDA = RAIZ / "experiments" / "analisis"


# ===== utilidades ==========================================================

def lee(nombre):
    ruta = RESULTADOS / f"{nombre}.csv"
    if not ruta.exists():
        print(f"  (falta {ruta.name}, se salta)")
        return []
    with open(ruta, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def escribe(nombre, filas, campos=None):
    SALIDA.mkdir(parents=True, exist_ok=True)
    if not filas:
        return
    campos = campos or list(filas[0])
    ruta = SALIDA / f"{nombre}.csv"
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos, restval="")
        w.writeheader()
        w.writerows(filas)
    print(f"  {len(filas)} filas -> {ruta.name}")


def etiqueta(r):
    """Nombre corto de la configuración, a partir de sus banderas."""
    partes = [x for x, c in (("facetas", "facetas"), ("lazy", "lazy"),
                             ("lazyparam", "lazy_constraints"),
                             ("blb", "blb"), ("maekawa", "maekawa"))
              if r.get(c) == "True"]
    return "+".join(partes) if partes else "base"


def tiempo(r):
    return (float(r["t_presolve"]) + float(r["t_construccion"])
            + float(r["t_optimize"]))


def resumen(ts):
    """Estadísticos de una lista de tiempos."""
    return {"n": len(ts),
            "total": round(sum(ts), 1),
            "media": round(statistics.mean(ts), 3),
            "mediana": round(statistics.median(ts), 3),
            "max": round(max(ts), 1)}


# ===== instancias comunes ==================================================

def comunes(filas, clave):
    """
    Instancias que TODAS las configuraciones resuelven, para una carpeta.

    Sin esto las medias no son comparables: cada configuración promediaría
    sobre un conjunto distinto, y la que menos resuelve saldría
    artificialmente rápida porque las difíciles no cuentan.
    """
    por_cfg = defaultdict(set)
    todas = set()
    for r in filas:
        por_cfg[clave(r)].add(r["instancia"])
        todas.add(r["instancia"])
    resueltas = defaultdict(set)
    for r in filas:
        if r["estado"] == "ok":
            resueltas[clave(r)].add(r["instancia"])
    if not resueltas:
        return set()
    return set.intersection(*resueltas.values())


# ===== tabla E1: el 2x2 del presolve =======================================

def tabla_E1(filas):
    clave = lambda r: (r["presolve"], r["presolve_gurobi"])
    salida = []

    for carpeta in sorted({r["carpeta"] for r in filas}):
        sub = [r for r in filas
               if r["carpeta"] == carpeta and etiqueta(r) == "base"]
        if not sub:
            continue
        comun = comunes(sub, clave)

        for cfg in sorted({clave(r) for r in sub}):
            de_cfg = [r for r in sub if clave(r) == cfg]
            ok = [r for r in de_cfg if r["estado"] == "ok"]
            en_comun = [tiempo(r) for r in ok if r["instancia"] in comun]
            if not en_comun:
                continue
            fila = {"carpeta": carpeta,
                    "presolve": cfg[0],
                    "presolve_gurobi": cfg[1],
                    "resueltas": len({r["instancia"] for r in ok}),
                    "comunes": len(comun)}
            fila.update(resumen(en_comun))
            vars_antes = [int(r["vars"]) for r in ok if r["vars"]]
            vars_desp = [int(r["vars_reducidas"]) for r in ok
                         if r["vars_reducidas"]]
            fila["vars"] = sum(vars_antes)
            fila["vars_reducidas"] = sum(vars_desp)
            fila["reduccion_pct"] = (
                round(100 * (1 - sum(vars_desp) / sum(vars_antes)), 1)
                if sum(vars_antes) else "")
            salida.append(fila)
    return salida


# ===== tabla de variantes: E2, E3, E4 ======================================

def tabla_variantes(filas):
    # solo la configuración base del presolve: presolve activo, Gurobi auto
    filas = [r for r in filas
             if r["presolve"] == "True" and r["presolve_gurobi"] == "-1"]
    clave = etiqueta
    salida = []

    for carpeta in sorted({r["carpeta"] for r in filas}):
        sub = [r for r in filas if r["carpeta"] == carpeta]
        comun = comunes(sub, clave)

        base = None
        for cfg in sorted({clave(r) for r in sub}):
            de_cfg = [r for r in sub if clave(r) == cfg]
            ok = [r for r in de_cfg if r["estado"] == "ok"]
            en_comun = [tiempo(r) for r in ok if r["instancia"] in comun]
            if not en_comun:
                continue
            fila = {"carpeta": carpeta, "variante": cfg,
                    "resueltas": len({r["instancia"] for r in ok}),
                    "comunes": len(comun)}
            fila.update(resumen(en_comun))
            nodos = [int(float(r["nodos"])) for r in ok if r["nodos"]]
            fila["nodos_total"] = sum(nodos)
            if cfg == "base":
                base = fila["total"]
            salida.append(fila)

        # razón frente a la base, una vez conocida
        for fila in salida:
            if fila["carpeta"] == carpeta and base:
                fila["vs_base"] = round(fila["total"] / base, 2)
    return salida


# ===== tablas de validación y del conjunto =================================

def tabla_validacion():
    salida = []
    for r in lee("EV_stamp"):
        if r["estado"] == "ok":
            salida.append({"prueba": "stamp_folding", "caso": f"n={r['n']}",
                           "esperado": r["esperados"],
                           "obtenido": r["obtenidos"],
                           "coincide": r["coincide"]})

    rec = lee("EV_recuentos")
    if rec:
        from collections import Counter
        c = Counter(r["estado"] for r in rec)
        for estado, n in sorted(c.items()):
            salida.append({"prueba": "recuentos", "caso": estado,
                           "esperado": "", "obtenido": n, "coincide": ""})

    fac = lee("EV_factibilidad")
    if fac:
        from collections import Counter
        c = Counter((r["carpeta"], r["estado"]) for r in fac)
        for (carpeta, estado), n in sorted(c.items()):
            salida.append({"prueba": "factibilidad",
                           "caso": f"{carpeta}/{estado}",
                           "esperado": "", "obtenido": n, "coincide": ""})
    return salida


def tabla_conjunto():
    filas = lee("E0")
    if not filas:
        return []
    from collections import Counter
    c = Counter((r["carpeta"], r["estado"]) for r in filas)
    salida = [{"carpeta": k[0], "estado": k[1], "instancias": n}
              for k, n in sorted(c.items())]

    for carpeta in sorted({r["carpeta"] for r in filas}):
        caras = sorted(int(r["caras"]) for r in filas
                       if r["carpeta"] == carpeta and r["estado"] == "ok")
        if caras:
            salida.append({"carpeta": carpeta, "estado": "_caras_min",
                           "instancias": caras[0]})
            salida.append({"carpeta": carpeta, "estado": "_caras_mediana",
                           "instancias": statistics.median(caras)})
            salida.append({"carpeta": carpeta, "estado": "_caras_max",
                           "instancias": caras[-1]})
    return salida


# ===== figuras =============================================================

def figuras(filas):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib no disponible, se saltan las figuras)")
        return

    SALIDA.mkdir(parents=True, exist_ok=True)

    # --- perfil de rendimiento de las variantes ---
    base = [r for r in filas
            if r["presolve"] == "True" and r["presolve_gurobi"] == "-1"
            and r["carpeta"] == "grids"]
    por_cfg = defaultdict(dict)
    for r in base:
        if r["estado"] == "ok":
            k = (r["instancia"], r["semilla"])
            por_cfg[etiqueta(r)][k] = tiempo(r)
    if por_cfg:
        comun = set.intersection(*(set(d) for d in por_cfg.values()))
        mejor = {k: min(d[k] for d in por_cfg.values()) for k in comun}

        fig, ax = plt.subplots(figsize=(6, 4))

        colores = [
            "tab:blue",
            "tab:orange",
            "tab:green",
            "tab:red",
            "tab:purple",
            "tab:brown",
            "tab:pink"
        ]

        for i, (cfg, d) in enumerate(sorted(por_cfg.items())):
            razones = sorted(d[k] / mejor[k] if mejor[k] > 0 else 1
                             for k in comun)
            ys = [(j + 1) / len(razones) for j in range(len(razones))]
            ax.step(
                razones,
                ys,
                where="post",
                label=cfg,
                linestyle="-",
                color=colores[i % len(colores)],
                linewidth=1.5
            )
        ax.set_xscale("log")
        # ax.set_xlabel("razón respecto del mejor tiempo")
        ax.set_xlabel("performance ratio")
        # ax.set_ylabel("fracción de instancias")
        ax.set_ylabel("fraction of instances")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(SALIDA / "perfil_variantes.pdf")
        plt.close(fig)
        print(f"  perfil_variantes.pdf ({len(comun)} puntos)")

    # --- tiempo frente a número de caras ---
    pts = [(int(r["caras"]), tiempo(r)) for r in filas
           if r["estado"] == "ok" and etiqueta(r) == "base"
           and r["presolve"] == "True" and r["presolve_gurobi"] == "-1"
           and r["caras"]]
    if pts:
        fig, ax = plt.subplots(figsize=(6.2, 4))
        ax.scatter([p[0] for p in pts], [p[1] for p in pts],
                   s=9, color="#1b6ca8", alpha=0.35, edgecolors="none")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("interior faces")
        ax.set_ylabel("time (s)")
        ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(SALIDA / "tiempo_vs_caras.pdf")
        plt.close(fig)
        print(f"  tiempo_vs_caras.pdf ({len(pts)} puntos)")



# ===== principal ===========================================================

def main():
    e1 = lee("E1")
    print(f"E1: {len(e1)} filas")

    print("\nTablas:")
    escribe("tabla_E1", tabla_E1(e1))
    escribe("tabla_variantes", tabla_variantes(e1))
    escribe("tabla_validacion", tabla_validacion())
    escribe("tabla_conjunto", tabla_conjunto())

    print("\nFiguras:")
    figuras(e1)


if __name__ == "__main__":
    main()