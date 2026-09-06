# -*- coding: utf-8 -*-
"""
Productor de resultados experimentales.

Una instancia, una configuración, una semilla -> una fila.
No agrega nada ni decide nada: eso es tarea de agrega.py y
genera_figuras.py.

Uso:
    python -O experiments/run_experimentos.py --experimento E0 \
        --instancia instances/flatfolder/grids/004_brown_3x3_A3.fold \
        --salida experiments/results/E0/

    python -O experiments/run_experimentos.py --experimento E0 \
        --todas instances/flatfolder/grids \
        --salida experiments/results/E0/
"""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path
import re
from flatfold.model import (analiza_completo, construye_modelo,
                                añade_restricciones)
from flatfold.presolve import Infactible

RAIZ = Path(__file__).resolve().parent.parent


# ===== metadatos de entorno ===============================================

def commit_actual():
    """El commit de HEAD, o cadena vacía si no hay repositorio."""
    try:
        r = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                           cwd=RAIZ, capture_output=True, text=True,
                           timeout=10)
        return r.stdout.strip() if r.returncode == 0 else ""
    except Exception:
        return ""


def version_gurobi():
    try:
        import gurobipy as gb
        return ".".join(str(v) for v in gb.gurobi.version())
    except Exception:
        return ""


def hash_fichero(ruta):
    return hashlib.sha256(Path(ruta).read_bytes()).hexdigest()[:12]


def metadatos(ruta):
    return {
        "instancia": Path(ruta).name,
        "carpeta": Path(ruta).parent.name,
        "hash_instancia": hash_fichero(ruta),
        "commit": commit_actual(),
        "gurobi": version_gurobi(),
        "python": platform.python_version(),
        "optimizado": int(sys.flags.optimize > 0),
        "fecha": time.strftime("%Y-%m-%d %H:%M:%S"),
    }


# ===== clasificación de rechazos ==========================================

def clasifica(mensaje):
    """Etiqueta un ValueError de la cadena según su texto."""
    if "retícula" in mensaje:
        return "reticula"
    if "pendiente" in mensaje:
        return "pendiente"
    if "Kawasaki" in mensaje or "grado impar" in mensaje:
        return "kawasaki"
    if "no convexa" in mensaje:
        return "no_convexa"
    if "se rompe" in mensaje or "inconsistencia" in mensaje:
        return "plegado"
    return "otro"


# ===== E0: caracterización del conjunto ===================================

def experimento_E0(ruta, config, semilla):
    """
    Recorre la cadena hasta las incidencias y registra qué se obtiene.
    No resuelve ningún modelo.
    """

    fila = metadatos(ruta)
    fila.update(estado="", error="", t_geometria=0.0,
                vertices=0, aristas=0, caras=0, subdividido=0,
                pares=0, tripletas=0, cruces=0,
                taco_taco=0, taco_tortilla=0, tortilla_tortilla=0)

    t = time.perf_counter()
    try:
        a = analiza_completo(ruta)
    except ValueError as e:
        fila["estado"] = clasifica(str(e))
        fila["error"] = str(e)[:200]
        fila["t_geometria"] = round(time.perf_counter() - t, 3)
        return fila
    except (Infactible, Exception) as e:
        fila["estado"] = "excepcion"
        fila["error"] = f"{type(e).__name__}: {e}"[:200]
        fila["t_geometria"] = round(time.perf_counter() - t, 3)
        return fila

    fila["t_geometria"] = round(time.perf_counter() - t, 3)
    fila["estado"] = "ok" if a.pares else "sin_pliegues"
    fila.update(
        vertices=len(a.fi.vertices), aristas=len(a.fi.aristas),
        caras=len(a.plegadas), subdividido=int(a.subdividido),
        pares=len(a.pares), tripletas=len(a.trip), cruces=len(a.cr),
        taco_taco=len(a.tt), taco_tortilla=len(a.tto),
        tortilla_tortilla=len(a.oo))
    return fila


EXPERIMENTOS = {"E0": experimento_E0}


# ===== escritura ==========================================================

def escribe(fila, carpeta, nombre):
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    destino = carpeta / f"{nombre}.csv"
    with open(destino, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fila))
        w.writeheader()
        w.writerow(fila)
    return destino


def nombre_salida(experimento, ruta, config, semilla):
    partes = [experimento, Path(ruta).parent.name, Path(ruta).stem]
    if config:
        partes.append(Path(config).stem)
    if semilla is not None:
        partes.append(f"s{semilla}")
    return "__".join(partes)



# ===== EV: validación =====================================================

# valores conocidos de la secuencia de stamp folding
STAMP = {2: 2, 3: 6, 4: 16, 5: 50, 6: 144, 7: 462}


def experimento_EV_stamp(ruta, config, semilla):
    """
    Número de plegados planos de la tira de n celdas, contra los valores
    conocidos de la literatura. Es el contraste externo más fuerte del
    proyecto: no depende de ninguna otra implementación.

    Solo se aplica a las instancias map_1xN; el resto se registra como
    'no_aplica' para que la fila exista y el recuento cuadre.
    """

    fila = metadatos(ruta)
    fila.update(n=0, esperados=0, obtenidos=0, coincide=0,
                estado="", error="", t_total=0.0)

    m_nombre = re.fullmatch(r"map_1x(\d+)", Path(ruta).stem)
    if not m_nombre:
        fila["estado"] = "no_aplica"
        return fila

    n = int(m_nombre.group(1))
    fila["n"] = n
    if n not in STAMP:
        fila["estado"] = "sin_referencia"
        return fila
    fila["esperados"] = STAMP[n]

    limite = config.get("limite_pool", 5000)
    t = time.perf_counter()
    try:
        a = analiza_completo(ruta)
        m, x = construye_modelo(a.pares)
        ind = {p: i for i, p in enumerate(a.pares)}
        añade_restricciones(m, x, ind, a.trip, a.cr, a.tt, a.tto, a.oo)
        m.Params.OutputFlag = 0
        m.Params.PoolSearchMode = 2
        m.Params.PoolSolutions = limite
        if semilla is not None:
            m.Params.Seed = semilla
        m.optimize()

        fila["obtenidos"] = m.SolCount
        if m.SolCount >= limite:
            fila["estado"] = "pool_lleno"
        else:
            fila["estado"] = "ok"
            fila["coincide"] = int(m.SolCount == STAMP[n])
    except Exception as e:
        fila["estado"] = "excepcion"
        fila["error"] = f"{type(e).__name__}: {e}"[:200]

    fila["t_total"] = round(time.perf_counter() - t, 3)
    return fila


def experimento_EV_factibilidad(ruta, config, semilla):
    """
    Las instancias admitidas deben resultar factibles: son patrones que
    Flat-Folder pliega. Una infactible indicaría que al modelo le sobra
    alguna restricción.

    No se fija la asignación M/V: se comprueba que existe algún orden de
    capas, no que la etiqueta del fichero sea realizable.
    """
    import gurobipy as gb
    from flatfold.model import (analiza_completo, construye_modelo,
                                añade_restricciones)

    fila = metadatos(ruta)
    fila.update(estado="", error="", pares=0, caras=0,
                t_geometria=0.0,
                t_construccion=0.0,
                t_optimize=0.0,
                t_solver=0.0,
                gurobi_status=0)


    limite = config.get("tiempo_limite", 300)

    t = time.perf_counter()
    try:
        a = analiza_completo(ruta)
    except ValueError as e:
        fila["estado"] = clasifica(str(e))
        fila["error"] = str(e)[:200]
        fila["t_geometria"] = round(time.perf_counter() - t, 3)
        return fila
    except Exception as e:
        fila["estado"] = "excepcion"
        fila["error"] = f"{type(e).__name__}: {e}"[:200]
        fila["t_geometria"] = round(time.perf_counter() - t, 3)
        return fila

    fila["t_geometria"] = round(time.perf_counter() - t, 3)
    fila["pares"] = len(a.pares)
    fila["caras"] = len(a.plegadas)

    if not a.pares:
        fila["estado"] = "sin_pliegues"
        return fila
    
    t = time.perf_counter()
    try:
        m, x = construye_modelo(a.pares)
        ind = {p: i for i, p in enumerate(a.pares)}
        añade_restricciones(m, x, ind, a.trip, a.cr, a.tt, a.tto, a.oo)
    
        m.Params.OutputFlag = 0
        m.Params.IntegralityFocus = 1
        m.Params.SolutionLimit = 1
        m.Params.TimeLimit = limite
    
        if semilla is not None:
            m.Params.Seed = semilla
    
        fila["t_construccion"] = round(time.perf_counter() - t, 3)
    
        t2 = time.perf_counter()
        m.optimize()
        fila["t_optimize"] = round(time.perf_counter() - t2, 3)
    
        fila["gurobi_status"] = m.Status
    
        if m.SolCount > 0:
            fila["estado"] = "factible"
        elif m.Status == gb.GRB.INFEASIBLE:
            fila["estado"] = "INFACTIBLE"
        else:
            fila["estado"] = "sin_resolver"
    
    except Exception as e:
        fila["estado"] = "excepcion"
        fila["error"] = f"{type(e).__name__}: {e}"[:200]
        
    fila["t_solver"] = round(
    fila["t_construccion"] + fila["t_optimize"], 3
    )
    return fila


_REFERENCIAS = {}


def carga_referencia(carpeta):
    """{number: fila} del CSV de referencia de esa carpeta, cacheado."""
    if carpeta not in _REFERENCIAS:
        ruta = RAIZ / "instances" / "flatfolder" / f"{carpeta}_data.csv"
        if not ruta.exists():
            _REFERENCIAS[carpeta] = None
        else:
            with open(ruta, encoding="utf-8") as f:
                _REFERENCIAS[carpeta] = {r["number"]: r
                                         for r in csv.DictReader(f)}
    return _REFERENCIAS[carpeta]


# nuestro campo -> campo del CSV de referencia
CAMPOS_REF = ["vertices", "edges", "faces", "variables", "taco-taco",
              "taco-tortilla", "tortilla-tortilla", "transitivity"]


def experimento_EV_recuentos(ruta, config, semilla):
    """
    Contrasta los recuentos de variables y restricciones contra la salida
    de referencia de Flat-Folder.

    AVISO: el CSV de referencia puede corresponder a una versión anterior
    de los ficheros. Por eso se comparan también 'vertices' y 'edges',
    que no deberían cambiar: si difieren, la discrepancia es de datos y
    no de modelo.
    """
    from flatfold.model import analiza_completo

    fila = metadatos(ruta)
    fila.update(estado="", error="", n_difieren=0, difieren="")
    for c in CAMPOS_REF:
        fila[f"mio_{c}"] = 0
        fila[f"ref_{c}"] = 0

    carpeta = Path(ruta).parent.name
    ref = carga_referencia(carpeta)
    if ref is None:
        fila["estado"] = "sin_referencia"
        return fila

    numero = Path(ruta).stem.split("_")[0]
    if numero not in ref:
        fila["estado"] = "sin_fila"
        return fila
    suya = ref[numero]

    try:
        a = analiza_completo(ruta)
    except ValueError as e:
        fila["estado"] = clasifica(str(e))
        fila["error"] = str(e)[:200]
        return fila
    except Exception as e:
        fila["estado"] = "excepcion"
        fila["error"] = f"{type(e).__name__}: {e}"[:200]
        return fila

    mio = {
        "vertices": len(a.fi_original.vertices),
        "edges": len(a.fi_original.aristas),
        "faces": len(a.plegadas) if a.plegadas else len(a.caras) - 1,
        "variables": len(a.pares),
        "taco-taco": len(a.tt),
        "taco-tortilla": len(a.cr) + len(a.tto),
        "tortilla-tortilla": len(a.cr_oo) + len(a.oo),
        "transitivity": len(a.trip),
    }

    difieren = []
    for c in CAMPOS_REF:
        fila[f"mio_{c}"] = mio[c]
        fila[f"ref_{c}"] = int(suya[c])
        if mio[c] != int(suya[c]):
            difieren.append(c)

    fila["n_difieren"] = len(difieren)
    fila["difieren"] = ";".join(difieren)
    fila["estado"] = "coincide" if not difieren else "difiere"
    return fila


EXPERIMENTOS = {"E0": experimento_E0,
                "EV_stamp": experimento_EV_stamp,
                "EV_factibilidad": experimento_EV_factibilidad,
                "EV_recuentos": experimento_EV_recuentos}



# ===== principal ==========================================================

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--experimento", required=True,
                   choices=sorted(EXPERIMENTOS))
    p.add_argument("--instancia", help="ruta de un fichero .fold")
    p.add_argument("--todas", help="carpeta: recorre todos sus .fold")
    p.add_argument("--config", help="fichero JSON de configuración")
    p.add_argument("--semilla", type=int, default=None)
    p.add_argument("--salida", required=True, help="carpeta de resultados")
    args = p.parse_args()

    if bool(args.instancia) == bool(args.todas):
        p.error("indica --instancia o --todas, pero no ambos")

    config = {}
    if args.config:
        config = json.loads(Path(args.config).read_text(encoding="utf-8"))

    funcion = EXPERIMENTOS[args.experimento]
    rutas = ([Path(args.instancia)] if args.instancia
             else sorted(Path(args.todas).glob("*.fold")))
    if not rutas:
        p.error(f"no hay ficheros .fold en {args.todas}")

    for ruta in rutas:
        fila = funcion(ruta, config, args.semilla)
        nombre = nombre_salida(args.experimento, ruta, args.config,
                               args.semilla)
        destino = escribe(fila, args.salida, nombre)
        print(f"{ruta.name}: {fila['estado']}  -> {destino.name}")


if __name__ == "__main__":
    main()