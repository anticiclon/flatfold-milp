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
import random
from flatfold.model import (analiza_completo, construye_modelo,
                                añade_restricciones, 
                                resuelve_con_lazy, 
                                resuelve_por_componentes)
from flatfold.presolve import Infactible
import gurobipy as gb
from flatfold.objective import (etiquetas_mv, energia, fija_etiquetas,
                                    restricciones_blb, restricciones_maekawa,
                                    etiquetas_por_pliegue_original)
from flatfold.presolve import (literales_a_ecuaciones, reduce_completo,
                                   componentes, reparte_por_componentes,
                                   energia_reducida, maekawa_reducido)


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


def subconjunto_prescrito(prescritas, porcentaje, semilla_sub):
    """
    Elige al azar un porcentaje de los pliegues para fijar.

    La semilla es INDEPENDIENTE de la de Gurobi: así, al variar esta
    última se mide solo la varianza del solver y no un problema distinto.
    """
    if porcentaje >= 100:
        return dict(prescritas)
    aristas = sorted(prescritas)
    k = round(len(aristas) * porcentaje / 100)
    elegidas = random.Random(semilla_sub).sample(aristas, k)
    return {e: prescritas[e] for e in elegidas}




def asignacion_perturbada(prescritas, k, semilla_pert):
    """
    Invierte k pliegues elegidos al azar sobre la asignación del fichero.

    La del fichero es realizable —son patrones que Flat-Folder pliega—,
    así que perturbarla produce una asignación que con alta probabilidad
    no lo es. El óptimo de la reparación queda entonces entre 1 y k:
    devolver esos k pliegues a su sitio es siempre una solución posible,
    y a veces bastan menos.

    Esa cota superior es una comprobación gratis de corrección: un óptimo
    mayor que k indicaría un error.

    La semilla es INDEPENDIENTE de la de Gurobi, para que al variar esta
    última se mida la varianza del solver y no un problema distinto.
    """
    aristas = sorted(prescritas)
    k = min(k, len(aristas))
    invertidas = set(random.Random(semilla_pert).sample(aristas, k))
    return {e: (1 - v if e in invertidas else v)
            for e, v in prescritas.items()}





def resuelve_con_config(a, config, semilla):
    """
    Construye y resuelve el modelo según las banderas de la configuración.

    Devuelve un diccionario con el resultado y los tiempos desglosados.
    Los dos caminos —con y sin presolve— quedan detrás de la misma
    interfaz, de modo que el experimento no tiene que distinguirlos.

    Los tiempos se miden abarcando las mismas fases en ambos caminos:
    t_presolve más t_construccion cubren hasta justo antes de optimizar,
    y t_optimize solo la resolución. Medirlos de otro modo en una rama
    invalidaría cualquier comparación.

    Casos: factibilidad, prescrita, parcial, optimizacion.
    """
    r = {"vars": len(a.pares), "vars_reducidas": 0, "componentes": 0,
         "t_presolve": 0.0, "t_construccion": 0.0, "t_optimize": 0.0,
         "estado": "", "valor": "", "nodos": 0, "gurobi_status": 0,
         "n_fijadas": 0, "k_perturbacion": ""}

    modo = config.get("caso", "factibilidad")
    limite = config.get("tiempo_limite", 300)

    etiq = etiquetas_mv(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext)
    orig = etiquetas_por_pliegue_original(etiq, a.representantes)
    blb = (restricciones_blb(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext,
                             etiq) if config.get("blb") else None)
    mk = (restricciones_maekawa(a.G, a.fi, a.iso, a.dic_aristas, a.idx_ext,
                                etiq) if config.get("maekawa") else None)

    # qué etiquetas se fijan, según el caso
    prescritas = None
    if modo in ("prescrita", "parcial"):
        prescritas = subconjunto_prescrito(
            a.fi_original.dicc_montana_valle,
            100 if modo == "prescrita" else config.get("porcentaje", 50),
            config.get("semilla_subconjunto", 0))
        r["n_fijadas"] = len(prescritas)

    # el objetivo de reparación: la asignación del fichero con k pliegues
    # invertidos. Invertirlos todos no sirve: suele seguir siendo
    # realizable y el óptimo sale 0, sin dificultad para el solver.
    objetivo = None
    if modo == "optimizacion":
        k = config.get("k_perturbacion", 3)
        r["k_perturbacion"] = k
        objetivo = asignacion_perturbada(
            a.fi_original.dicc_montana_valle, k,
            config.get("semilla_perturbacion", 0))

    modo_solver = "optimizacion" if modo == "optimizacion" else "factibilidad"

    # ---------------- camino con presolve ----------------
    if config.get("presolve"):
        t = time.perf_counter()
        try:
            ig, par, fij = literales_a_ecuaciones(
                a.pares, a.cr + a.cr_oo, a.tto, a.oo, a.tt,
                etiquetas=orig if prescritas else None,
                prescritas=prescritas, blb=blb)
            mapeo, restantes, vivas, _ = reduce_completo(
                len(a.pares), ig, par, fij, a.trip, a.pares)
            mkr = maekawa_reducido(mapeo, a.pares, mk) if mk else None
            comps = componentes(mapeo, restantes, a.trip, a.pares, mkr)

            terminos, const = ({}, 0)
            if objetivo is not None:
                terminos, const = energia_reducida(mapeo, a.pares, orig,
                                                   objetivo)

            datos = reparte_por_componentes(comps, vivas, restantes,
                                            terminos, mkr)
        except Infactible:
            # el presolve detecta la contradicción sin llegar al solver
            r["t_presolve"] = round(time.perf_counter() - t, 3)
            r["estado"] = "infactible"
            return r

        r["t_presolve"] = round(time.perf_counter() - t, 3)
        r["vars_reducidas"] = sum(len(d["raices"]) for d in datos)
        r["componentes"] = len(comps)

        t = time.perf_counter()
        valor, _, nodos = resuelve_por_componentes(
            datos, const, modo=modo_solver, tiempo=limite,
            presolve_gurobi=config.get("presolve_gurobi", -1),
            facetas=config.get("facetas", False),
            lazy=config.get("lazy", False),
            semilla=semilla)
        r["t_optimize"] = round(time.perf_counter() - t, 3)
        r["nodos"] = int(nodos)

        if valor is None or valor is False:
            r["estado"] = "sin_resolver"
        else:
            r["estado"] = "ok"
            if modo_solver == "optimizacion":
                r["valor"] = int(round(valor))
        return r

    # ---------------- camino sin presolve ----------------
    t = time.perf_counter()
    m, x = construye_modelo(a.pares)
    ind = {p: i for i, p in enumerate(a.pares)}
    añade_restricciones(m, x, ind, a.trip, a.cr + a.cr_oo, a.tt, a.tto,
                        a.oo, blb=blb, maekawa=mk,
                        con_blb=bool(blb), con_maekawa=bool(mk),
                        facetas_paridad_flag=config.get("facetas", False),
                        transitividad_lazy=config.get("lazy", False))
    m.Params.OutputFlag = 0
    m.Params.IntegralityFocus = 1
    m.Params.MIPGap = 0
    m.Params.Threads = 1
    m.Params.TimeLimit = limite
    m.Params.Presolve = config.get("presolve_gurobi", -1)
    if semilla is not None:
        m.Params.Seed = semilla

    if prescritas:
        fija_etiquetas(m, x, orig, prescritas)

    if modo == "optimizacion":
        m.setObjective(energia(x, orig, objetivo), gb.GRB.MINIMIZE)
    else:
        m.Params.SolutionLimit = 1

    r["vars_reducidas"] = len(a.pares)
    r["componentes"] = 1
    r["t_construccion"] = round(time.perf_counter() - t, 3)

    t = time.perf_counter()
    if config.get("lazy"):
        resuelve_con_lazy(m)
    else:
        m.optimize()
    r["t_optimize"] = round(time.perf_counter() - t, 3)

    r["nodos"] = int(m.NodeCount)
    r["gurobi_status"] = m.Status
    if m.SolCount > 0:
        r["estado"] = "ok"
        if modo == "optimizacion":
            r["valor"] = int(round(m.ObjVal))
    elif m.Status == gb.GRB.INFEASIBLE:
        r["estado"] = "infactible"
    else:
        r["estado"] = "sin_resolver"
    return r





def experimento_E1(ruta, config, semilla):
    """
    Presolve sobre F2, cruzado con el presolve de Gurobi.

    Las cuatro celdas del 2x2 se definen en los ficheros de
    configuración; aquí solo se ejecuta lo que digan. El mismo código
    sirve para E2, E3 y E4 cambiando las banderas, de modo que esas
    variantes no necesitan función propia.
    """

    fila = metadatos(ruta)
    fila.update(config)
    fila.update(semilla=semilla if semilla is not None else "",
                estado="", error="", t_geometria=0.0, caras=0)

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
    fila["caras"] = len(a.plegadas)

    if not a.pares:
        fila["estado"] = "sin_pliegues"
        return fila

    try:
        fila.update(resuelve_con_config(a, config, semilla))
    except Exception as e:
        fila["estado"] = "excepcion"
        fila["error"] = f"{type(e).__name__}: {e}"[:200]

    return fila


EXPERIMENTOS = {"E0": experimento_E0,
                "EV_stamp": experimento_EV_stamp,
                "EV_factibilidad": experimento_EV_factibilidad,
                "EV_recuentos": experimento_EV_recuentos,
                "E1": experimento_E1}


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