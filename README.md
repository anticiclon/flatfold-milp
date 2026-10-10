# flatfold-milp

A mixed-integer linear programming model for the **global flat foldability** of
*box-pleated* crease patterns. Given a pattern, the program decides whether it
admits any valid flat-folded state and, if so, returns one: the complete layer
order, face by face.

The code accompanies the paper in preparation by G. González Domínguez and
M. C. V. Nascimento (Instituto Tecnológico de Aeronáutica, Brazil).

## What it solves, and what it does not

The model covers patterns whose creases have slopes 0, infinite or ±1 on a
rational grid, and which satisfy Kawasaki's condition at every interior vertex.
Those are the modelling hypotheses; a file that fails them is rejected when it
is read or validated, not later.

Within that scope the answer is exact. All coordinates are scaled to integers,
and the three predicates that decide whether a constraint is generated —overlap,
proper crossing and collinear overlap— are comparisons between integers. A
constraint is present or absent for an exact geometric reason, never because of
a rounding decision in a degenerate configuration, which is precisely where
these constraints arise.

There is one quantity the model **cannot** optimise: the thickness of the stack
over a point of the plane. It is determined by the folded geometry rather than
by the layer order, and is the same for every valid order of the pattern.

## Installation

```bash
pip install -e .
pip install networkx gurobipy      # not declared in pyproject.toml
```

Python 3.10 or later and a Gurobi licence are required. The size-limited licence
that ships with `pip install gurobipy` is enough for the small instances
(`instances/toy/`), but not for those in `instances/flatfolder/`: they abort
with *Model too large for size-limited license*.

`matplotlib` is needed only by `experiments/genera_figuras.py`, which warns and
carries on if it is missing.

## Minimal usage

```python
from flatfold.model import (analiza_completo, construye_modelo,
                            añade_restricciones, resuelve)
from flatfold.salida import verifica_solucion, exporta_estado_fold

a = analiza_completo("instances/toy/waterbomb4.fold")

m, x = construye_modelo(a.pares)
ind = {p: i for i, p in enumerate(a.pares)}
añade_restricciones(m, x, ind, a.trip, a.cr + a.cr_oo, a.tt, a.tto, a.oo)

estado, valores = resuelve(m, x)           # estado 2 = GRB.OPTIMAL
print(verifica_solucion(a, valores))       # [] if the order is valid
exporta_estado_fold(a, valores, "estado.fold")
```

`valores` is `{(f, g): 0 or 1}` with `f < g`: it is 1 when face `f` lies above
`g`. For comparisons in the other direction there is `literal(x, u, v)`, which
returns `x[u,v]` or `1 - x[v,u]` according to the order of the indices.

The resulting `estado.fold` file opens in
[Flat-Folder](https://origamimagiro.github.io/flat-folder/), which draws the
stack and lets you check it by eye.

Two things about `analiza_completo` worth keeping in mind:

- It returns an `Analisis`, not just the constraint lists. If the pattern has
  non-convex faces it subdivides them, and then `a.fi` is no longer the
  `FoldInput` of the file —that one is `a.fi_original`— and the edge indices
  change with it. That is what `a.origen` and `a.representantes` are for: they
  translate between the two.
- The crossing incidences come in **two** lists: `a.cr` (taco-tortilla, the
  edge's two faces overlap each other) and `a.cr_oo` (tortilla-tortilla, they do
  not). They are passed to the model concatenated, as above, because they
  generate the same equality; the distinction matters when checking counts
  against an independent implementation.

## Layout

```
src/flatfold/        the model and the whole chain that generates it
tests/               401 pytest tests
scripts/             single-use diagnostic and comparison scripts
experiments/         reproducible experiment infrastructure
instances/toy/       small, hand-made cases
instances/flatfolder/  the Flat-Folder examples (see FUENTE.txt)
```

### The modules, in the order they are used

| Module | What it does |
|---|---|
| `io_fold.py` | reads the FOLD format and normalises the coordinates to the coarsest rational grid that explains them |
| `exact.py` | exact arithmetic: isometries, reflections, scaling to integers, `sigma` |
| `validate.py` | checks the hypotheses: admissible slopes, Kawasaki, sectors and local minima |
| `pattern.py` | traces the faces, builds the dual graph and folds the pattern by composing reflections |
| `geometry.py` | predicates on convex polygons: overlap, crossing, collinearity |
| `subdivision.py` | convex subdivision by vertical cuts, when some face is not convex |
| `incidences.py` | generates the five families of incidences: triples, crossings, taco-taco, taco-tortilla, tortilla-tortilla |
| `model.py` | builds the program, solves it, and exposes the variants and the derived objectives |
| `objective.py` | M/V labels, energy, prescribed assignments, Maekawa, big-little-big |
| `presolve.py` | algebraic presolve over F2 with decomposition into components |
| `salida.py` | solution verification and export to FOLD |

### Variants of the formulation

They are flags of `añade_restricciones` and of `construye_reducido`, not
separate code paths, and they work the same way with and without presolve:

- `facetas_paridad_flag` — describes each parity condition by the facets of its
  polytope instead of with an auxiliary carry variable.
- `transitividad_lazy` — does not write the triples; it separates them by
  callback. The caller must then set `m.Params.LazyConstraints = 1` and optimise
  with `resuelve_con_lazy(m)`.
- `con_maekawa`, `con_blb` — add the Maekawa and big-little-big conditions as
  redundant constraints.

Which one is preferable depends on the instance, and the answer changes with
difficulty: on instances solved without branching, the variants that strengthen
the relaxation contribute nothing, and on the hard ones they do.

### Solution verification

`salida.py` has two checkers, deliberately different from each other:

- `verifica_solucion` recomputes the incidences from the folded geometry instead
  of reusing the lists in the `Analisis`. It detects errors in the functions
  that generate them and sign errors in the translation into the model. It
  shares the predicates of `geometry.py` with the generator.
- `verifica_punto_a_punto` uses **no** function from `geometry.py` at all: its
  only predicate is a point-in-convex-polygon test by cross-product signs, in
  integer arithmetic. It is a sampling over a lattice rather than an exhaustive
  check, but it is the only genuinely independent verification path.

## Tests

```bash
python -m pytest -q
```

Those that touch instances in `instances/flatfolder/` need a full Gurobi
licence; with the size-limited one they fail because of size, not because of the
code.

## Experiments

Three pieces, kept separate on purpose: one produces results, another joins
them, the third consumes them. None decides anything the others need to know.

```bash
# one instance, one configuration, one seed -> one CSV
python -O experiments/run_experimentos.py --experimento E1 \
    --instancia instances/flatfolder/grids/004_brown_3x3_A3.fold \
    --config experiments/configs/E4_base.json \
    --salida experiments/results/E1 --semilla 0

# or a whole folder
python -O experiments/run_experimentos.py --experimento E0 \
    --todas instances/flatfolder/grids --salida experiments/results/E0

# join the individual rows
python experiments/agrega.py E1

# tables and figures of the paper, into experiments/analisis/
python experiments/genera_figuras.py
```

One file per run is what makes it possible to parallelise without the tasks
treading on each other; `agrega.py` concatenates them afterwards. Every row
carries the commit, the Gurobi version, the platform and the hash of the
`.fold`, so that a result can be placed in time.

The available experiments are `E0` (geometry and size, without solving),
`EV_stamp`, `EV_factibilidad` and `EV_recuentos` (external validation) and `E1`
(the solving experiment, parameterised by the configuration).

Configurations are JSON files in `experiments/configs/`. The keys that are read:

| Key | Effect |
|---|---|
| `caso` | `factibilidad`, `prescrita`, `parcial` or `optimizacion` |
| `presolve` | enables the algebraic presolve and solving by components |
| `presolve_gurobi` | Gurobi's `Presolve` parameter (`-1` automatic) |
| `facetas`, `lazy`, `lazy_constraints`, `blb`, `maekawa` | the variants |
| `tiempo_limite`, `hilos` | global budget, not per component |
| `k_perturbacion`, `semilla_perturbacion` | how many creases are inverted to build the repair objective |
| `porcentaje`, `semilla_subconjunto` | what fraction of creases is prescribed in the `parcial` case |

Times are measured spanning the same phases on both paths —with presolve and
without it—: `t_presolve + t_construccion` covers everything up to just before
optimising, and `t_optimize` only the solving, and the number of threads is the
same. Measuring it differently on one branch invalidates the comparison.

### On a PBS cluster

```bash
python experiments/genera_jobs.py --experimento E1 --carpeta grids \
    --configs E4_base E4_maekawa --semillas 0 1 2
qsub -v JOBS=experiments/jobs.txt experiments/tanda.pbs
```

`genera_jobs.py` writes one line of arguments per run into
`experiments/jobs.txt`, and `tanda.pbs` is an array job that runs the line
corresponding to its index. In the `.pbs` file, adjust the `-J` range, which
must match the number of lines, and the resources requested, which depend on
your cluster's queues.

## Diagnostic scripts

`scripts/` holds single-use tools, written to answer one specific question and
kept because they document how it was answered. They are not part of the
library and are not covered by the tests. The most useful ones:

- `valida_conjunto.py` — walks a set and classifies how far the chain gets:
  reading, hypotheses, topology, folding, convexity. It solves nothing.
- `resuelve_conjunto.py` — solves a folder and checks the expected feasibility.
- `compara_csv.py` — checks the constraint counts against the Flat-Folder
  reference CSVs, field by field.
- `compara_presolve.py` — measures the full model against the reduced one, and
  checks that both give the same result.
- `diagnostico.py` — reimplements Flat-Folder's `ExE_fill_BT` on our data, to
  isolate discrepancies in the collinear families.

## Instances

`instances/toy/` are small, hand-built cases, several of them deliberately
invalid (`bad_kawasaki`, `bad_slope`), which is what makes them useful as tests.

`instances/flatfolder/` are the examples from Jason S. Ku's
[flat-folder](https://github.com/origamimagiro/flat-folder) repository,
redistributed under the MIT licence; the full text is in
`LICENSE-flat-folder.txt`. **Read `instances/flatfolder/FUENTE.txt` before using
them**: it documents their provenance, what fraction of each folder falls within
the model's hypotheses, the field-by-field comparison with the reference CSVs,
and the explanation of the 24 discrepancies. Each pattern carries a
`file_author` field; when publishing results on this set it is appropriate to
credit both the source and the designers.
