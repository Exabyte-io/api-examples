# mat3ra.notebooks_utils — Developer Guide

## Architecture

Four progressive layers, each building on the previous:

```
primitive/  →  core/  →  ipython/  →  pyodide/
```

Top-level files (`auth.py`, `io.py`, `ui.py`, `plot.py`, `settings.py`, `material.py`) are thin
routing/re-export adapters that keep notebook code environment-agnostic.

---

## Layers

### `primitive/`
Only Python stdlib. No third-party packages, no domain knowledge.
Enums, environment detection, logger, CLI prompt helpers.

### `core/`
Third-party packages allowed (`mat3ra.api_client`, `requests`, `numpy`, …).
No IPython, no ipywidgets, no browser/pyodide APIs.

```
core/api/              — Platform credentials and OIDC auth.
core/io.py             — Plain-Python file IO and HTTP.
core/entity/<domain>/  — One folder per domain (material, workflow, job, compute, property).
    api.py             — REST API calls only.
    io.py              — Local filesystem reads/writes only.
    analysis.py        — Pure computation (numpy/pymatgen/ase). No API, no display.
    job.py             — Cross-entity helpers (e.g. property derived from job data).
```

### `ipython/`
Imports `IPython.display`, `ipywidgets`, or generates HTML/JS for notebook output.
Must work in JupyterLab, Colab, and VS Code notebooks — not browser/pyodide specific.

```
ipython/ui.py              — Generic cell output: display_JSON, image grid, viewer HTML/JS.
ipython/io.py              — Browser file-download helpers.
ipython/plot/              — Domain-agnostic plot primitives (_plotly.py, _matplotlib.py).
ipython/entity/<domain>/   — Domain-specific display.
    visualize.py           — Renders domain objects into notebook cells.
    plot.py                — Domain-specific charts (RDF, strain, EOS, …).
```

### `pyodide/`
Calls `micropip`, `pyodide.http.pyfetch`, `BroadcastChannel`, or `js` (Emscripten/WASM APIs).
JupyterLite-specific overrides of `core/` or `ipython/` capabilities.

```
pyodide/io.py          — JS data bridge: kernel↔host shared state, file writes.
pyodide/ui.py          — Async input widgets (PyodideFuture / BroadcastChannel).
pyodide/runtime.py     — Interruptible loop and abort controller.
pyodide/packages/      — micropip installation and config.yml parsing.
pyodide/api/           — Host-to-kernel token injection.
```

---

## Dependency Rules

```
primitive/  →  (nothing from this package)
core/       →  primitive/
ipython/    →  primitive/, core/
pyodide/    →  primitive/, core/, ipython/
top-level   →  any layer (routing only, no new logic)
```

Within `core/entity/`: domains must not import each other; `api.py`, `io.py`, and `analysis.py`
must stay separate (no cross-imports within the same entity folder).

Within `ipython/entity/`: may import from `core/entity/<same domain>/` and `ipython/` primitives,
but not from other entity domains.

---

## Where does my code go?

| Question | Destination |
|---|---|
| Uses only stdlib? | `primitive/` |
| Calls a 3rd-party package, produces no output? | `core/` |
| Reads/writes domain objects from disk? | `core/entity/<domain>/io.py` |
| Calls the REST API? | `core/entity/<domain>/api.py` |
| Computes/transforms domain data (no display)? | `core/entity/<domain>/analysis.py` |
| Renders into a notebook cell (IPython/HTML)? | `ipython/` |
| Domain-specific chart? | `ipython/entity/<domain>/plot.py` |
| Domain-specific visualiser? | `ipython/entity/<domain>/visualize.py` |
| Calls micropip / pyfetch / BroadcastChannel? | `pyodide/` |
| Branches on environment to pick implementation? | top-level routing file |
