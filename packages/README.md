# Packages for Pyodide/JupyterLite

This directory contains pre-built Python wheels used by the JupyterLite notebooks. Wheels are loaded at runtime via `emfs:/drive/packages/` in [config.yml](../config.yml).

## Wheel Sources

### Built by [Exabyte-io/build-pyodide](https://github.com/Exabyte-io/build-pyodide)

The following wheels are compiled for the Emscripten/WASM target (`cp311-emscripten_3_1_4x_wasm32`) using the CI workflows in the [build-pyodide](https://github.com/Exabyte-io/build-pyodide) repository:

| Wheel | Notes |
|-------|-------|
| `torch-2.1.0a0-cp311-cp311-emscripten_3_1_45_wasm32.whl` | PyTorch for Pyodide (WASM build) |
| `pymatgen-2024.*.whl` (emscripten) | Pymatgen with C extensions compiled for WASM |
| `spglib-2.5.0-cp311-cp311-emscripten_3_1_46_wasm32.whl` | Spglib with C extensions compiled for WASM |
| `icet-3.0-cp311-cp311-emscripten_3_1_45_wasm32.whl` | ICET for WASM |
| `numba-0.60.0-cp311-cp311-emscripten_3_1_45_wasm32.whl` | Numba for WASM |
| `pydantic_core-2.18.2-py3-none-any.whl` | Pydantic core (Rust → pure-Python stub) |

### Downloaded from PyPI (unmodified)

These are standard pure-Python wheels downloaded directly from [PyPI](https://pypi.org) with no modifications:

| Wheel | PyPI |
|-------|------|
| `antlr4_python3_runtime-4.9.3-py3-none-any.whl` | [antlr4-python3-runtime](https://pypi.org/project/antlr4-python3-runtime/4.9.3/) |
| `paginate-0.5.6-py3-none-any.whl` | [paginate](https://pypi.org/project/paginate/0.5.6/) |
| `pydantic-2.7.1-py3-none-any.whl` | [pydantic](https://pypi.org/project/pydantic/2.7.1/) |
| `pymatgen-2024.4.13-py3-none-any.whl` | [pymatgen](https://pypi.org/project/pymatgen/2024.4.13/) (pure-Python variant) |
| `ruamel.yaml-0.17.32-py3-none-any.whl` | [ruamel.yaml](https://pypi.org/project/ruamel.yaml/0.17.32/) |
| `spglib-2.0.2-py3-none-any.whl` | [spglib](https://pypi.org/project/spglib/2.0.2/) (pure-Python variant) |
| `watchdog-2.3.1-py3-none-any.whl` | [watchdog](https://pypi.org/project/watchdog/2.3.1/) |

### Custom-built wheels (not in build-pyodide)

These wheels required manual modifications to work in Pyodide and are **not** built by the build-pyodide CI. Instructions to reproduce each are below.

---

#### `mattersim-1.1.2-py3-none-any.whl` (70 KB)

**Source**: [microsoft/mattersim](https://github.com/microsoft/mattersim) v1.1.2

**Why custom**: The upstream MatterSim package uses Cython-compiled extensions (`_m3gnet_cython.pyx`) that cannot run in Pyodide's WASM runtime. The custom wheel replaces these with pure NumPy equivalents.

**How to reproduce**:

```bash
# 1. Download source
pip download mattersim==1.1.2 --no-deps --no-binary :all:
tar xf mattersim-1.1.2.tar.gz && cd mattersim-1.1.2

# 2. Remove Cython build requirement
#    Edit pyproject.toml: remove "cython" from build-system.requires

# 3. Replace Cython extensions with NumPy implementations
#    In mattersim/forcefield/m3gnet/modules/:
#    - Remove _m3gnet_cython.pyx and any .so/.c compiled files
#    - The M3GNet module code falls back to pure NumPy when Cython is unavailable

# 4. Build pure-Python wheel
pip wheel . --no-deps --no-build-isolation
# Output: mattersim-1.1.2-py3-none-any.whl
```

**Runtime patches**: Requires `apply_all_patches(include_mattersim=True)` in [torch.py](../src/py/mat3ra/notebooks_utils/pyodide/packages/torch.py) which stubs `loguru`, `azure`, `e3nn` JIT, `torch_geometric`, and `torch_runstats`.

---

#### `sevenn-0.12.1-py3-none-any.whl` (9.2 MB)

**Source**: [MDIL-SNU/SevenNet](https://github.com/MDIL-SNU/SevenNet) v0.12.1 (PyPI)

**Why custom**: The upstream wheel is 42.6 MB and includes 4 pretrained models + a CUDA-only native library (`pair_d3.so`). We strip it down to the single model needed for inference.

**How to reproduce**:

```bash
# 1. Download the upstream wheel
pip download sevenn==0.12.1 --no-deps
# sevenn-0.12.1-py3-none-any.whl (~42.6 MB)

# 2. Extract, strip, and repack
python3 << 'EOF'
import zipfile, os, shutil

src_whl = "sevenn-0.12.1-py3-none-any.whl"
work_dir = "sevenn_work"

# Extract
with zipfile.ZipFile(src_whl, 'r') as z:
    z.extractall(work_dir)

pkg = os.path.join(work_dir, "sevenn")

# Remove CUDA D3 dispersion library (not needed for CPU inference)
for f in ["pair_d3.so"]:
    p = os.path.join(pkg, f)
    if os.path.exists(p):
        os.remove(p)

# Remove CUDA D3 source code
shutil.rmtree(os.path.join(pkg, "pair_e3gnn"), ignore_errors=True)

# Keep only 7net-0 (latest), remove other pretrained models
pp = os.path.join(pkg, "pretrained_potentials")
for name in os.listdir(pp):
    if name != "SevenNet_0__11Jul2024":
        shutil.rmtree(os.path.join(pp, name))

# Repack
out_whl = "sevenn-0.12.1-py3-none-any.whl"
with zipfile.ZipFile(out_whl, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(work_dir):
        for f in files:
            full = os.path.join(root, f)
            arcname = os.path.relpath(full, work_dir)
            zf.write(full, arcname)

shutil.rmtree(work_dir)
EOF
# Output: sevenn-0.12.1-py3-none-any.whl (~9.2 MB)
```

**What was removed** (saves ~33 MB):
- `pair_d3.so` — CUDA-only D3 dispersion native library (1.8 MB)
- `pair_e3gnn/` — C source for D3 dispersion (3.2 MB)
- `pretrained_potentials/SevenNet_0__22May2024/` — older 7net-0 checkpoint (9.8 MB)
- `pretrained_potentials/SevenNet_MF_0/` — multi-fidelity model (9.8 MB)
- `pretrained_potentials/SevenNet_l3i5/` — large model variant (14 MB)

**What was kept**:
- `pretrained_potentials/SevenNet_0__11Jul2024/checkpoint_sevennet_0.pth` — 7net-0 model (9.8 MB)
- All Python source code (nn/, train/, calculator.py, etc.)

**Runtime patches**: Requires `apply_all_patches(include_sevennet=True)` in [torch.py](../src/py/mat3ra/notebooks_utils/pyodide/packages/torch.py) which stubs `pandas`, `braceexpand`, `tqdm`, provides a `torch_geometric.data.Data` replacement, and sets `e3nn` to eager mode (no JIT).

---

## Models

The `models/` subdirectory contains pretrained model checkpoint files:

| File | Model | Size |
|------|-------|------|
| `mattersim-v1.0.0-1M.pth` | MatterSim M3GNet (1M params) | ~4 MB |

> **Note**: The SevenNet 7net-0 model is bundled inside the `sevenn` wheel itself under `pretrained_potentials/`.
