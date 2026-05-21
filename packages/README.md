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

#### `chgnet-0.3.8-py3-none-any.whl` (4.3 MB)

**Source**: [CederGroupHub/chgnet](https://github.com/CederGroupHub/chgnet) v0.3.8 (PyPI)

**Why custom**: The upstream wheel includes a Cython-compiled graph converter extension (`cygraph.so`) that cannot run in Pyodide's WASM runtime. CHGNet has a built-in `algorithm="legacy"` fallback that uses pure Python when the C extension is unavailable, so we simply strip the `.so` file.

**How to reproduce**:

```bash
# 1. Download the upstream wheel
pip download chgnet==0.3.8 --no-deps
# chgnet-0.3.8-*.whl (~4.5 MB)

# 2. Extract, strip the Cython extension, and repack
python3 << 'EOF'
import zipfile, os, shutil

src_whl = next(f for f in os.listdir('.') if f.startswith('chgnet-0.3.8'))
work_dir = "chgnet_work"

# Extract
with zipfile.ZipFile(src_whl, 'r') as z:
    z.extractall(work_dir)

# Remove the Cython .so extension (runtime falls back to pure Python)
pkg = os.path.join(work_dir, "chgnet", "graph")
for f in os.listdir(pkg):
    if f.startswith("cygraph") and f.endswith(".so"):
        os.remove(os.path.join(pkg, f))
        print(f"Removed: {f}")

# Repack
out_whl = "chgnet-0.3.8-py3-none-any.whl"
with zipfile.ZipFile(out_whl, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(work_dir):
        for f in files:
            full = os.path.join(root, f)
            arcname = os.path.relpath(full, work_dir)
            zf.write(full, arcname)

shutil.rmtree(work_dir)
EOF
# Output: chgnet-0.3.8-py3-none-any.whl (~4.3 MB)
```

**What was removed**:
- `chgnet/graph/cygraph.*.so` — Cython-compiled graph converter (runtime falls back to `algorithm="legacy"`)

**What was kept**:
- All Python source code
- Pretrained model v0.3.0 checkpoint (~4 MB, bundled in `chgnet/pretrained/`)

**Runtime patches**: Requires `apply_all_patches(include_chgnet=True)` in [torch.py](../src/py/mat3ra/notebooks_utils/pyodide/packages/torch.py) which stubs `nvidia_smi` (GPU memory detection), `Cython` (build-time dep), and `palettable` (color palette lib imported via `pymatgen.util.plotting`).

---

#### `nequip-0.15.0-py3-none-any.whl` (254 KB)

**Source**: [mir-group/nequip](https://github.com/mir-group/nequip) v0.15.0 (PyPI)

**Why custom**: The upstream NequIP wheel declares heavy training dependencies (`hydra-core`, `lightning`, `torchmetrics`, `lmdb`) that are not needed for inference and would fail to install in Pyodide. The custom wheel strips these from `Requires-Dist` metadata.

**How to reproduce**:

```bash
# 1. Create a venv with build tools
python3 -m venv /tmp/nequip_build && source /tmp/nequip_build/bin/activate
pip install build wheel setuptools

# 2. Download and extract source
pip download nequip==0.15.0 --no-deps --no-binary :all:
tar xf nequip-0.15.0.tar.gz && cd nequip-0.15.0

# 3. Strip heavy dependencies from pyproject.toml
#    Remove these lines from [project].dependencies:
#      "hydra-core", "lightning", "torchmetrics>=1.6.0",
#      "lmdb", "tqdm", "requests"
#    Keep: torch, numpy, matscipy, ase, e3nn, pyyaml

# 4. Build the wheel
python -m build --wheel
# Output: dist/nequip-0.15.0-py3-none-any.whl (~254 KB)
```

**Runtime patches**: Requires `apply_all_patches(include_nequip=True)` in [torch.py](../src/py/mat3ra/notebooks_utils/pyodide/packages/torch.py) which stubs `hydra` (with working `instantiate`), `lightning`, `pytorch_lightning`, `torchmetrics`, `lmdb`, `matscipy` (uses ASE neighbor lists via `NEQUIP_NL=ase`), and patches `e3nn` for eager mode (no JIT).

---

## Models

The `models/` subdirectory contains pretrained model checkpoint files:

| File | Model | Size | Notes |
|------|-------|------|-------|
| `mattersim-v1.0.0-1M.pth` | MatterSim M3GNet (1M params) | ~4 MB | Direct checkpoint |
| `nequip-oam-s-config-sd.pth` | NequIP-OAM-S (618K params) | ~2.4 MB | Config + state_dict extracted from `.nequip.zip` package |

> **Note**: The SevenNet 7net-0 model is bundled inside the `sevenn` wheel under `pretrained_potentials/`. The CHGNet v0.3.0 model is bundled inside the `chgnet` wheel under `chgnet/pretrained/`.

### NequIP-OAM-S model extraction

The `nequip-oam-s-config-sd.pth` file contains the model architecture config and state_dict extracted from the NequIP package format (`.nequip.zip`). This is necessary because NequIP's standard loading paths use `torch.package.PackageImporter` or `torch.jit.load`, neither of which work in Pyodide.

**How to reproduce**:

```bash
# 1. Set up venv with NequIP
python3 -m venv /tmp/nequip_venv && source /tmp/nequip_venv/bin/activate
pip install nequip==0.15.0

# 2. Download the NequIP-OAM-S model from Zenodo
wget https://zenodo.org/records/18775904/files/NequIP-OAM-S-0.1.nequip.zip

# 3. Extract config + state_dict
python3 << 'PYEOF'
import torch, io, yaml

pkg_path = "NequIP-OAM-S-0.1.nequip.zip"
imp = torch.package.PackageImporter(pkg_path)

# Patch for CPU-only loading
orig = torch.storage._load_from_bytes
def _load_cpu(b):
    return torch.load(io.BytesIO(b), map_location="cpu", weights_only=False)
torch.storage._load_from_bytes = _load_cpu
pkg_model = imp.load_pickle(package="model", resource="eager_model.pkl", map_location="cpu")
torch.storage._load_from_bytes = orig

# Extract state_dict and metadata
inner_sd = pkg_model.sole_model.state_dict()
metadata = yaml.safe_load(imp.load_text(package="model", resource="package_metadata.txt"))
type_names = [metadata["atom_types"][i] for i in range(len(metadata["atom_types"]))]

# Extract per-type energy scales and shifts
scales_t = pkg_model.sole_model.model.func.per_type_energy_scale_shift.scales.data
shifts_t = pkg_model.sole_model.model.func.per_type_energy_scale_shift.shifts.data
scales_dict = {tn: scales_t[i, 0].item() for i, tn in enumerate(type_names)}
shifts_dict = {tn: shifts_t[i, 0].item() for i, tn in enumerate(type_names)}

# Extract avg_num_neighbors from scatter_norm_factor
snf = pkg_model.sole_model.model.func.layer0_convnet.conv.scatter_norm_factor
avg_nn = 1.0 / (snf ** 2)

# Save config + state_dict
torch.save({
    "type_names": type_names,
    "r_max": 4.5,
    "irreps_edge_sh": "1x0e+1x1e",
    "type_embed_num_features": 32,
    "feature_irreps_hidden": ["128x0e+64x1e", "128x0e"],
    "radial_mlp_depth": [1, 1],
    "radial_mlp_width": [128, 128],
    "avg_num_neighbors": avg_nn,
    "per_type_energy_scales": scales_dict,
    "per_type_energy_shifts": shifts_dict,
    "has_zbl": True,
    "state_dict": inner_sd,
    "metadata": metadata,
}, "nequip-oam-s-config-sd.pth")
PYEOF
# Output: nequip-oam-s-config-sd.pth (~2.4 MB)
```

The model is then loaded at runtime using `load_nequip_model()` from [torch.py](../src/py/mat3ra/notebooks_utils/pyodide/packages/torch.py), which rebuilds the architecture using `FullNequIPGNNModel` + manually adds the ZBL pair potential, then loads the state_dict.
