"""
Patches for MACE to work in Pyodide environment.

This module stubs out missing modules and C-extension packages that are not available
in Pyodide's WASM build but are required by MACE imports.

Usage:
    from utils.patch_mace import apply_patches
    apply_patches()
"""

import sys
import types


def _matscipy_neighbour_list_compat(quantities, atoms=None, cutoff=None, positions=None, cell=None, pbc=None, **_):
    """Translate matscipy-style keyword-arg call into an ASE neighbor_list call."""
    from ase import Atoms as _Atoms
    from ase.neighborlist import neighbor_list as _ase_neighbor_list

    if atoms is None:
        atoms = _Atoms(positions=positions, cell=cell, pbc=pbc if pbc is not None else [False, False, False])
    return _ase_neighbor_list(quantities, atoms, cutoff)


class _LoggingTensorModeStub:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _capture_logs_stub(*args, **kwargs):
    return _LoggingTensorModeStub()


def patch_mace():
    """Stub modules absent in pyodide's stripped torch build and missing C-extension packages."""

    # Stub torch.testing._internal modules
    _internal = types.ModuleType("torch.testing._internal")
    _internal.__path__ = []
    _internal.__package__ = "torch.testing._internal"

    _common_utils = types.ModuleType("torch.testing._internal.common_utils")
    _common_utils.dtype_abbrs = {}

    _logging_tensor = types.ModuleType("torch.testing._internal.logging_tensor")
    _logging_tensor.LoggingTensorMode = _LoggingTensorModeStub
    _logging_tensor.capture_logs = _capture_logs_stub

    _internal.common_utils = _common_utils
    _internal.logging_tensor = _logging_tensor
    sys.modules["torch.testing._internal"] = _internal
    sys.modules["torch.testing._internal.common_utils"] = _common_utils
    sys.modules["torch.testing._internal.logging_tensor"] = _logging_tensor

    # Stub matscipy with ASE neighbor_list compatibility
    _matscipy = types.ModuleType("matscipy")
    _matscipy.__path__ = []
    _matscipy.__package__ = "matscipy"
    _matscipy_neighbours = types.ModuleType("matscipy.neighbours")
    _matscipy_neighbours.neighbour_list = _matscipy_neighbour_list_compat
    _matscipy.neighbours = _matscipy_neighbours
    sys.modules["matscipy"] = _matscipy
    sys.modules["matscipy.neighbours"] = _matscipy_neighbours

    # lmdb, orjson, h5py are C-extension packages that cannot be compiled in
    # Emscripten (no subprocess / cc support).  They are only used by the
    # training/dataset-loading paths (mace.data.lmdb_dataset,
    # mace.data.hdf5_dataset) which are pulled in by mace.data.__init__ but
    # are never exercised during Pyodide inference.  Stubs let the import
    # succeed; any actual call to these APIs would raise at runtime, which is
    # the correct behaviour for an unsupported operation.
    for _pkg in ("lmdb", "h5py"):
        if _pkg not in sys.modules:
            sys.modules[_pkg] = types.ModuleType(_pkg)

    # mace.data.atomic_data accesses torch_geometric.data.Data at class-definition
    # time.  In Pyodide the submodule attribute can be unset when mace.tools.__init__
    # is still mid-import (circular via mace.cli.visualise_train) at the point where
    # train.py runs `from . import torch_geometric` – so the subpackage exists in
    # sys.modules but .data has not yet been stamped onto it.  Pre-importing here,
    # with the matscipy stub already in place, forces the full init to complete and
    # guarantees the attribute is set before any mace.data import runs.
    try:
        import importlib as _importlib

        _tg = _importlib.import_module("mace.tools.torch_geometric")
        _tg_data = _importlib.import_module("mace.tools.torch_geometric.data")
        _tg.data = _tg_data
    except Exception:
        pass
