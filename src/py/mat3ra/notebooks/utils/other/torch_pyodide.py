"""
Patches for PyTorch and related packages to work in Pyodide environment.

This module provides patches for various torch-related packages that don't work
in Pyodide's WASM environment, organized by functionality.

Usage:
    from mat3ra.notebooks.utils.other.torch_pyodide import (
        patch_torch_linalg,
        patch_torch_testing,
        patch_matscipy,
        patch_lmdb_h5py,
        patch_mace_tools,
    )

    # Apply all patches
    patch_torch_linalg()
    patch_torch_testing()
    patch_matscipy()
    patch_lmdb_h5py()
    patch_mace_tools()
"""

import sys
import types
from collections import namedtuple

import numpy as np
import torch

# Define return types to mimic PyTorch's named tuples
EigRet = namedtuple("linalg_eig", ["eigenvalues", "eigenvectors"])  # type: ignore
EighRet = namedtuple("linalg_eigh", ["eigenvalues", "eigenvectors"])  # type: ignore
LUFactorReturn = namedtuple("LUFactorReturn", ["LU", "pivots"])


def _to_np(tensor):
    return tensor.detach().cpu().numpy()


def _to_torch(array, device, dtype=None):
    if dtype:
        return torch.tensor(array, device=device, dtype=dtype)
    return torch.tensor(array, device=device)


# ==============================================================================
# Torch linalg patches
# ==============================================================================


def _patch_solve(A, B, *args, **kwargs):
    return _to_torch(np.linalg.solve(_to_np(A), _to_np(B)), A.device)


def _patch_inv(A, *args, **kwargs):
    return _to_torch(np.linalg.inv(_to_np(A)), A.device)


def _patch_det(A, *args, **kwargs):
    return _to_torch(np.linalg.det(_to_np(A)), A.device)


def _patch_cholesky(A, *args, **kwargs):
    return _to_torch(np.linalg.cholesky(_to_np(A)), A.device)


def _patch_eig(A, *args, **kwargs):
    vals, vecs = np.linalg.eig(_to_np(A))
    return EigRet(_to_torch(vals, A.device), _to_torch(vecs, A.device))


def _patch_eigh(A, UPLO="L", *args, **kwargs):
    vals, vecs = np.linalg.eigh(_to_np(A), UPLO=UPLO)
    return EighRet(_to_torch(vals, A.device), _to_torch(vecs, A.device))


def _patch_lu_factor(A, *args, **kwargs):
    import scipy.linalg

    A_np = _to_np(A)
    if A_np.ndim > 2:
        orig_shape = A_np.shape
        A_reshaped = A_np.reshape(-1, orig_shape[-2], orig_shape[-1])
        lu_list, piv_list = [], []
        for mat in A_reshaped:
            lu, piv = scipy.linalg.lu_factor(mat)
            lu_list.append(lu)
            piv_list.append(piv)
        LU_np = np.stack(lu_list).reshape(orig_shape)
        piv_np = np.stack(piv_list).reshape(orig_shape[:-2] + (-1,))
    else:
        LU_np, piv_np = scipy.linalg.lu_factor(A_np)
    return LUFactorReturn(_to_torch(LU_np, A.device, A.dtype), _to_torch(piv_np, A.device, torch.int32))


def _patch_lu_solve(LU, pivots, B, *args, **kwargs):
    import scipy.linalg

    LU_np, piv_np, B_np = _to_np(LU), _to_np(pivots), _to_np(B)
    if LU_np.ndim > 2:
        orig_shape_B = B_np.shape
        LU_reshaped = LU_np.reshape(-1, LU_np.shape[-2], LU_np.shape[-1])
        piv_reshaped = piv_np.reshape(-1, piv_np.shape[-1])
        is_vector = B_np.ndim == LU_np.ndim - 1
        B_reshaped = (
            B_np.reshape(-1, B_np.shape[-1], 1) if is_vector else B_np.reshape(-1, B_np.shape[-2], B_np.shape[-1])
        )
        X_list = [scipy.linalg.lu_solve((lu, piv), b) for lu, piv, b in zip(LU_reshaped, piv_reshaped, B_reshaped)]
        X_np = np.stack(X_list).reshape(orig_shape_B)
    else:
        X_np = scipy.linalg.lu_solve((LU_np, piv_np), B_np)
    return _to_torch(X_np, B.device, B.dtype)


def _tensor_array_compat(self, dtype=None):
    """Replacement for Tensor.__array__ in Pyodide where tensor.numpy() is unavailable."""
    arr = np.array(self.tolist())
    return arr.astype(dtype) if dtype is not None else arr


def patch_torch_linalg():
    """
    Patch torch.linalg functions to use NumPy/SciPy implementations.

    This fixes LAPACK functions that are not available in Pyodide's WASM torch build.
    """
    torch.linalg.solve = _patch_solve
    torch.linalg.inv = _patch_inv
    torch.inverse = _patch_inv  # Alias
    torch.linalg.det = _patch_det
    torch.det = _patch_det  # Alias
    torch.linalg.cholesky = _patch_cholesky
    torch.linalg.eig = _patch_eig
    torch.linalg.eigh = _patch_eigh
    torch.linalg.lu_factor = _patch_lu_factor
    torch.linalg.lu_solve = _patch_lu_solve

    # Fix numpy compatibility
    torch.Tensor.__array__ = _tensor_array_compat
    torch.Tensor.numpy = lambda self: np.array(self.detach().tolist())

    # Fix torch.compiler.is_compiling for Pyodide
    if not hasattr(torch, "compiler"):
        torch.compiler = types.ModuleType("torch.compiler")
    if not hasattr(torch.compiler, "is_compiling"):
        torch.compiler.is_compiling = lambda: False

    print("✓ Torch linalg patches applied")


# ==============================================================================
# Torch testing patches
# ==============================================================================


class _LoggingTensorModeStub:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


def _capture_logs_stub(*args, **kwargs):
    return _LoggingTensorModeStub()


def patch_torch_testing():
    """
    Stub torch.testing._internal modules that are missing in Pyodide.

    These modules are used by MACE imports but not needed for inference.
    """
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

    print("✓ Torch testing patches applied")


# ==============================================================================
# Matscipy patches
# ==============================================================================


def _matscipy_neighbour_list_compat(quantities, atoms=None, cutoff=None, positions=None, cell=None, pbc=None, **_):
    """Translate matscipy-style keyword-arg call into an ASE neighbor_list call."""
    from ase import Atoms as _Atoms
    from ase.neighborlist import neighbor_list as _ase_neighbor_list

    if atoms is None:
        atoms = _Atoms(positions=positions, cell=cell, pbc=pbc if pbc is not None else [False, False, False])
    return _ase_neighbor_list(quantities, atoms, cutoff)


def patch_matscipy():
    """
    Stub matscipy package with ASE neighbor_list compatibility.

    Matscipy is a C-extension package that cannot be compiled for WASM.
    MACE uses it for neighbor calculations, which we redirect to ASE.
    """
    _matscipy = types.ModuleType("matscipy")
    _matscipy.__path__ = []
    _matscipy.__package__ = "matscipy"
    _matscipy_neighbours = types.ModuleType("matscipy.neighbours")
    _matscipy_neighbours.neighbour_list = _matscipy_neighbour_list_compat
    _matscipy.neighbours = _matscipy_neighbours
    sys.modules["matscipy"] = _matscipy
    sys.modules["matscipy.neighbours"] = _matscipy_neighbours

    print("✓ Matscipy patches applied")


# ==============================================================================
# LMDB and HDF5 patches
# ==============================================================================


def patch_mace_training():
    """
    Stub lmdb and h5py packages.

    These are C-extension packages used by MACE's training/dataset code
    but not needed for inference. Stubs allow imports to succeed.
    """
    for _pkg in ("lmdb", "h5py"):
        if _pkg not in sys.modules:
            sys.modules[_pkg] = types.ModuleType(_pkg)

    print("✓ LMDB and HDF5 stubs applied")


# ==============================================================================
# MACE tools patches
# ==============================================================================


def patch_mace_tools():
    """
    Fix MACE's torch_geometric import order issues in Pyodide.

    In Pyodide, torch_geometric.data may not be set during circular imports.
    Pre-importing ensures the attribute is available when MACE needs it.
    """
    try:
        import importlib as _importlib

        _tg = _importlib.import_module("mace.tools.torch_geometric")
        _tg_data = _importlib.import_module("mace.tools.torch_geometric.data")
        _tg.data = _tg_data
        print("✓ MACE tools patches applied")
    except Exception as e:
        print(f"⚠ MACE tools patches skipped: {e}")


# ==============================================================================
# Convenience function to apply all patches
# ==============================================================================


def apply_all_patches():
    """Apply all torch and MACE patches for Pyodide in one call."""
    patch_torch_linalg()
    patch_torch_testing()
    patch_matscipy()
    patch_mace_training()
    patch_mace_tools()
    print("\n✅ All Pyodide patches applied successfully!")
