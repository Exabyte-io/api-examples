from collections import namedtuple

import numpy as np
import torch

# Define return types to mimic PyTorch's named tuples
EigRet = namedtuple("linalg_eig", ["eigenvalues", "eigenvectors"])  # type: ignore
EighRet = namedtuple("linalg_eigh", ["eigenvalues", "eigenvectors"])  # type: ignore


def _to_np(tensor):
    return tensor.detach().cpu().numpy()


def _to_torch(array, device, dtype=None):
    if dtype:
        return torch.tensor(array, device=device, dtype=dtype)
    return torch.tensor(array, device=device)


# --- Define the Patches ---


def patch_solve(A, B, *args, **kwargs):
    return _to_torch(np.linalg.solve(_to_np(A), _to_np(B)), A.device)


def patch_inv(A, *args, **kwargs):
    return _to_torch(np.linalg.inv(_to_np(A)), A.device)


def patch_det(A, *args, **kwargs):
    return _to_torch(np.linalg.det(_to_np(A)), A.device)


def patch_cholesky(A, *args, **kwargs):
    # NumPy's cholesky operates on the lower triangle by default, same as PyTorch
    return _to_torch(np.linalg.cholesky(_to_np(A)), A.device)


def patch_eig(A, *args, **kwargs):
    vals, vecs = np.linalg.eig(_to_np(A))
    return EigRet(_to_torch(vals, A.device), _to_torch(vecs, A.device))


def patch_eigh(A, UPLO="L", *args, **kwargs):
    vals, vecs = np.linalg.eigh(_to_np(A), UPLO=UPLO)
    return EighRet(_to_torch(vals, A.device), _to_torch(vecs, A.device))


# --- Apply the Patches to PyTorch ---

torch.linalg.solve = patch_solve
torch.linalg.inv = patch_inv
torch.inverse = patch_inv  # Alias
torch.linalg.det = patch_det
torch.det = patch_det  # Alias
torch.linalg.cholesky = patch_cholesky
torch.linalg.eig = patch_eig
torch.linalg.eigh = patch_eigh


# Fix numpy
def _tensor_array_compat(self, dtype=None):
    """Replacement for Tensor.__array__ in Pyodide where tensor.numpy() is unavailable."""
    arr = np.array(self.tolist())
    return arr.astype(dtype) if dtype is not None else arr


torch.Tensor.__array__ = _tensor_array_compat
torch.Tensor.numpy = lambda self: np.array(self.detach().tolist())


# Keep the SciPy LU patches we made earlier just in case
import scipy.linalg  # noqa: E402

LUFactorReturn = namedtuple("LUFactorReturn", ["LU", "pivots"])


def patch_lu_factor(A, *args, **kwargs):
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


def patch_lu_solve(LU, pivots, B, *args, **kwargs):
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


torch.linalg.lu_factor = patch_lu_factor
torch.linalg.lu_solve = patch_lu_solve

print("All major torch.linalg functions successfully patched to use NumPy/SciPy!")
