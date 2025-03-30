import numpy as np
import scipy.sparse as sp
import numba
import utils.consts
from typing import Optional
from sklearn.utils import sparsefuncs

# A series of support functions modified from the scanpy API,
# to allow functionality with expression data stored
# with smaller precision than int64 or float64

# Support functions from preprocessing/_utils.py that allow sparse functionality
def sparse_mean_var_major_axis(data, indices, indptr, major_len, minor_len, dtype):
    """
    Computes mean and variance for a sparse array for the major axis.
    Given arrays for a csr matrix, returns the means and variances for each
    row back.
    """
    means = np.zeros(major_len, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)
    for i in range(major_len):
        startptr = indptr[i]
        endptr = indptr[i + 1]
        counts = endptr - startptr
        for j in range(startptr, endptr):
            means[i] += data[j]
        means[i] /= minor_len
        for j in range(startptr, endptr):
            diff = data[j] - means[i]
            variances[i] += diff * diff
        variances[i] += (minor_len - counts) * means[i] ** 2
        variances[i] /= minor_len
    return means, variances


def sparse_mean_var_minor_axis(data, indices, major_len, minor_len, dtype):
    """
    Computes mean and variance for a sparse matrix for the minor axis.
    Given arrays for a csr matrix, returns the means and variances for each
    column back.
    """
    non_zero = indices.shape[0]
    means = np.zeros(minor_len, dtype=dtype)
    variances = np.zeros_like(means, dtype=dtype)
    counts = np.zeros(minor_len, dtype=np.int64)
    for i in range(non_zero):
        col_ind = indices[i]
        means[col_ind] += data[i]
    for i in range(minor_len):
        means[i] /= major_len
    for i in range(non_zero):
        col_ind = indices[i]
        diff = data[i] - means[col_ind]
        variances[col_ind] += diff * diff
        counts[col_ind] += 1
    for i in range(minor_len):
        variances[i] += (major_len - counts[i]) * means[i] ** 2
        variances[i] /= major_len
    return means, variances


def sparse_mean_variance_axis(mtx: sp.spmatrix,
        axis: int,
        running_mean_var_dtype: str = 'float64'):
    """
    This code and internal functions are based on sklearns
    `sparsefuncs.mean_variance_axis`.
    Modifications:
    * allow deciding on the output type, which can increase accuracy when calculating the mean and variance of 32bit floats.
    * This doesn't currently implement support for null values, but could.
    * Uses numba not cython
    """
    assert axis in (0, 1)
    if isinstance(mtx, sp.csr_matrix):
        ax_minor = 1
        shape = mtx.shape
    elif isinstance(mtx, sp.csc_matrix):
        ax_minor = 0
        shape = mtx.shape[::-1]
    else:
        raise ValueError("This function only works on sparse csr and csc matrices")
    if axis == ax_minor:
        return sparse_mean_var_major_axis(
            mtx.data, mtx.indices, mtx.indptr, *shape, running_mean_var_dtype)
    else:
        return sparse_mean_var_minor_axis(mtx.data, mtx.indices, *shape, running_mean_var_dtype)


# Copy of preprocessing/_utils.py, _get_mean_var, modified to run in this
# script as a function and to accept the mean/var dtype as an input argument
def get_mean_var__dtype_mod(X,
        axis=0,
        running_mean_var_dtype= 'float64'):
    if sp.issparse(X):
        mean, var = sparse_mean_variance_axis(X, axis=axis, running_mean_var_dtype=running_mean_var_dtype)
    else:
        mean = np.mean(X, axis=axis, dtype=running_mean_var_dtype)
        mean_sq = np.multiply(X, X).mean(axis=axis, dtype=running_mean_var_dtype)
        var = mean_sq - mean**2
    # enforce R convention (unbiased estimator) for variance
    var *= X.shape[axis] / (X.shape[axis] - 1)
    return mean, var


# Copy of preprocessing/_simple.py, scale_array, modified to work with
# a specified dtype
def scale_array__dtype_mod(
    X,
    zero_center: bool = True,
    max_value: Optional[float] = None,
    copy: bool = False,
    return_mean_std: bool = False,
    working_dtype_float = 'float64',
    running_mean_var_dtype = 'float64'):
    if copy:
        X = X.copy()
    if not zero_center and max_value is not None:
        print( # Be careful of what? This should be more specific
            "... be careful when using `max_value` " "without `zero_center`.")
    if np.issubdtype(X.dtype, np.integer):
        print('... as scaling leads to float results, integer '
            'input is cast to float, returning copy.')
        X = X.astype(working_dtype_float)

    # Compute means and variances per column using higher
    # precision to avoid overflow during the accumulator
    # steps (where the max value of e.g. the sum can go much
    # higher than the values of individual elements)
    mean, var = get_mean_var__dtype_mod(X,
            running_mean_var_dtype = running_mean_var_dtype)
    # Set means and variances back to working precision (because
    # once computed, these *should* be less than the max element
    # value) for the scaling computation
    mean = mean.astype(working_dtype_float)
    var = var.astype(working_dtype_float)
    std = np.sqrt(var)
    std[std == 0] = 1
    if sp.issparse(X):
        if zero_center:
            # Zero-centering creates many non-zero elements,
            # so the matrix with these values shouldn't be sparse.
            print('Cannot zero-center sparse matrix: converting to dense for scaling.')
            X = X.copy().astype(working_dtype_float).todense()
            X -= mean
            X /= std
        else:
            sparsefuncs.inplace_column_scale(X, 1 / std)
    else:
        if zero_center:
            X -= mean
        X /= std
    # do the clipping
    if max_value is not None:
        print(f"... clipping at max_value {max_value}")
        X[X > max_value] = max_value
    if return_mean_std:
        return X.astype(working_dtype_float), mean, std
    else:
        return X.astype(working_dtype_float)

