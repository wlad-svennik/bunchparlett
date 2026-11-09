# Bunch–Parlett LDL decomposition

A small SymPy implementation of LDL decomposition for indefinite Hermitian matrices, using the Bunch–Parlett pivoting strategy described in Section 4.4.3 of *Matrix Computations*, fourth edition.

The function `LDL(M)` returns `L`, `D` and `P` with

$$
M = P^T L D L^* P,
$$

where `P` is a permutation matrix, `L` is unit lower triangular, and `D` is block diagonal with blocks of size 1 or 2. Here $L^*$ denotes the conjugate transpose.

## Example

Install SymPy, then run the following from the repository directory:

```python
from sympy import Matrix, simplify, zeros
from bunchparlett import LDL

M = Matrix([[4, 2, 0],
            [2, -1, 1],
            [0, 1, 3]])
L, D, P = LDL(M)

assert simplify(M - P.T * L * D * L.H * P) == zeros(3)
```

## Scope

This is a research prototype intended for small matrices. It follows the recursive factorisation directly and uses matrix inverses for the block calculations.

Inputs must be Hermitian, with entries whose magnitudes SymPy can compare. The code does not validate the Hermitian condition. Zero active submatrices are handled by returning identity factors.

The implementation is in [`bunchparlett.py`](bunchparlett.py).
