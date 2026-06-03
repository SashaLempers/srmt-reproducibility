"""
NP-4 : Bareiss exact sur Q (en fait sur Z, fraction-free) pour Y_4(det3) et Y_4(perm3).
Matrice 1134 x 1134, entrees int dans {-1, 0, 1} ou tres petites.
Bareiss : O(n^3) operations, valeurs intermediaires explosent mais restent rationnelles.
"""
import sys
import json
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, '/home/user/workspace/srmt')

SRMT = Path('/home/user/workspace/srmt')

M_d = np.load(SRMT / 'np3_M_det_a4.npy')
M_p = np.load(SRMT / 'np3_M_perm_a4.npy')

print(f"M_det shape={M_d.shape}, nnz={np.count_nonzero(M_d)}")
print(f"M_perm shape={M_p.shape}, nnz={np.count_nonzero(M_p)}")
print(f"M_det max abs val: {np.abs(M_d).max()}")
print(f"M_perm max abs val: {np.abs(M_p).max()}")


def bareiss_rank(M):
    """Bareiss fraction-free, retourne le rang sur Q. Modifie M (utiliser copy avant)."""
    # M est array d'objects Python (int)
    nr, nc = M.shape
    A = [[int(x) for x in M[i]] for i in range(nr)]
    rank = 0
    prev_pivot = 1
    sign = 1
    row = 0
    for col in range(nc):
        if row >= nr:
            break
        # Cherche pivot
        pivot_row = -1
        for r in range(row, nr):
            if A[r][col] != 0:
                pivot_row = r
                break
        if pivot_row == -1:
            continue
        if pivot_row != row:
            A[row], A[pivot_row] = A[pivot_row], A[row]
            sign = -sign
        pivot = A[row][col]
        # Bareiss : pour r > row, A[r][c] = (pivot * A[r][c] - A[r][col] * A[row][c]) // prev_pivot
        for r in range(row + 1, nr):
            if A[r][col] == 0:
                continue
            arc = A[r][col]
            # Vectorisons par liste comprehension Python
            row_old = A[r]
            pivot_row_vec = A[row]
            new_row = [(pivot * row_old[c] - arc * pivot_row_vec[c]) // prev_pivot for c in range(nc)]
            A[r] = new_row
        prev_pivot = pivot
        rank += 1
        row += 1
    return rank


print("\n=== Bareiss exact Y_4(det3) ===")
t0 = time.time()
r_d = bareiss_rank(M_d.copy())
t1 = time.time()
print(f"  rank Q = {r_d} (en {t1-t0:.1f}s)")

print("\n=== Bareiss exact Y_4(perm3) ===")
t0 = time.time()
r_p = bareiss_rank(M_p.copy())
t1 = time.time()
print(f"  rank Q = {r_p} (en {t1-t0:.1f}s)")

print(f"\n=== VERDICT ===")
print(f"  rank Q(Y_4(det3))  = {r_d}")
print(f"  rank Q(Y_4(perm3)) = {r_p}")
print(f"  DIFFERS Q : {r_d != r_p}")
print(f"  Delta = {abs(r_d - r_p)}")

with open(SRMT / 'np4_bareiss_a4.json', 'w') as f:
    json.dump({
        'a': 4,
        'matrix_shape': list(M_d.shape),
        'rank_Q_det3': r_d,
        'rank_Q_perm3': r_p,
        'differs_Q': r_d != r_p,
        'delta': abs(r_d - r_p),
    }, f, indent=2)

print(f"\nSauvegarde np4_bareiss_a4.json")
