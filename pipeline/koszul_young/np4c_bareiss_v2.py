"""
Bareiss v2 : utilise sympy pour les operations exactes en grands entiers.
Plus lent mais correct.
"""
import sys
import time
import numpy as np
from pathlib import Path
import sympy as sp

SRMT = Path('/home/user/workspace/srmt')

M_d = np.load(SRMT / 'np3_M_det_a4.npy')
M_p = np.load(SRMT / 'np3_M_perm_a4.npy')


def bareiss_rank_v2(M_np):
    """Bareiss fraction-free version Python int pures (pas de numpy int64)."""
    nr, nc = M_np.shape
    # Conversion vers liste de listes Python int (illimite)
    A = [[int(M_np[i, j]) for j in range(nc)] for i in range(nr)]
    
    rank = 0
    prev_pivot = 1
    row = 0
    for col in range(nc):
        if row >= nr:
            break
        pivot_row = -1
        for r in range(row, nr):
            if A[r][col] != 0:
                pivot_row = r
                break
        if pivot_row == -1:
            continue
        if pivot_row != row:
            A[row], A[pivot_row] = A[pivot_row], A[row]
        pivot = A[row][col]
        # Verification : prev_pivot doit diviser exactement le determinant 2x2
        for r in range(row + 1, nr):
            if A[r][col] == 0:
                continue
            arc = A[r][col]
            for c in range(nc):
                num = pivot * A[r][c] - arc * A[row][c]
                if prev_pivot == 1:
                    A[r][c] = num
                else:
                    # Division exacte garantie en Bareiss
                    q, rem = divmod(num, prev_pivot)
                    if rem != 0:
                        # Erreur d'invariant Bareiss : potentiellement signe ou ordre
                        # On gere les nombres negatifs : divmod(-5, 3) = (-2, 1) en Python
                        # Verifier que c'est bien 0 quand on attend.
                        # Si non-zero, c'est un bug.
                        raise ValueError(f"Bareiss invariant violation: num={num}, prev_pivot={prev_pivot}, rem={rem} at row={r}, col={c}")
                    A[r][c] = q
        prev_pivot = pivot
        rank += 1
        row += 1
    return rank


print(f"M_det shape={M_d.shape}, max abs val: {np.abs(M_d).max()}")

print(f"\n=== Bareiss v2 (Python int) Y_4(det3) ===")
t0 = time.time()
r_d = bareiss_rank_v2(M_d)
t1 = time.time()
print(f"  rank Q = {r_d} (en {t1-t0:.1f}s)")

print(f"\n=== Bareiss v2 Y_4(perm3) ===")
t0 = time.time()
r_p = bareiss_rank_v2(M_p)
t1 = time.time()
print(f"  rank Q = {r_p} (en {t1-t0:.1f}s)")

print(f"\n=== VERDICT v2 ===")
print(f"  rank Q(Y_4(det3))  = {r_d}")
print(f"  rank Q(Y_4(perm3)) = {r_p}")
print(f"  DIFFERS Q : {r_d != r_p}")
