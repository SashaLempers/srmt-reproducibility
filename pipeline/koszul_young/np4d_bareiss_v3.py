"""
Bareiss v3 : reformulation correcte pour rang d'une matrice rectangulaire/singuliere.

Approche : utiliser la methode de fraction-free Gauss-Jordan SANS l'optimisation Bareiss
(qui suppose pivotage diagonal sans skip). On garde des entiers et on divise par PGCD
de la ligne a la fin pour eviter l'explosion.

Alternative : sympy.Matrix.rank() qui utilise des algos optimises pour matrices creuses
(plus lent mais correct).

Pour 1134x1134 avec 2520 nnz tres creuse, sympy devrait gerer en temps raisonnable.
"""
import sys
import time
import numpy as np
from pathlib import Path

SRMT = Path('/home/user/workspace/srmt')

M_d = np.load(SRMT / 'np3_M_det_a4.npy')
M_p = np.load(SRMT / 'np3_M_perm_a4.npy')


def gauss_rank_exact(M_np):
    """Gauss-Jordan exact sur Q via Python int + GCD pour normaliser.
    Lent mais correct."""
    from math import gcd
    from fractions import Fraction
    nr, nc = M_np.shape
    A = [[Fraction(int(M_np[i, j])) for j in range(nc)] for i in range(nr)]
    
    rank = 0
    row = 0
    for col in range(nc):
        if row >= nr:
            break
        # Trouve pivot
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
        # Normalise ligne pivot
        A[row] = [x / pivot for x in A[row]]
        # Elimine
        for r in range(nr):
            if r == row:
                continue
            if A[r][col] != 0:
                factor = A[r][col]
                A[r] = [A[r][c] - factor * A[row][c] for c in range(nc)]
        row += 1
        rank += 1
    return rank


def gauss_rank_sparse(M_np):
    """Gauss exact mais en exploitant la sparsity : on stocke chaque ligne comme dict {col: Fraction}."""
    from fractions import Fraction
    nr, nc = M_np.shape
    # Listes de dicts
    rows = []
    for i in range(nr):
        d = {}
        for j in range(nc):
            v = int(M_np[i, j])
            if v != 0:
                d[j] = Fraction(v)
        rows.append(d)
    
    rank = 0
    used_rows = [False] * nr
    for col in range(nc):
        # Trouve une ligne non encore utilisee avec entree non-nulle en col
        pivot_r = -1
        for r in range(nr):
            if not used_rows[r] and col in rows[r]:
                pivot_r = r
                break
        if pivot_r == -1:
            continue
        used_rows[pivot_r] = True
        rank += 1
        pivot_val = rows[pivot_r][col]
        # Normalise pas necessaire si on n'a pas besoin du resultat reduit.
        # Mais on doit eliminer la colonne col des autres lignes (used ET unused).
        for r in range(nr):
            if r == pivot_r:
                continue
            if col not in rows[r]:
                continue
            factor = rows[r][col] / pivot_val
            # rows[r] = rows[r] - factor * rows[pivot_r]
            for c, v in rows[pivot_r].items():
                if c in rows[r]:
                    new_v = rows[r][c] - factor * v
                    if new_v == 0:
                        del rows[r][c]
                    else:
                        rows[r][c] = new_v
                else:
                    rows[r][c] = -factor * v
    return rank


print(f"M_det: shape={M_d.shape}, nnz={np.count_nonzero(M_d)}")

print(f"\n=== Gauss sparse Y_4(det3) ===")
t0 = time.time()
r_d = gauss_rank_sparse(M_d)
t1 = time.time()
print(f"  rank Q = {r_d} (en {t1-t0:.1f}s)")

print(f"\n=== Gauss sparse Y_4(perm3) ===")
t0 = time.time()
r_p = gauss_rank_sparse(M_p)
t1 = time.time()
print(f"  rank Q = {r_p} (en {t1-t0:.1f}s)")

print(f"\n=== VERDICT EXACT Q ===")
print(f"  rank Q(Y_4(det3))  = {r_d}")
print(f"  rank Q(Y_4(perm3)) = {r_p}")
print(f"  DIFFERS Q : {r_d != r_p}")
print(f"  Delta = {abs(r_d - r_p)}")

import json
with open(SRMT / 'np4d_gauss_exact.json', 'w') as f:
    json.dump({
        'a': 4,
        'method': 'Sparse Fraction Gauss-Jordan over Q',
        'matrix_shape': list(M_d.shape),
        'rank_Q_det3': r_d,
        'rank_Q_perm3': r_p,
        'differs_Q': r_d != r_p,
        'delta': abs(r_d - r_p),
    }, f, indent=2)
print("Sauvegarde np4d_gauss_exact.json")
