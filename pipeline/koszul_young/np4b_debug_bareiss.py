"""
Debug : verifier Bareiss en comparant a un rang fiable.
On utilise sympy.Matrix.rank() (rationnel exact) sur un sous-bloc.

Strategie : 
  1. Verifier d'abord que les ranks det/perm sont coherent en testant un solveur exact different : sympy.
  2. Si sympy confirme rank_det = 950, alors mon bareiss est buggue.
  3. Si sympy donne 975, alors les 4 primes mentaient (impossible statistiquement, donc tres improbable).
"""
import sys
import time
import numpy as np
from pathlib import Path
import sympy as sp

SRMT = Path('/home/user/workspace/srmt')

M_d = np.load(SRMT / 'np3_M_det_a4.npy')
M_p = np.load(SRMT / 'np3_M_perm_a4.npy')

print(f"M_det shape={M_d.shape}, dtype={M_d.dtype}")

# Sympy Matrix - mais 1134x1134 est trop gros pour sympy.rank() en raisonnable.
# Plan B : utiliser scipy.sparse + galois pour rang en F_p avec une lib externe.
# Plan C : refaire Bareiss en sympy avec Matrix.rref(), au moins sur un bloc reduit.

# D'abord testons : extraire les colonnes/lignes non-nulles, et reduire la taille du probleme.
nnz_rows_d = np.any(M_d != 0, axis=1)
nnz_cols_d = np.any(M_d != 0, axis=0)
print(f"M_det: {nnz_rows_d.sum()} non-zero rows, {nnz_cols_d.sum()} non-zero cols")
nnz_rows_p = np.any(M_p != 0, axis=1)
nnz_cols_p = np.any(M_p != 0, axis=0)
print(f"M_perm: {nnz_rows_p.sum()} non-zero rows, {nnz_cols_p.sum()} non-zero cols")

# On a 2520 nnz entries dans une 1134x1134, donc ~2 par row. 
# Une matrice tres creuse : beaucoup de lignes/cols pourraient etre nulles.

# Reduire a la sous-matrice non-zero
M_d_red = M_d[np.ix_(nnz_rows_d, nnz_cols_d)]
M_p_red = M_p[np.ix_(nnz_rows_p, nnz_cols_p)]
print(f"M_det reduced: {M_d_red.shape}")
print(f"M_perm reduced: {M_p_red.shape}")

# Rang sur Q via sympy.Matrix de la reduite - mais reste gros.
# Essayons via numpy.linalg.matrix_rank en float (informatif).
import numpy as np
print(f"\nRangs numpy.linalg (float64) :")
print(f"  M_det:  {np.linalg.matrix_rank(M_d.astype(np.float64))}")
print(f"  M_perm: {np.linalg.matrix_rank(M_p.astype(np.float64))}")

# Verification mod p avec deux nouvelles primes pour eliminer toute coincidence
import sys
sys.path.insert(0, '/home/user/workspace/srmt')
from np2_koszul_matrix import rank_dense_mod_p

for p in [104729, 5000011, 7919, 31337]:
    r_d = rank_dense_mod_p(M_d.copy(), p)
    r_p = rank_dense_mod_p(M_p.copy(), p)
    print(f"  p={p}: rk_det={r_d}, rk_perm={r_p}")

# Alors avec quel p Bareiss buggue-t-il ? Verifions Bareiss sur une matrice simple connue.
print(f"\n=== Test Bareiss sur matrice simple ===")
M_test = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.int64)  # rank 2
print(f"Test 3x3 rank 2:")
sys.path.insert(0, '/home/user/workspace/srmt')
from np4_bareiss_a4 import bareiss_rank
print(f"  bareiss: {bareiss_rank(M_test.copy())}")
print(f"  numpy: {np.linalg.matrix_rank(M_test.astype(float))}")

M_test2 = np.eye(5, dtype=np.int64)
print(f"Test eye(5) rank 5: bareiss = {bareiss_rank(M_test2.copy())}")

M_test3 = np.zeros((5, 5), dtype=np.int64)
M_test3[0, 0] = 1
M_test3[1, 1] = 1
M_test3[3, 3] = 1  # rank 3, mais avec ligne/col 2 vide
print(f"Test sparse rank 3: bareiss = {bareiss_rank(M_test3.copy())}")
