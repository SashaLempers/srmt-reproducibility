"""
Rang modulaire RAPIDE via numpy int32 + élimination Gauss optimisée par lignes.
Utilisable pour matrices entières jusqu'a ~10^8 cellules (en mod p < 2^15).

Strategie: p = 65521 (premier < 2^16). On stocke A en int64 (numpy), apres chaque
modif on prend modulo p. Vectorise par numpy.
"""
import numpy as np
import sys, json, time, os, hashlib
sys.path.insert(0, '/home/user/workspace/srmt')
from step2_schur import partitions, dim_Schur, SSYT
from step3_pieri import horizontal_strips, mono_basis, mono_index
from step3b_gamma import horizontal_strip_pairs

with open('/home/user/workspace/srmt/step1_detperm.json') as F:
    artifact = json.load(F)
det3 = artifact['det3']
perm3 = artifact['perm3']

P = 65521

def build_gamma_sparse_arrays(lam, mu, f_vec, n=9):
    """Construit Gamma_{lam,mu}(f) en format COO (row, col, val) pour vectorisation."""
    cells_added = horizontal_strip_pairs(lam, mu)
    if cells_added is None or len(cells_added) != 3:
        raise ValueError("not HS3")
    lam_ssyts = list(SSYT(lam, n))
    mu_ssyts = list(SSYT(mu, n))
    lam_idx = {T:k for k,T in enumerate(lam_ssyts)}
    
    rows = []; cols = []; vals = []
    for r_idx, T_prime in enumerate(mu_ssyts):
        T = []
        for i in range(len(lam)):
            T.append(tuple(T_prime[i][:lam[i]]) if lam[i] > 0 else ())
        T = tuple(T)
        if T not in lam_idx:
            continue
        c_idx = lam_idx[T]
        added_vals = sorted(T_prime[i][j] for (i,j) in cells_added)
        a,b,c = added_vals[0]-1, added_vals[1]-1, added_vals[2]-1
        coef = f_vec[mono_index[(a,b,c)]]
        if coef != 0:
            rows.append(r_idx); cols.append(c_idx); vals.append(coef)
    return np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32), np.array(vals, dtype=np.int64), len(mu_ssyts), len(lam_ssyts)


def rank_mod_p_dense_numpy(rows, cols, vals, R, C, p=P):
    """Construit la matrice dense in-place, fait Gauss en mod p, retourne rang."""
    # Memoire: R*C * 4 bytes en int32. Pour 1e8 cellules = 400 MB. OK pour <= 2e7. 
    # Le sandbox a 8GB.
    if R * C > 5e7:
        return None  # trop grand pour stockage dense
    A = np.zeros((R, C), dtype=np.int64)
    np.add.at(A, (rows, cols), vals)
    A %= p
    
    rank = 0
    row = 0
    for col in range(C):
        if row >= R: break
        # trouve pivot
        piv = -1
        for r in range(row, R):
            if A[r, col] != 0: piv = r; break
        if piv < 0: continue
        if piv != row:
            A[[row, piv]] = A[[piv, row]]
        inv = pow(int(A[row, col]), p-2, p)
        for r in range(R):
            if r == row: continue
            if A[r, col] != 0:
                factor = (A[r, col] * inv) % p
                A[r] = (A[r] - factor * A[row]) % p
        row += 1; rank += 1
    return rank


def rank_mod_p_sparse_smith(rows, cols, vals, R, C, p=P):
    """Approche: utiliser scipy sparse + LU mod p. Pas dispo direct. On fait Gauss
    sur la matrice dense, avec Numpy vectorise pour les operations sur lignes."""
    return rank_mod_p_dense_numpy(rows, cols, vals, R, C, p)


# Test rapide sur un petit candidat
if __name__ == "__main__":
    import time
    # (1,1) -> (4,1): 36 x 3168
    print("=== Test rapide (1,1) -> (4,1) ===")
    t0 = time.time()
    rd = build_gamma_sparse_arrays((1,1), (4,1), det3, n=9)
    rp = build_gamma_sparse_arrays((1,1), (4,1), perm3, n=9)
    print(f"build {time.time()-t0:.2f}s; nnz det={len(rd[2])} perm={len(rp[2])} shape={rd[3]}x{rd[4]}")
    t0 = time.time()
    rk_d = rank_mod_p_dense_numpy(*rd)
    print(f"rank det = {rk_d} in {time.time()-t0:.2f}s")
    t0 = time.time()
    rk_p = rank_mod_p_dense_numpy(*rp)
    print(f"rank perm = {rk_p} in {time.time()-t0:.2f}s")
    
    # Test gros: (3) -> (4,2): 165 x 10692 = 1.76M
    print("\n=== Test gros (3) -> (4,2) ===")
    t0 = time.time()
    rd = build_gamma_sparse_arrays((3,), (4,2), det3, n=9)
    rp = build_gamma_sparse_arrays((3,), (4,2), perm3, n=9)
    print(f"build {time.time()-t0:.2f}s; nnz det={len(rd[2])} perm={len(rp[2])} shape={rd[3]}x{rd[4]}")
    t0 = time.time()
    rk_d = rank_mod_p_dense_numpy(*rd)
    print(f"rank det = {rk_d} in {time.time()-t0:.2f}s")
    t0 = time.time()
    rk_p = rank_mod_p_dense_numpy(*rp)
    print(f"rank perm = {rk_p} in {time.time()-t0:.2f}s")
