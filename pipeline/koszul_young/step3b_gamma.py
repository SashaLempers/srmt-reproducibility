"""
Construction concrete de Gamma_{lambda,mu}(f) comme matrice entiere via Pieri.
mu/lambda doit etre une bande horizontale (HS) de taille 3.
"""
import json, sys, time
from itertools import permutations
sys.path.insert(0, '/home/user/workspace/srmt')
from step2_schur import dim_Schur, SSYT, count_SSYT
from step3_pieri import horizontal_strips, mono_basis, mono_index, det3, perm3


def horizontal_strip_pairs(lam, mu):
    """Verifie si mu/lam est une bande horizontale et retourne les cases ajoutees, sinon None."""
    L = max(len(lam), len(mu))
    lam_e = list(lam) + [0]*(L - len(lam))
    mu_e = list(mu) + [0]*(L - len(mu))
    cells = []
    for i in range(L):
        if mu_e[i] < lam_e[i]: return None
        for j in range(lam_e[i], mu_e[i]):
            cells.append((i,j))
    for i in range(L-1):
        if mu_e[i+1] > lam_e[i]: return None
    return cells


def build_gamma_matrix(lam, mu, f_vec, n=9):
    """Construit Gamma_{lam,mu}(f) en base SSYT.
       f_vec: vecteur 165 entiers (coefs sur base monomiale 0-indexée mono_basis).
       Output: list of lists (rows = SSYT(mu,n), cols = SSYT(lam,n)).
    """
    cells_added = horizontal_strip_pairs(lam, mu)
    if cells_added is None or len(cells_added) != 3:
        raise ValueError(f"mu/lam not HS of size 3: lam={lam}, mu={mu}")
    
    lam_ssyts = list(SSYT(lam, n))
    mu_ssyts = list(SSYT(mu, n))
    lam_idx = {T:k for k,T in enumerate(lam_ssyts)}
    
    nrows = len(mu_ssyts)
    ncols = len(lam_ssyts)
    M = [[0]*ncols for _ in range(nrows)]
    
    for row_idx, T_prime in enumerate(mu_ssyts):
        T = []
        for i in range(len(lam)):
            T.append(tuple(T_prime[i][:lam[i]]) if lam[i] > 0 else ())
        T = tuple(T)
        added_vals = sorted(T_prime[i][j] for (i,j) in cells_added)
        if T not in lam_idx:
            continue
        col = lam_idx[T]
        a,b,c = added_vals[0]-1, added_vals[1]-1, added_vals[2]-1
        key = (a,b,c)
        coef = f_vec[mono_index[key]]
        M[row_idx][col] += coef
    return M, lam_ssyts, mu_ssyts


def rank_int_matrix_bareiss(M):
    """Rang exact sur Q via reduction par lignes en arithmetique exacte (Fraction)."""
    from fractions import Fraction
    if not M or not M[0]: return 0
    A = [[Fraction(x) for x in row] for row in M]
    R = len(A); C = len(A[0])
    rank = 0
    row = 0
    for col in range(C):
        if row >= R: break
        # pivot
        piv = -1
        for r in range(row, R):
            if A[r][col] != 0:
                piv = r; break
        if piv < 0: continue
        A[row], A[piv] = A[piv], A[row]
        for r in range(R):
            if r == row: continue
            if A[r][col] != 0:
                factor = A[r][col] / A[row][col]
                for c in range(col, C):
                    A[r][c] -= factor * A[row][c]
        row += 1
        rank += 1
    return rank


def rank_mod_p(M, p):
    """Rang modulo p via Gauss en F_p."""
    if not M or not M[0]: return 0
    A = [[x % p for x in row] for row in M]
    R = len(A); C = len(A[0])
    rank = 0; row = 0
    for col in range(C):
        if row >= R: break
        piv = -1
        for r in range(row, R):
            if A[r][col] != 0:
                piv = r; break
        if piv < 0: continue
        A[row], A[piv] = A[piv], A[row]
        inv = pow(A[row][col], p-2, p)
        for r in range(R):
            if r == row: continue
            if A[r][col] != 0:
                factor = (A[r][col] * inv) % p
                for c in range(col, C):
                    A[r][c] = (A[r][c] - factor * A[row][c]) % p
        row += 1; rank += 1
    return rank


if __name__ == "__main__":
    # Test trivial: Gamma_{(), (3)}
    M_det, _, _ = build_gamma_matrix((), (3,), det3, n=9)
    assert len(M_det) == 165 and len(M_det[0]) == 1
    M_perm, _, _ = build_gamma_matrix((), (3,), perm3, n=9)
    r_det = rank_int_matrix_bareiss(M_det)
    r_perm = rank_int_matrix_bareiss(M_perm)
    print(f"Gamma_(){{(3)}}: rank det={r_det}, perm={r_perm}")
    assert r_det == 1 and r_perm == 1
    
    # Test Pieri legit: Gamma_{(1),(4)}
    M_det2, _, _ = build_gamma_matrix((1,), (4,), det3, n=9)
    print(f"Gamma_(1)(4) shape: {len(M_det2)} x {len(M_det2[0])} (expect 495 x 9)")
    assert len(M_det2) == 495 and len(M_det2[0]) == 9
    M_perm2, _, _ = build_gamma_matrix((1,), (4,), perm3, n=9)
    r_det2 = rank_int_matrix_bareiss(M_det2)
    r_perm2 = rank_int_matrix_bareiss(M_perm2)
    print(f"Gamma_(1)(4): rank det={r_det2}, perm={r_perm2}")
    
    # Test Pieri legit: Gamma_{(1),(3,1)}
    M_det3, _, _ = build_gamma_matrix((1,), (3,1), det3, n=9)
    print(f"Gamma_(1)(3,1) shape: {len(M_det3)} x {len(M_det3[0])} (expect 990 x 9)")
    M_perm3, _, _ = build_gamma_matrix((1,), (3,1), perm3, n=9)
    r_det3 = rank_int_matrix_bareiss(M_det3)
    r_perm3 = rank_int_matrix_bareiss(M_perm3)
    print(f"Gamma_(1)(3,1): rank det={r_det3}, perm={r_perm3}")
    
    # Test Pieri legit: Gamma_{(1),(2,1,1)}? mu/lam HS de 3: (2,1,1) - (1) = cells (0,1),(1,0),(2,0). 
    # Check HS: mu_1 = 1 <= lam_0 = 1, mu_2 = 1 <= lam_1 = 0? NO -> not HS
    # try Gamma_{(1),(1,1,1,1)}: HS? (1,1,1,1)-(1) = (1,0),(2,0),(3,0). mu_1=1<=lam_0=1 OK, mu_2=1<=lam_1=0 FAIL
    # try Gamma_{(2),(5)}: HS yes, 3 cases en ligne 0.
    M_det4, _, _ = build_gamma_matrix((2,), (5,), det3, n=9)
    print(f"Gamma_(2)(5) shape: {len(M_det4)} x {len(M_det4[0])}")
    M_perm4, _, _ = build_gamma_matrix((2,), (5,), perm3, n=9)
    r_det4 = rank_int_matrix_bareiss(M_det4)
    r_perm4 = rank_int_matrix_bareiss(M_perm4)
    print(f"Gamma_(2)(5): rank det={r_det4}, perm={r_perm4}")
