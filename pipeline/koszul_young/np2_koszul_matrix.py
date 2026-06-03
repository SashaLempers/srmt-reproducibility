"""
NP-2 : construction de la matrice Koszul Y_a(f) selon LO11 §1.2 Remark 4.1.3,
specialisee a delta=1 (cas d=3).

YF_a(f) : V^* (x) Lambda^a V  ->  V (x) Lambda^{a+1} V
   alpha (x) e_J  --maps to--  sum_{k not in J}  ( partial_k f restricted to V (as catalecticant on alpha) ) * sign_k(J)  (x)  e_k wedge e_J

Plus precisement, en bases :
   M[(j, K), (i, J)] = sum over k in {0..n-1} :
     [ alpha = e_i^* coefficient ] * [ catalecticant : df/dx_k transforme e_i en quel x_j coef ] * [coef Koszul]
     
La forme la plus simple, suivant L17 §7.8 / LO11 §1.2 :
Le bloc (J,K) de la matrice (vue comme matrice n x n a coef polynomiaux en x) est :
  - non nul ssi K = {k} cup J avec k not in J,
  - dans ce cas vaut sign(k -> J,K) * df/dx_k VU COMME element de S^2 V.
  
Pour delta=1, on doit aplatir df/dx_k (polynome de S^2 V) en matrice 9x9 (catalecticant) :
   Cat(df/dx_k)[i][j] = coef de "x_i x_j" dans df/dx_k.
   
ATTENTION conventions : on travaille avec f comme element de Sym^3 V.
   d^2 f / dx_k dx_i = polynome lineaire en x.
   coef de x_j dans d^2 f / dx_k dx_i = entier (independent de x).
Donc Cat(df/dx_k)[i,j] = coef de x_i x_j dans df/dx_k. Or
   df/dx_k = sum_{a<=b} c_{kab} x_a x_b   (forme triee).
   coef de "x_i x_j" en convention SYMETRIQUE : si i != j, c'est c_{k, min(i,j), max(i,j)}.
   si i = j : c'est c_{k,i,i}.
   
Donc Cat(df/dx_k) est juste la *matrice de coefs* avec ligne i, col j -> c_{k, sort(i,j)}.
Note : c'est SYMETRIQUE en (i,j).

Plus precisement, en utilisant la formule de polarisation symetrique :
   df/dx_k vu comme map V* -> V (Sym^2 V = applications symetriques) est notre H_k.
   H_k[i,j] = coef de x_i x_j dans df/dx_k avec convention sym.

Construction explicite :
- cols : pour chaque a-tuple J (combination), pour chaque i in {0..n-1} : index col = (i, J)
- rows : pour chaque (a+1)-tuple K, pour chaque j in {0..n-1} : index row = (j, K)
- M[(j,K), (i,J)] = sum over k : (K = sort({k} cup J)) * sign_k(J,K) * H_k[i, j]

Comme la condition K = sort({k} u J) determine k de maniere unique a partir de J, K (k = unique elt of K \ J),
le calcul est simple si J subset K (sinon 0).
"""
import json
import sys
import pickle
import itertools
import time
from pathlib import Path
from fractions import Fraction
import numpy as np

sys.path.insert(0, '/home/user/workspace/srmt')

n = 9
SRMT = Path('/home/user/workspace/srmt')

with open(SRMT / 'np1_partials.pkl', 'rb') as f:
    pdata = pickle.load(f)
det3_partials = pdata['det3_partials']
perm3_partials = pdata['perm3_partials']
# det3_partials[k] = dict {(i,j) trie : coef de x_i x_j dans d det3/dx_k}.


def hess_block(partials_k):
    """Matrice 9x9 H_k telle que H_k[i,j] = coef de x_i x_j dans df/dx_k (convention symetrique)."""
    H = np.zeros((n, n), dtype=np.int64)
    for (i, j), v in partials_k.items():
        if i == j:
            H[i, j] = v
        else:
            H[i, j] = v
            H[j, i] = v
    return H


def koszul_sign(k, J):
    """Signe pour insertion de k dans J trie : (-1)^{position d'insertion}.
       J est un tuple trie de longueur a, k not in J. K = sorted(J + (k,)) = (..., k inserse, ...).
       Le signe est (-1)^{nb d'elts de J plus petits que k} = position de k dans K."""
    pos = sum(1 for j in J if j < k)
    return 1 if pos % 2 == 0 else -1


def build_koszul_matrix(partials, a):
    """Construit la matrice Y_a(f) : V* (x) Lambda^a V -> V (x) Lambda^{a+1} V.
    
    Forme : (n * C(n, a+1)) lignes par (n * C(n, a)) colonnes.
    Index col : (i, J) ou i in [0,n), J in combinations(range(n), a).
    Index row : (j, K) ou j in [0,n), K in combinations(range(n), a+1).
    """
    # Precalcule Hess_k pour chaque k
    H = [hess_block(partials[k]) for k in range(n)]
    
    cols_J = list(itertools.combinations(range(n), a))
    rows_K = list(itertools.combinations(range(n), a+1))
    ncols = n * len(cols_J)
    nrows = n * len(rows_K)
    
    # Build sparse en COO
    rows_idx = []
    cols_idx = []
    vals = []
    
    idx_J = {J: jj for jj, J in enumerate(cols_J)}
    idx_K = {K: kk for kk, K in enumerate(rows_K)}
    
    for J_idx, J in enumerate(cols_J):
        # k parcourt les elements not in J
        J_set = set(J)
        for k in range(n):
            if k in J_set:
                continue
            K = tuple(sorted(J + (k,)))
            K_idx = idx_K[K]
            sgn = koszul_sign(k, J)
            Hk = H[k]
            # Pour chaque (i, j), si Hk[i,j] != 0, ajoute entree
            for i in range(n):
                col = i * len(cols_J) + J_idx
                for j in range(n):
                    v = Hk[i, j]
                    if v == 0:
                        continue
                    row = j * len(rows_K) + K_idx
                    rows_idx.append(row)
                    cols_idx.append(col)
                    vals.append(sgn * v)
    
    return nrows, ncols, rows_idx, cols_idx, vals


def rank_mod_p(nrows, ncols, rows_idx, cols_idx, vals, p):
    """Rang sur F_p via numpy dense int64."""
    M = np.zeros((nrows, ncols), dtype=np.int64)
    for r, c, v in zip(rows_idx, cols_idx, vals):
        M[r, c] = (M[r, c] + v) % p
    return rank_dense_mod_p(M, p)


def rank_dense_mod_p(M, p):
    """Elimination de Gauss mod p, M modifie en place."""
    M = M.copy()
    nrows, ncols = M.shape
    rank = 0
    row = 0
    for col in range(ncols):
        if row >= nrows:
            break
        # trouve pivot
        pivot = -1
        for r in range(row, nrows):
            if M[r, col] % p != 0:
                pivot = r
                break
        if pivot == -1:
            continue
        if pivot != row:
            M[[row, pivot]] = M[[pivot, row]]
        # normalise
        inv = pow(int(M[row, col]) % p, p - 2, p)
        M[row] = (M[row] * inv) % p
        # elimine
        for r in range(nrows):
            if r != row and M[r, col] % p != 0:
                factor = M[r, col] % p
                M[r] = (M[r] - factor * M[row]) % p
        row += 1
        rank += 1
    return rank


# Test pour a=1, 2, 3, 4
PRIME = 65521

results = []
for a in [1, 2, 3, 4]:
    ncols_J = len(list(itertools.combinations(range(n), a)))
    nrows_K = len(list(itertools.combinations(range(n), a+1)))
    dim_src = n * ncols_J
    dim_tgt = n * nrows_K
    prod = dim_src * dim_tgt
    print(f"\n=== a={a} ===  src=V*(x)L^{a}V dim={dim_src}, tgt=V(x)L^{a+1}V dim={dim_tgt}, prod={prod}")
    
    t0 = time.time()
    nr_d, nc_d, ri_d, ci_d, v_d = build_koszul_matrix(det3_partials, a)
    t1 = time.time()
    nr_p, nc_p, ri_p, ci_p, v_p = build_koszul_matrix(perm3_partials, a)
    t2 = time.time()
    
    r_d = rank_mod_p(nr_d, nc_d, ri_d, ci_d, v_d, PRIME)
    t3 = time.time()
    r_p = rank_mod_p(nr_p, nc_p, ri_p, ci_p, v_p, PRIME)
    t4 = time.time()
    
    print(f"  build det: {t1-t0:.2f}s, build perm: {t2-t1:.2f}s")
    print(f"  rank det: {r_d} (mod {PRIME}) in {t3-t2:.2f}s")
    print(f"  rank perm: {r_p} (mod {PRIME}) in {t4-t3:.2f}s")
    print(f"  DIFFERS: {r_d != r_p}")
    
    results.append({
        'a': a,
        'dim_src': dim_src,
        'dim_tgt': dim_tgt,
        'prod': prod,
        'rank_det_mod': r_d,
        'rank_perm_mod': r_p,
        'differs': r_d != r_p,
        'prime': PRIME,
    })

with open(SRMT / 'np2_koszul_principal.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n\n=== RESUME Koszul principal LO11 §1.2 sur Sym^3(C^9) ===")
for r in results:
    print(f"  a={r['a']}: rk_det={r['rank_det_mod']}, rk_perm={r['rank_perm_mod']}, DIFFERS={r['differs']}")
