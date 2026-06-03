"""
NP-5 : validation de la construction Koszul.

Test 1 : f = x_0^3 (cube d'une seule variable). 
  Cette forme cubique a rang symetrique 1 (c'est w^3 avec w=e_0).
  D'apres LO11 Prop 4.1.1 : pour [phi] in v_d(PV) (i.e. rang sym = 1),
  rank Y_a(phi) <= ? Selon LO11 sec 1.2 : la matrice est en block-form et chaque bloc
  est la matrice K_n (Koszul) avec un coef ; le rang d'une rang-1 puissance.
  En fait pour phi = w^d : Y_a(w^d)(alpha (x) omega) = alpha(w)^delta * w^delta * w (x) (w wedge omega)
  Donc l'image est dans : w (x) w wedge Lambda^a V, donc rang <= C(n-1, a).
  Pour n=9, a=4 : C(8,4) = 70.

Test 2 : equivariance sous une permutation 9->9 simple.
"""
import sys
import time
import numpy as np
from pathlib import Path
from fractions import Fraction
import itertools

sys.path.insert(0, '/home/user/workspace/srmt')
from np2_koszul_matrix import build_koszul_matrix, rank_dense_mod_p
from np4d_bareiss_v3 import gauss_rank_sparse

# Pour matrices denses, on utilise rank modulaire p=65521 (suffisant pour sanity).

SRMT = Path('/home/user/workspace/srmt')
n = 9


def poly_partials(poly_dict):
    """Renvoie liste de 9 dicts {(b,c) trie: coef} representant df/dx_i."""
    result = []
    for i in range(n):
        out = {}
        for (a, b, c), v in poly_dict.items():
            mult = (a == i) + (b == i) + (c == i)
            if mult == 0:
                continue
            lst = [a, b, c]
            lst.remove(i)
            key = tuple(sorted(lst))
            out[key] = out.get(key, 0) + v * mult
        result.append({k: v for k, v in out.items() if v != 0})
    return result


# Test 1 : f = x_0^3
print("=== TEST 1 : f = x_0^3 ===")
f1 = {(0, 0, 0): 1}
p1 = poly_partials(f1)
print(f"  df/dx_0 = {p1[0]}")  # should be {(0,0): 3}
assert p1[0] == {(0, 0): 3}
for i in range(1, n):
    assert p1[i] == {}, f"df/dx_{i} should be 0, got {p1[i]}"
print("  partials OK")

a = 4
nr, nc, ri, ci, vs = build_koszul_matrix(p1, a)
M1 = np.zeros((nr, nc), dtype=np.int64)
for r, c, v in zip(ri, ci, vs):
    M1[r, c] += v
print(f"  shape: {M1.shape}, nnz: {np.count_nonzero(M1)}")
r1 = rank_dense_mod_p(M1.copy(), 65521)
print(f"  rank Y_4(x_0^3) = {r1}")
print(f"  predit borne sup : C(8,4)=70 (pour decomposable)")
print(f"  -> {'OK (<=70)' if r1 <= 70 else 'PROBLEME, doit etre <= 70'}")

# Test 2 : f = x_0^3 + x_1^3 (rang sym 2)
print("\n=== TEST 2 : f = x_0^3 + x_1^3 ===")
f2 = {(0, 0, 0): 1, (1, 1, 1): 1}
p2 = poly_partials(f2)
nr, nc, ri, ci, vs = build_koszul_matrix(p2, a)
M2 = np.zeros((nr, nc), dtype=np.int64)
for r, c, v in zip(ri, ci, vs):
    M2[r, c] += v
r2 = rank_dense_mod_p(M2.copy(), 65521)
print(f"  rank Y_4(x_0^3 + x_1^3) = {r2}")
print(f"  predit borne sup : 2 * 70 = 140 (LO Prop 4.1.1 avec t=1, r=2)")
print(f"  -> {'OK (<=140)' if r2 <= 140 else 'PROBLEME, doit etre <= 140'}")

# Test 3 : equivariance par swap x_0 <-> x_1.
# det3 sous swap des deux premieres lignes -> -det3 (signe), perm3 -> perm3.
# Comme on calcule le rang, le rang de Y_4(-f) = rang de Y_4(f). Donc rang reste identique.
print("\n=== TEST 3 : equivariance par swap x_0 <-> x_3 (echange lignes 0 et 1 de la matrice 3x3) ===")
# Note : x_{i,j} = idx 3*i+j. Echange ligne 0 et 1 echange (idx 0,1,2) <-> (idx 3,4,5).
# perm sigma sur indices : sigma(0)=3, sigma(1)=4, sigma(2)=5, sigma(3)=0, sigma(4)=1, sigma(5)=2, sigma(6)=6, etc.
sigma = [3, 4, 5, 0, 1, 2, 6, 7, 8]
import json
with open(SRMT / 'step1_detperm.json') as fh:
    s1 = json.load(fh)
det3_vec = s1['det3']
mono_basis = []
for i in range(n):
    for j in range(i, n):
        for k in range(j, n):
            mono_basis.append((i, j, k))

def apply_perm(vec, sigma):
    """Applique une permutation des variables au polynome."""
    out = {}
    for idx, v in enumerate(vec):
        if v == 0:
            continue
        (a, b, c) = mono_basis[idx]
        new = tuple(sorted((sigma[a], sigma[b], sigma[c])))
        out[new] = out.get(new, 0) + v
    return out

det3_swapped = apply_perm(det3_vec, sigma)
p_swap = poly_partials(det3_swapped)
nr, nc, ri, ci, vs = build_koszul_matrix(p_swap, a)
M_sw = np.zeros((nr, nc), dtype=np.int64)
for r, c, v in zip(ri, ci, vs):
    M_sw[r, c] += v
r_sw = rank_dense_mod_p(M_sw.copy(), 65521)
print(f"  rank Y_4(det3 swap rows 01) = {r_sw}")
print(f"  comparaison avec rank Y_4(det3) = 950")
print(f"  -> {'OK equivariance' if r_sw == 950 else 'PROBLEME, devrait etre 950'}")

# Test 4 : Cubique generique aleatoire (rang sym 9 attendu). 
# Pour comparer la valeur 950 a la borne attendue.
print("\n=== TEST 4 : cubique aleatoire (sanity) ===")
import random
random.seed(42)
f_rand = {m: random.randint(-3, 3) for m in mono_basis}
p_rand = poly_partials(f_rand)
nr, nc, ri, ci, vs = build_koszul_matrix(p_rand, a)
M_r = np.zeros((nr, nc), dtype=np.int64)
for r, c, v in zip(ri, ci, vs):
    M_r[r, c] += v
r_rand = rank_dense_mod_p(M_r.copy(), 65521)
print(f"  rank Y_4(random) = {r_rand}")
print(f"  shape: {M_r.shape}, max abs: {np.abs(M_r).max()}")
print(f"  note : pour cubique generique, attendu rang plein ou proche du min de dim")
