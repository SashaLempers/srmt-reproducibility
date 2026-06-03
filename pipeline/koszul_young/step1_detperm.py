"""
Etape 1: Construction explicite de det_3 et perm_3 dans Sym^3(C^9).

C^9 = M_3(C), base e_{ij} indexée par (i,j) in {0,1,2}^2.
Sym^3(C^9) a dimension C(9+3-1, 3) = C(11,3) = 165.
Une base monomiale est { e_{i1,j1} * e_{i2,j2} * e_{i3,j3} } avec multi-index trié.

Convention: variables x_{ij} = e_{ij}. det_3 et perm_3 sont des polynômes homogènes de degré 3 en ces 9 variables.

det_3 = sum_{sigma in S_3} sgn(sigma) * x_{0,sigma(0)} * x_{1,sigma(1)} * x_{2,sigma(2)}
perm_3 = sum_{sigma in S_3} x_{0,sigma(0)} * x_{1,sigma(1)} * x_{2,sigma(2)}

On les exprime comme vecteurs dans Q^165, base = multiset monomes de degré 3 en 9 variables, ordonnée lexicographiquement.
"""
import json, itertools, hashlib
from fractions import Fraction

# 9 variables, indexées 0..8 via (i,j) -> 3*i+j
def idx(i,j): return 3*i+j

# Base monomiale de Sym^3(C^9) = multisets de taille 3 sur {0,...,8}
basis = []
for a in range(9):
    for b in range(a,9):
        for c in range(b,9):
            basis.append((a,b,c))
assert len(basis) == 165, len(basis)
basis_index = {m:k for k,m in enumerate(basis)}

def monomial_key(triple):
    return tuple(sorted(triple))

def build_poly(coeffs_dict):
    """coeffs_dict : {tuple sorted of 3 indices -> coef int}. Renvoie vecteur 165 entiers."""
    v = [0]*165
    for k,c in coeffs_dict.items():
        v[basis_index[k]] += c
    return v

from itertools import permutations

det3_terms = {}
perm3_terms = {}
for sigma in permutations(range(3)):
    triple = tuple(sorted([idx(i, sigma[i]) for i in range(3)]))
    sgn = 1
    # signature
    inv = sum(1 for i in range(3) for j in range(i+1,3) if sigma[i] > sigma[j])
    sgn = -1 if inv % 2 else 1
    det3_terms[triple] = det3_terms.get(triple, 0) + sgn
    perm3_terms[triple] = perm3_terms.get(triple, 0) + 1

det3 = build_poly(det3_terms)
perm3 = build_poly(perm3_terms)

# Sanity checks
nnz_det = sum(1 for x in det3 if x != 0)
nnz_perm = sum(1 for x in perm3 if x != 0)
sum_det = sum(det3)
sum_perm = sum(perm3)
print(f"det3 nnz = {nnz_det}, sum = {sum_det}")
print(f"perm3 nnz = {nnz_perm}, sum = {sum_perm}")
# det3: 6 termes, somme = 0 (3 positifs + 3 négatifs)
# perm3: 6 termes, somme = 6
assert nnz_det == 6 and sum_det == 0
assert nnz_perm == 6 and sum_perm == 6

# Évaluation sur matrice identité: det(I)=1, perm(I)=1
# On évalue à x_{ij} = delta_{ij}, i.e. x_0=1,x_4=1,x_8=1, sinon 0
def eval_poly(vec, vals):
    """vec dans Q^165, vals dans Q^9. Évalue le poly homogène de degré 3."""
    from math import factorial
    s = Fraction(0)
    for k,(a,b,c) in enumerate(basis):
        if vec[k] == 0: continue
        # le monome x_a x_b x_c (multiset) a un coefficient multinomial 1 sur cette base
        s += Fraction(vec[k]) * vals[a] * vals[b] * vals[c]
    return s

I = [0]*9
I[idx(0,0)] = 1; I[idx(1,1)] = 1; I[idx(2,2)] = 1
print(f"det3(I) = {eval_poly(det3, I)}  (attendu 1)")
print(f"perm3(I) = {eval_poly(perm3, I)}  (attendu 1)")
assert eval_poly(det3, I) == 1
assert eval_poly(perm3, I) == 1

# Évaluation sur matrice J = ones: det(J)=0, perm(J)=6
J = [1]*9
print(f"det3(J) = {eval_poly(det3, J)}  (attendu 0)")
print(f"perm3(J) = {eval_poly(perm3, J)}  (attendu 6)")
assert eval_poly(det3, J) == 0
assert eval_poly(perm3, J) == 6

# Archivage + SHA-256
def sha(v):
    return hashlib.sha256(json.dumps(v).encode()).hexdigest()

artifact = {
    'basis_size': 165,
    'det3': det3,
    'perm3': perm3,
    'det3_sha256': sha(det3),
    'perm3_sha256': sha(perm3),
    'basis_sha256': sha([list(b) for b in basis]),
}
with open('/home/user/workspace/srmt/step1_detperm.json','w') as f:
    json.dump(artifact, f, indent=2)
print("\nSHA256 det3 =", artifact['det3_sha256'])
print("SHA256 perm3 =", artifact['perm3_sha256'])
print("SHA256 basis =", artifact['basis_sha256'])
print("Step 1 OK")
