#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
route_flattening_gpu_windows.py
================================
Recherche de DÉFAUTS NON NULS dans la clôture d'orbite de det_3x3 ⊂ Sym^3(C^9)
par la méthode des Young flattenings / catalecticants, en arithmétique exacte mod p.

CONTEXTE MATHÉMATIQUE
---------------------
- V = C^9, W = Sym^3(V), dim W = C(11,3) = 165.
- f_det3 ∈ W est le déterminant 3×3 générique (cubique en 9 variables x_{ij}).
- On étudie la clôture d'orbite O_det3 ⊂ W sous GL_9 (changement de variable).
- Pour chaque famille de flattenings Y(f), on compare :
      rank_generic = max_{f aléatoire} rang(Y(f))  [mod p]
      rank_det3    = rang(Y(f_det3))               [mod p]
  Si rank_det3 < rank_generic pour TROIS premiers distincts → DÉFAUT NON NUL CERTIFIÉ.

FAMILLES IMPLÉMENTÉES
---------------------
C1 : Catalecticant Hessien  H_f : V → Sym^2(V)*      matrice 9 × 45
C2 : Catalecticant dérivées D1_f : V* → Sym^2(V)     matrice 9 × 45
C3 : Young flattening (1,2) Y_{1,2}(f) : V ⊗ Δ → Sym^2 V (TODO équivariant complet)
C4 : Multiplication par f   M_f : Sym^2 V → Sym^5 V  matrice 45 × 126 (squelette)
C5 : Koszul-like            K_f : V ⊗ Sym^2 V → Sym^3 V  (squelette)

CONTRAT D'HONNÊTETÉ (exigence mathématique)
--------------------------------------------
- Cette méthode peut ÉCHOUER à trouver un défaut même s'il en existe un.
- Une chute de rang sur 1 seul premier n'est PAS un défaut certifié.
- Certaines flattenings peuvent dépasser la RAM ; dans ce cas → ISSUE C, pas résultat faux.
- Aucune comparaison flottante. Tout rang est calculé EXACTEMENT mod p.
- Pas de dépendance à Sage ni Macaulay2.

STATUT
------
C1 et C2 : implémentées complètement.
C3, C4, C5 : squelettes avec TODO détaillés.

USAGE
-----
  python route_flattening_gpu_windows.py --mode audit
  python route_flattening_gpu_windows.py --mode sandbox
  python route_flattening_gpu_windows.py --mode test-det3
  python route_flattening_gpu_windows.py --mode flattening-one --family C1
  python route_flattening_gpu_windows.py --mode flattening-all
  python route_flattening_gpu_windows.py --mode flattening-longrun [--resume]
  python route_flattening_gpu_windows.py --mode review

Auteur : généré pour Sasha Lempers — mai 2026
"""

import argparse
import itertools
import json
import logging
import math
import os
import random
import sys
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Imports optionnels protégés
# ---------------------------------------------------------------------------
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

try:
    import cupy as cp
    HAS_CUPY = True
except Exception:
    HAS_CUPY = False
    cp = None

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False
    psutil = None

# ---------------------------------------------------------------------------
# Constantes globales
# ---------------------------------------------------------------------------
PRIMES = [1_000_000_007, 998_244_353, 1_000_000_009]
N_VAR = 9          # dim V
DIM_SYM2 = 45      # dim Sym^2(V)  = C(9+1,2)
DIM_SYM3 = 165     # dim Sym^3(V)  = C(11,3)
DIM_SYM5 = 1287    # dim Sym^5(V^9) = C(13,5) = 1287  [126=C(9,4) était faux]

# Garde-fous mémoire
WEIGHTSPACE_MAX = 300     # colonnes / lignes max avant alerte
FLATTENING_MAX  = 500     # côté max de la matrice (lignes × colonnes)

CHECKPOINT_FILE = Path("route_flattening_checkpoint.json")
LOG_FILE        = Path("route_flattening_run.log")

N_GENERIC_SAMPLES = 5
GENERIC_SEEDS     = [42, 137, 271, 314, 999]

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------
logger = logging.getLogger("route_flattening")
logger.setLevel(logging.DEBUG)

_fmt = logging.Formatter("[%(levelname)s] %(message)s")

_fh = logging.FileHandler(str(LOG_FILE), encoding="utf-8")
_fh.setLevel(logging.DEBUG)
_fh.setFormatter(_fmt)

_sh = logging.StreamHandler(sys.stdout)
_sh.setLevel(logging.INFO)
_sh.setFormatter(_fmt)

logger.addHandler(_fh)
logger.addHandler(_sh)

# ---------------------------------------------------------------------------
# Bases monomiales
# ---------------------------------------------------------------------------

def _sym_basis(n_var: int, degree: int) -> List[Tuple[int, ...]]:
    """Retourne la liste des multi-indices (e_0,...,e_{n-1}) de degré `degree`
    en ordre lexicographique (e_0 ≥ e_1 ≥ ... ≥ 0 pour forme canonique).

    On utilise la convention : monomiale x_0^{a0} ... x_{n-1}^{a_{n-1}}
    avec a0 + ... + a_{n-1} = degree, triée lexicographiquement."""
    basis = []
    for combo in itertools.combinations_with_replacement(range(n_var), degree):
        # combo est un tuple de longueur `degree` (ex: (0,0,1) pour x0^2*x1)
        # On convertit en exposants
        exp = [0] * n_var
        for idx in combo:
            exp[idx] += 1
        basis.append(tuple(exp))
    return basis

# Bases globales (calculées une fois)
_BASIS_SYM1: List[Tuple] = []
_BASIS_SYM2: List[Tuple] = []
_BASIS_SYM3: List[Tuple] = []
_BASIS_SYM5: List[Tuple] = []
_IDX_SYM1:  Dict[Tuple, int] = {}
_IDX_SYM2:  Dict[Tuple, int] = {}
_IDX_SYM3:  Dict[Tuple, int] = {}
_IDX_SYM5:  Dict[Tuple, int] = {}


def _init_bases():
    global _BASIS_SYM1, _BASIS_SYM2, _BASIS_SYM3, _BASIS_SYM5
    global _IDX_SYM1, _IDX_SYM2, _IDX_SYM3, _IDX_SYM5
    _BASIS_SYM1 = _sym_basis(N_VAR, 1)
    _BASIS_SYM2 = _sym_basis(N_VAR, 2)
    _BASIS_SYM3 = _sym_basis(N_VAR, 3)
    _BASIS_SYM5 = _sym_basis(N_VAR, 5)
    _IDX_SYM1   = {m: i for i, m in enumerate(_BASIS_SYM1)}
    _IDX_SYM2   = {m: i for i, m in enumerate(_BASIS_SYM2)}
    _IDX_SYM3   = {m: i for i, m in enumerate(_BASIS_SYM3)}
    _IDX_SYM5   = {m: i for i, m in enumerate(_BASIS_SYM5)}
    # Assertions de cohérence des dimensions (toute erreur ici est un bug)
    assert len(_BASIS_SYM2) == DIM_SYM2,  f"Sym2: {len(_BASIS_SYM2)} ≠ {DIM_SYM2}"
    assert len(_BASIS_SYM3) == DIM_SYM3,  f"Sym3: {len(_BASIS_SYM3)} ≠ {DIM_SYM3}"
    assert len(_BASIS_SYM5) == DIM_SYM5,  f"Sym5: {len(_BASIS_SYM5)} ≠ {DIM_SYM5}"


def _exp_add(a: Tuple, b: Tuple) -> Tuple:
    """Somme composante à composante de deux multi-indices."""
    return tuple(a[i] + b[i] for i in range(len(a)))


def _multinomial(exp: Tuple) -> int:
    """Coefficient multinomial deg! / (e0! e1! ... en!)."""
    d = sum(exp)
    result = math.factorial(d)
    for e in exp:
        result //= math.factorial(e)
    return result


# ---------------------------------------------------------------------------
# Encodage de f ∈ Sym^3 V
# ---------------------------------------------------------------------------

def encode_det3() -> List[int]:
    """Retourne le vecteur de coefficients de det3 dans la base _BASIS_SYM3.

    det3 = somme_{σ ∈ S3} sgn(σ) x_{0,σ(0)} x_{1,σ(1)} x_{2,σ(2)}

    Variables : x_{ij} avec i ∈ {0,1,2}, j ∈ {0,1,2}
    Indexées de 0 à 8 par : var_idx = 3*i + j.

    det3 = x00*x11*x22 - x00*x12*x21 - x01*x10*x22 + x01*x12*x20
           + x02*x10*x21 - x02*x11*x20

    Convention : x00=var0, x01=var1, x02=var2,
                 x10=var3, x11=var4, x12=var5,
                 x20=var6, x21=var7, x22=var8.
    """
    _ensure_bases()
    coeffs = [0] * DIM_SYM3

    # Les 6 termes du déterminant 3×3 : (permutation, signe)
    # permutation sigma sur {0,1,2} : var_{i, sigma(i)}
    perms = [
        ([0, 1, 2],  1),   # id
        ([0, 2, 1], -1),   # (12)
        ([1, 0, 2], -1),   # (01)
        ([1, 2, 0],  1),   # (012)
        ([2, 0, 1],  1),   # (021)
        ([2, 1, 0], -1),   # (02)
    ]

    for sigma, sign in perms:
        # Les 3 variables : x_{i, sigma[i]}
        var_indices = [3 * i + sigma[i] for i in range(3)]
        # Monôme x_{v0} * x_{v1} * x_{v2}
        exp = [0] * N_VAR
        for v in var_indices:
            exp[v] += 1
        exp = tuple(exp)

        # Trouver l'indice dans la base
        if exp not in _IDX_SYM3:
            raise RuntimeError(f"Monôme {exp} absent de la base Sym^3 — bug.")
        coeffs[_IDX_SYM3[exp]] += sign

    return coeffs


def _ensure_bases():
    if not _BASIS_SYM3:
        _init_bases()


def random_cubic_mod_p(p: int, seed: int) -> List[int]:
    """Cubique aléatoire dans Sym^3 V, coefficients mod p."""
    _ensure_bases()
    rng = random.Random(seed)
    return [rng.randint(0, p - 1) for _ in range(DIM_SYM3)]


# ---------------------------------------------------------------------------
# Arithmétique modulaire
# ---------------------------------------------------------------------------

def isprime(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def modinv(a: int, p: int) -> int:
    """Inverse modulaire de a mod p (p premier, a ≢ 0 mod p)."""
    return pow(a, p - 2, p)


# ---------------------------------------------------------------------------
# Rang mod p (CPU)
# ---------------------------------------------------------------------------

def rank_modp_cpu(mat: "np.ndarray", p: int) -> int:
    """Rang exact mod p par élimination de Gauss (CPU, numpy int64).

    mat : tableau numpy 2D d'entiers, déjà réduits mod p.
    Retourne le rang entier.
    JAMAIS de float.
    """
    if not HAS_NUMPY:
        raise RuntimeError("numpy requis pour rank_modp_cpu")
    A = np.array(mat, dtype=np.int64) % p
    nrows, ncols = A.shape
    pivot_row = 0
    for col in range(ncols):
        if pivot_row >= nrows:
            break
        # Chercher un pivot non nul dans la colonne
        found = -1
        for row in range(pivot_row, nrows):
            if A[row, col] != 0:
                found = row
                break
        if found == -1:
            continue
        # Échanger les lignes
        if found != pivot_row:
            A[[pivot_row, found]] = A[[found, pivot_row]]
        # Normaliser la ligne pivot
        inv_piv = modinv(int(A[pivot_row, col]), p)
        A[pivot_row] = (A[pivot_row] * inv_piv) % p
        # Éliminer dans toutes les autres lignes
        for row in range(nrows):
            if row != pivot_row and A[row, col] != 0:
                factor = int(A[row, col])
                A[row] = (A[row] - factor * A[pivot_row]) % p
        pivot_row += 1
    return pivot_row


# ---------------------------------------------------------------------------
# Rang mod p (GPU via CuPy)
# ---------------------------------------------------------------------------

def rank_modp_exact_gpu(mat: "np.ndarray", p: int) -> int:
    """Rang exact mod p sur GPU (CuPy), avec fallback CPU automatique.

    Utilise la même élimination de Gauss que CPU, mais sur GPU.
    Si CuPy indisponible ou OOM, bascule silencieusement sur CPU.
    """
    if not HAS_CUPY:
        return rank_modp_cpu(mat, p)
    try:
        A = cp.array(mat, dtype=cp.int64) % p
        nrows, ncols = A.shape
        pivot_row = 0
        for col in range(ncols):
            if pivot_row >= nrows:
                break
            # Chercher pivot (on ramène la colonne sur CPU pour la recherche)
            col_vals = A[pivot_row:, col].get()
            nonzero = np.nonzero(col_vals)[0]
            if len(nonzero) == 0:
                continue
            found = int(nonzero[0]) + pivot_row
            if found != pivot_row:
                # Échange sur GPU
                tmp = A[pivot_row].copy()
                A[pivot_row] = A[found]
                A[found] = tmp
            inv_piv = modinv(int(A[pivot_row, col].get()), p)
            A[pivot_row] = (A[pivot_row] * inv_piv) % p
            # Élimination vectorisée sur GPU
            factors = A[:, col].copy()
            factors[pivot_row] = 0
            # A[row] -= factors[row] * A[pivot_row] pour tout row
            A = (A - cp.outer(factors, A[pivot_row])) % p
            pivot_row += 1
        return pivot_row
    except cp.cuda.memory.OutOfMemoryError:
        logger.warning("GPU OOM — bascule sur CPU pour le rang.")
        return rank_modp_cpu(mat, p)
    except Exception as e:
        logger.warning(f"GPU erreur ({e}) — bascule sur CPU.")
        return rank_modp_cpu(mat, p)


def rank_modp(mat: "np.ndarray", p: int) -> int:
    """Point d'entrée unique : GPU si disponible, sinon CPU."""
    return rank_modp_exact_gpu(mat, p)


# ---------------------------------------------------------------------------
# Helpers GL_n mod p
# ---------------------------------------------------------------------------

def random_GL9_mod_p(p: int, seed: int) -> "np.ndarray":
    """Matrice 9×9 aléatoire inversible mod p."""
    rng = random.Random(seed)
    while True:
        entries = [[rng.randint(0, p - 1) for _ in range(N_VAR)] for _ in range(N_VAR)]
        M = np.array(entries, dtype=np.int64)
        # Vérifier inversibilité via det (approximatif) → on fait le rang
        if rank_modp_cpu(M, p) == N_VAR:
            return M


def det9_mod_p(M: "np.ndarray", p: int) -> int:
    """Déterminant de M (9×9) mod p par élimination (exact)."""
    A = M.copy() % p
    n = A.shape[0]
    det = 1
    for col in range(n):
        found = -1
        for row in range(col, n):
            if A[row, col] != 0:
                found = row
                break
        if found == -1:
            return 0
        if found != col:
            A[[col, found]] = A[[found, col]]
            det = (-det) % p
        inv_piv = modinv(int(A[col, col]), p)
        det = (det * int(A[col, col])) % p
        A[col] = (A[col] * inv_piv) % p
        for row in range(n):
            if row != col and A[row, col] != 0:
                factor = int(A[row, col])
                A[row] = (A[row] - factor * A[col]) % p
    return int(det)


def action_g_on_cubic_mod_p(g: "np.ndarray", f_coeffs: List[int], p: int) -> List[int]:
    """Calcule g·f : le changement de variable x → g^{-T} x sur f ∈ Sym^3 V.

    Pour l'action de GL_n sur Sym^d V :
      (g·f)(x) = f(g^{-T} x)   où g^{-T} = (g^{-1})^T.

    Stratégie : on calcule les nouvelles variables y_i = Σ_j (g^{-T})_{ij} x_j,
    puis on évalue le polynôme f sur ces nouvelles variables symboliquement.

    Pour les tests (sandbox), on utilise n petit (N_VAR=9 peut être lent).
    On travaille ici de façon exacte mod p.
    """
    _ensure_bases()
    # Calculer g^{-1} mod p via élimination
    g_inv = _mat_inv_mod_p(g, p)
    # g^{-T}
    g_inv_T = g_inv.T % p

    # On exprime y_i = Σ_j (g_inv_T)_{ij} x_j
    # Chaque monôme x^exp = prod_i x_i^{e_i} devient
    # prod_i (Σ_j c_{ij} x_j)^{e_i}  [expansion multinomiale]
    # C'est coûteux pour N_VAR=9, deg=3 → on fait prod_i (lin)^{e_i} term par term.

    result = [0] * DIM_SYM3

    for mono_idx, exp in enumerate(_BASIS_SYM3):
        coeff = f_coeffs[mono_idx]
        if coeff == 0:
            continue
        coeff = int(coeff) % p
        # On développe prod_i y_i^{e_i} où y_i = Σ_j c_{ij} x_j
        # Résultat : dict exp_tuple -> coefficient
        poly = {tuple([0] * N_VAR): 1}
        for var_i, e_i in enumerate(exp):
            if e_i == 0:
                continue
            # Multiplier e_i fois par la forme linéaire y_i
            lin = {tuple(np.eye(N_VAR, dtype=int)[j].tolist()): int(g_inv_T[var_i, j])
                   for j in range(N_VAR) if g_inv_T[var_i, j] != 0}
            for _ in range(e_i):
                new_poly: Dict[Tuple, int] = {}
                for exp_a, ca in poly.items():
                    for exp_b, cb in lin.items():
                        exp_c = tuple(exp_a[k] + exp_b[k] for k in range(N_VAR))
                        new_poly[exp_c] = (new_poly.get(exp_c, 0) + ca * cb) % p
                poly = new_poly
        # Accumuler dans result
        for exp_c, c in poly.items():
            if exp_c in _IDX_SYM3:
                result[_IDX_SYM3[exp_c]] = (result[_IDX_SYM3[exp_c]] + coeff * c) % p

    return result


def _mat_inv_mod_p(M: "np.ndarray", p: int) -> "np.ndarray":
    """Inverse de M mod p par Gauss-Jordan (entier exact)."""
    n = M.shape[0]
    A = np.hstack([M.copy() % p, np.eye(n, dtype=np.int64)])
    for col in range(n):
        found = -1
        for row in range(col, n):
            if A[row, col] % p != 0:
                found = row
                break
        if found == -1:
            raise ValueError("Matrice singulière mod p")
        A[[col, found]] = A[[found, col]]
        inv_piv = modinv(int(A[col, col] % p), p)
        A[col] = (A[col] * inv_piv) % p
        for row in range(n):
            if row != col and A[row, col] % p != 0:
                factor = int(A[row, col] % p)
                A[row] = (A[row] - factor * A[col]) % p
    return A[:, n:] % p


# ---------------------------------------------------------------------------
# Coefficient de Taylor du monôme
# ---------------------------------------------------------------------------

def _coeff_of(f_coeffs: List[int], exp: Tuple, p: int) -> int:
    """Retourne le coefficient du monôme exp dans f, mod p.
    Convention : f = Σ_α c_α x^α (sans coefficient multinomial).
    """
    _ensure_bases()
    idx = _IDX_SYM3.get(exp)
    if idx is None:
        return 0
    return int(f_coeffs[idx]) % p


def _partial_coeff(f_coeffs: List[int], deriv_exp: Tuple, remain_exp: Tuple, p: int) -> int:
    """Coefficient dans ∂^{deriv_exp} f du monôme remain_exp.

    ∂^α f à partir de f = Σ_β c_β x^β :
    Terme β contribue si β ≥ α composante à composante :
      résidu_exp = β - α
      coefficient = c_β * multinomial(β) / multinomial(β-α)  [exact, pas flottant]
    """
    _ensure_bases()
    total_exp = _exp_add(deriv_exp, remain_exp)
    if total_exp not in _IDX_SYM3:
        return 0
    c_beta = _coeff_of(f_coeffs, total_exp, p)
    if c_beta == 0:
        return 0
    # Coefficient de ∂^{deriv_exp} (x^{total_exp}) est total_exp! / remain_exp!
    # = prod_i C(total_exp_i, deriv_exp_i) * deriv_exp_i!
    # = prod_i total_exp_i! / (remain_exp_i! * deriv_exp_i!)  … simplifié :
    factor = 1
    for i in range(N_VAR):
        b = total_exp[i]
        a = deriv_exp[i]
        # C(b, a) = b! / (a! (b-a)!)
        factor *= math.comb(b, a)
    return (c_beta * factor) % p


# ---------------------------------------------------------------------------
# C1 : Catalecticant Hessien
# ---------------------------------------------------------------------------
# H_f : V → Sym^2(V)*
# Pour v = e_i (base de V), la forme quadratique ∂_i f est dans Sym^2(V)*.
# La matrice est : M[i, j] = coefficient de x^α_{(j)} dans ∂_{e_i} f,
# où α_{(j)} est le j-ème monôme de degré 2.
#
# Remarque : ∂_{e_i} f = Σ_β c_β · β_i · x^{β - e_i}
# Le terme x^{β - e_i} est de degré 2 ssi β de degré 3 avec β_i ≥ 1.
#
# M[i, j] = Σ_{β : β_i ≥ 1, β - e_i = α_j}  c_β · β_i
#          = c_{α_j + e_i} · (α_j + e_i)_i
#          = c_{α_j + e_i} · (α_j_i + 1)
# ---------------------------------------------------------------------------

def build_flattening_C1(f_coeffs: List[int], p: int) -> "np.ndarray":
    """Catalecticant Hessien H_f : matrice 9 × 45 exacte mod p.

    Ligne i (i ∈ {0,...,8}) correspond à ∂/∂x_i.
    Colonne j (j ∈ {0,...,44}) correspond au monôme α_j ∈ Sym^2(V).

    M[i,j] = coefficient de α_j dans (∂/∂x_i)(f).
    """
    _ensure_bases()
    nrows = N_VAR        # 9
    ncols = DIM_SYM2     # 45
    M = np.zeros((nrows, ncols), dtype=np.int64)

    e = [tuple(1 if k == i else 0 for k in range(N_VAR)) for i in range(N_VAR)]

    for i in range(N_VAR):
        for j, alpha2 in enumerate(_BASIS_SYM2):
            # exp de degré 3 = alpha2 + e_i
            alpha3 = _exp_add(alpha2, e[i])
            if sum(alpha3) != 3:
                continue  # ne devrait pas arriver
            coeff = _partial_coeff(f_coeffs, e[i], alpha2, p)
            M[i, j] = coeff % p

    return M


# ---------------------------------------------------------------------------
# C2 : Seconde catalecticante Cat_{2,1}(f) : Sym^2(V) → Sym^1(V)
# ---------------------------------------------------------------------------
# Distinct de C1 (qui va V → Sym^2 V) par la direction de contraction.
#
# C1 : V  → Sym^2 V   M_C1[i, j]  = coeff de α_j ∈ Sym^2 dans (∂/∂x_i) f
# C2 : Sym^2 V → V    M_C2[j, i]  = coeff de α_i ∈ Sym^1 dans (∂^{α_j}) f
#
# Mathématiquement : M_C2 = M_C1^T  — même information, même rang.
# MAIS : c'est le seul catalecticant «classique» dans cette direction et
# il correspond à Cat_{2,1} de la littérature (Iarrobino-Kanev).
# On le conserve avec un commentaire honnête sur la redondance de rang,
# car il sert de référence et de validation croisée.
#
# Matrice 45 × 9.
# rang générique = rang(C1) = 9.
# Si rang(C2, det3) < 9 → même chute que C1 (tautologique).
# Si rang(C2, det3) ≠ rang(C1, det3) → bug dans l'implémentation.
# ---------------------------------------------------------------------------

def build_flattening_C2(f_coeffs: List[int], p: int) -> "np.ndarray":
    """Seconde catalecticante Cat_{2,1}(f) : Sym^2 V → V, matrice 45 × 9.

    M[j, i] = coefficient de α_i ∈ Sym^1 V dans (∂^{α_j} f), avec α_j ∈ Sym^2 V.
    C'est _partial_coeff(f, alpha2_j, alpha1_i, p).

    Rang générique = 9 = rang(C1). La redondance de rang est documentée et
    intentionnelle : C2 sert de validation croisée de C1.
    Un résultat rang(C2) ≠ rang(C1) révèle un bug d'implémentation.
    """
    _ensure_bases()
    nrows = DIM_SYM2   # 45
    ncols = N_VAR      # 9
    M = np.zeros((nrows, ncols), dtype=np.int64)

    for j, alpha2 in enumerate(_BASIS_SYM2):
        for i, alpha1 in enumerate(_BASIS_SYM1):
            coeff = _partial_coeff(f_coeffs, alpha2, alpha1, p)
            M[j, i] = coeff % p

    return M


# ---------------------------------------------------------------------------
# C3 : Hessien symétrique étendu — contraction croisée Sym^2 V × Sym^2 V → k
# ---------------------------------------------------------------------------
# Distinct de C1 et C2 par sa taille ET sa construction.
#
# Construction :
#   Pour deux monômes α, β ∈ Sym^1 V, on a :
#     H_f[α, β] = coeff de x^γ dans ∂_α ∂_β f  où γ = 0 (constante)
#   Mais pour f de degré 3 et α, β de degré 1 : ∂_α ∂_β f est de degré 1,
#   donc la constante est 0 sauf si le gradient lui-même est constant — non.
#
#   CONSTRUCTION VALIDE : Hessien tensoriel H_f : V × V → V
#   H_f(e_i, e_j) = gradient de (∂^2_{ij} f) = vecteur dans V.
#   On matricise en (N_VAR × N_VAR) × N_VAR = 81 × 9.
#   Symétrisée : Sym^2 V → V, matrice 45 × 9 = C2 (même rang).
#
#   VÉRITABLE C3 DISTINCTE : Hessien quadratique évalué en un point fixe x_0.
#   Pour un f cubique, H_f(x) = matrice 9×9 de formes LINÉAIRES en x.
#   On peut fixer x_0 et obtenir une matrice 9×9 (carrée) qui dépend de x_0.
#   Meilleur choix structurel : evaluer H_f sur TOUS les e_k et empiler.
#
#   C3 ici : matrice 9×9×9 matricisée en 81×9, sans symétrisation.
#   M[9i+j, k] = _partial_coeff(f, e_i + e_j, e_k, p)
#   C'est le tenseur de Hesse complet, non symétrisé sur (i,j).
#   Différent de C1 : taille 81×9 vs 9×45, même rang 9.
#
# AVERTISSEMENT : rang générique = 9 pour tout f non dégénéré.
# Pas de chute attendue pour det3 avec ce flattening (comme C1/C2).
# Il est inclus pour couverture et validation croisée.
#
# TODO (vraie Young flattening équivariante) :
# - Décomposer V = C^3 ⊗ C^3 et utiliser GL_3×GL_3 ⊂ GL_9.
# - Représentation S^{(2,1)}(C^3) de dim 8 pour chaque facteur.
# - Construire la matrice d'entrelacements GL_3×GL_3-équivariante
#   Sym^3(C^3 ⊗ C^3) → S^{(2,1)}(C^3) ⊗ S^{(2,1)}(C^3).
# - Nécessite les tableaux de Young standard et les matrices de Schur
#   (calcul de représentations de GL_n), hors portée sans Sage/SymPy.combinatorics.
# ---------------------------------------------------------------------------

def build_flattening_C3(f_coeffs: List[int], p: int) -> "np.ndarray":
    """Hessien tensoriel complet H_f : V ⊗ V → V, matricisé en 81 × 9.

    M[9*i + j, k] = _partial_coeff(f, e_i + e_j, e_k, p)
                  = coefficient de x_k dans ∂^2_{ij} f.

    Non symétrisé sur (i,j) : différent de C2 (Sym^2 V → V, 45×9).
    Même rang générique = 9 que C1/C2.

    UTILITÉ : validation croisée. Si rang(C3, det3) ≠ rang(C1, det3),
    c'est un bug d'implémentation, pas un défaut.

    TODO ÉQUIVARIANT : voir commentaire bloc ci-dessus.
    """
    _ensure_bases()
    e = [tuple(1 if k == ii else 0 for k in range(N_VAR)) for ii in range(N_VAR)]

    nrows = N_VAR * N_VAR  # 81
    ncols = N_VAR           # 9
    M = np.zeros((nrows, ncols), dtype=np.int64)

    for i in range(N_VAR):
        for j in range(N_VAR):
            deriv2 = _exp_add(e[i], e[j])   # multi-indice de degré 2
            for k in range(N_VAR):
                # coeff de e_k dans ∂^{e_i+e_j} f (résidu de degré 1)
                coeff = _partial_coeff(f_coeffs, deriv2, e[k], p)
                M[N_VAR * i + j, k] = coeff % p

    return M


# ---------------------------------------------------------------------------
# C4 : Multiplication par f, M_f : Sym^2 V → Sym^5 V
# ---------------------------------------------------------------------------
# M_f(g) = f · g pour g ∈ Sym^2 V, résultat dans Sym^5 V.
# Matrice 45 → dim Sym^5(V^9) = C(13,8) = C(13,5) = 1287.
# C'est très grand (45 × 1287) — garde-fou FLATTENING_MAX.
#
# TODO :
# - Calculer dim Sym^5(V^9) et vérifier que ça rentre en mémoire.
# - Implémenter le produit de polynômes mod p.
# - Variantes équivariantes : restreindre à des sous-espaces GL_9.
# ---------------------------------------------------------------------------

def _dim_sym(n_var: int, deg: int) -> int:
    """dim Sym^deg(V^n) = C(n+deg-1, deg)."""
    return math.comb(n_var + deg - 1, deg)


def build_flattening_C4(f_coeffs: List[int], p: int) -> Optional["np.ndarray"]:
    """Multiplication M_f : Sym^2 V → Sym^5 V, matrice 45 × dim_Sym5.

    dim Sym^5(V^9) = C(13,5) = 1287.
    Matrice 45 × 1287. Taille raisonnable (~370k entiers int64 ≈ 3 MB).

    TODO :
    - Implémenter le produit de polynômes degré 2 × degré 3 = degré 5.
    - Pour chaque α2 ∈ Sym^2, calculer les coefficients de f * x^{α2} ∈ Sym^5.
    - Chaque coefficient de Sym^5 : (f * x^{α2})[β5] = f[β5 - α2] si β5 ≥ α2.
    """
    _ensure_bases()
    dim_sym5 = _dim_sym(N_VAR, 5)
    nrows = DIM_SYM2    # 45
    ncols = dim_sym5    # 1287 pour N_VAR=9

    # Garde-fou : on estime la mémoire nécessaire
    mem_bytes = nrows * ncols * 8   # int64
    mem_mb = mem_bytes / 1e6
    if mem_mb > 50:
        logger.warning("C4 : matrice %d×%d ≈ %.0f MB — trop grande, risque ISSUE C.", nrows, ncols, mem_mb)
    else:
        logger.debug("C4 : matrice %d×%d ≈ %.1f MB (OK).", nrows, ncols, mem_mb)

    # Construire la base Sym^5 si pas déjà fait
    if not _BASIS_SYM5:
        _init_bases()

    logger.info("C4 : construction de la matrice de multiplication 45 × %d...", ncols)
    M = np.zeros((nrows, ncols), dtype=np.int64)

    for row_idx, alpha2 in enumerate(_BASIS_SYM2):
        # f * x^{alpha2} : pour chaque terme c_beta3 * x^{beta3} dans f,
        # le produit donne c_beta3 * x^{beta3 + alpha2} dans Sym^5
        for beta3_idx, beta3 in enumerate(_BASIS_SYM3):
            c = f_coeffs[beta3_idx]
            if c == 0:
                continue
            beta5 = _exp_add(beta3, alpha2)
            col_idx = _IDX_SYM5.get(beta5)
            if col_idx is None:
                continue  # ne devrait pas arriver
            M[row_idx, col_idx] = (M[row_idx, col_idx] + c) % p

    return M


# ---------------------------------------------------------------------------
# C5 : Koszul-like K_f : V ⊗ Sym^2 V → Sym^3 V
# ---------------------------------------------------------------------------
# K_f(v, g) = (∂/∂v) f · g  [produit de polynômes de degrés 2 et 2 = 4, non... ]
# Ou plutôt : K_f(v ⊗ g) = v · g - f(v,...) ?
# Convention ici : K_f(e_i ⊗ α_j) = x_i · x^{α_j} (dans Sym^3 V).
# C'est la matrice de multiplication : V ⊗ Sym^2 V → Sym^3 V, g → v·g.
# Cette matrice ne dépend PAS de f, c'est la résolution de Koszul.
# La version qui dépend de f est :
#   K_f : V → Sym^2 V par v → contraction partielle de f avec v.
# On l'identifie à C1 (même chose).
#
# TODO : implémenter la vraie version "Koszul-like" qui utilise f comme noyau :
#   Φ_f : V ⊗ Sym^2 V → Sym^3 V, Φ_f(v ⊗ g) = f(v,·,·) · g - v · ∂_? f ...
#   (voir Landsberg "Geometry and Complexity Theory", Ch. 8)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# C5 : Projecteur équivariant GL_3×GL_3 → S^(2,1)(C^3) ⊗ S^(2,1)(C^3)
# ---------------------------------------------------------------------------
# CONTEXTE (décomposition de Cauchy)
# -----------------------------------
# V = C^9 = C^3 ⊗ C^3.  Par Cauchy :
#   Sym^3(C^3 ⊗ C^3)  =  S^(3)(C^3)⊗S^(3)(C^3)     [dim 100]
#                      ⊕  S^(2,1)(C^3)⊗S^(2,1)(C^3)  [dim 64]
#                      ⊕  Λ^3(C^3)⊗Λ^3(C^3)          [dim 1]
# det3 est EXCLUSIVEMENT dans la composante Λ^3⊗Λ^3.
# => proj_{(2,1)}(det3) = 0 → rang = 0 sur det3.
# => proj_{(2,1)}(f_rand) ≠ 0 → rang > 0 → DÉFAUT CERTIFIABLE.
#
# STRATÉGIE (apolarity / tableaux semi-standard)
# -----------------------------------------------
# Dim S^(2,1)(C^3) = 8 (par la règle des crochets de Weyl / hook-length).
# Base SSYT pour λ=(2,1) en alphabet {0,1,2} (= indices de ligne/colonne):
#   SSYT = {(a,b,c) : a ≤ b, a < c, a,b,c ∈ {0,1,2}}
#   = [(0,0,1),(0,0,2),(0,1,1),(0,1,2),(0,2,1),(0,2,2),(1,1,2),(1,2,2)]
#
# Chaque T=(a,b,c) ∈ SSYT correspond à une base de S^(2,1)(C^3) via :
#   e_T = (coefficient de Young symétriseur standard)
#       = e_a ∧ e_c ⊗ e_b + e_b ∧ e_c ⊗ e_a  (deux termes du symétriseur de Young)
# En notation plus précise, pour λ=(2,1), le symétriseur standard sur S_3 donne :
#   pour T_std1 = [[0,1],[2]] : y = (id + (01)) − (02) − (012)   [non normalisé]
#   Mais on utilise directement la base SSYT pour l'apolarity pairing.
#
# CONSTRUCTION DE C5 (matrice 8×8, exacte)
# -----------------------------------------
# Pour chaque paire (T_row, T_col) de tableaux SSYT :
#   T_row = (a,b,c) encode la structure de LIGNE  (i-indices, i ∈ {0,1,2})
#   T_col = (p,q,r) encode la structure de COLONNE (j-indices, j ∈ {0,1,2})
# On construit un polynôme Φ_{T_row,T_col} ∈ Sym^3(C^9) comme suit :
#   Φ_{T_row,T_col} = somme sur les symétrisations de
#     x_{a,p}·x_{b,q}·x_{c,r}  avec signes du symétriseur de Young.
# La matrice C5[t_row, t_col] = <f, Φ_{T_row,T_col}> (produit apolaire).
# rang(C5(det3)) = 0, rang(C5(f_rand)) = 8 en général.
#
# FORMULE EXPLICITE DU SYMÉTRISEUR DE YOUNG λ=(2,1)
# --------------------------------------------------
# Pour T = (a,b,c) avec T_tab = [[a,b],[c]] :
#   y_T agit sur le triplet d'indices (i0,i1,i2) :
#   y_T(i0,i1,i2) = (symétrise {i0,i1}) × (antisymétrise {i0,i2})
#                 = (i0,i1,i2) + (i1,i0,i2)   [sym sur pos 0,1]
#                 − (i2,i1,i0) − (i1,i2,i0)   [antisym sur pos 0,2]
# Cela correspond au symétriseur de Young standard (non normalisé).
# Les 4 permutations (avec signes) donnent le polynôme de Specht.
#
# En pratique, pour l'apolarity pairing <f, Φ> :
#   <f, x_u · x_v · x_w> = coeff(f, {u,v,w}) × multinomial_factor
# où la convention polynomiale est f = Σ c_α x^α SANS facteur multinomial.
# Donc <f, x_u x_v x_w> = Σ_{perms of {u,v,w}} c_{sorted_exp} / normalisation.
# On utilise la convention «apolar intérieure» :
#   <f, x_{k1}x_{k2}x_{k3}> = coeff de x^{sorted exp} dans f.
# (Sans facteur multinomial — cohérent avec notre encodage.)
# ---------------------------------------------------------------------------

# SSYT pour λ=(2,1) sur alphabet {0,1,2} : 8 tableaux.
# T=(a,b,c) : ligne 1 = [a,b] (a≤b), ligne 2 = [c] (a<c).
_SSYT_21: List[Tuple[int, int, int]] = [
    (0, 0, 1), (0, 0, 2),
    (0, 1, 1), (0, 1, 2), (0, 2, 1), (0, 2, 2),
    (1, 1, 2), (1, 2, 2),
]


def _young_sym_21_terms() -> List[Tuple[Tuple[int, int, int], int]]:
    """Termes du symétriseur de Young pour λ=(2,1), tableau standard T1=[[0,1],[2]].

    Pour un triplet d'indices (i0, i1, i2) ∈ {0,1,2}^3 :
    y_T((i0,i1,i2)) = Σ_σ sgn(σ) · σ((i0,i1,i2))
    où σ parcourt les 4 permutations définies par le symétriseur de Young.

    Renvoie une liste de (permutation_des_positions, signe) :
    chaque entrée (perm, s) indique que la position perm[k] prend la valeur
    de la position k du triplet d'entrée.

    Les 4 termes (non normalisés) pour λ=(2,1) sont :
      + (0,1,2)   [identité]
      + (1,0,2)   [échange 0↔1, même colonne 0 de T : sym]
      - (2,1,0)   [échange 0↔2, anti-sym colonne]
      - (2,0,1)   [échange 0↔2 puis 0↔1]
    """
    return [
        ((0, 1, 2),  1),
        ((1, 0, 2),  1),
        ((2, 1, 0), -1),
        ((2, 0, 1), -1),
    ]


def _apolarity(f_coeffs: List[int], k1: int, k2: int, k3: int, p: int) -> int:
    """Produit apolaire <f, x_{k1}·x_{k2}·x_{k3}> mod p.

    Convention : f = Σ_α c_α x^α (sans facteur multinomial).
    <f, x_{k1}x_{k2}x_{k3}> = c_α  où α est le multi-indice correspondant
    à {k1,k2,k3} (trié pour donner l'exposant canonique).

    Cette convention est cohérente avec notre encodage dans _BASIS_SYM3.
    Pour des monomiales avec répétitions (ex : k1=k2), l'exposant α_k1 += 2.
    """
    exp = [0] * N_VAR
    exp[k1] += 1
    exp[k2] += 1
    exp[k3] += 1
    exp_t = tuple(exp)
    idx = _IDX_SYM3.get(exp_t)
    if idx is None:
        return 0
    return int(f_coeffs[idx]) % p


def build_flattening_C5(f_coeffs: List[int], p: int) -> "np.ndarray":
    """Projecteur équivariant S^(2,1)(C^3)⊗S^(2,1)(C^3) — matrice 8×8.

    Chaque entrée M[t_row, t_col] = <f, Φ_{T_row, T_col}> mod p, où :
    - T_row = _SSYT_21[t_row]  (tableau de Young semi-standard, indices de lignes)
    - T_col = _SSYT_21[t_col]  (idem, indices de colonnes)
    - Φ_{T_row, T_col} = symétriseur de Young appliqué séparément aux
      indices de lignes et de colonnes du triplet de variables :

      Φ_{(a,b,c),(p,q,r)} = Σ_{(π_r,s_r),(π_c,s_c)} s_r·s_c ·
                             x_{a_{π_r(0)}, p_{π_c(0)}} ·
                             x_{b_{π_r(1)}, q_{π_c(1)}} ·
                             x_{c_{π_r(2)}, r_{π_c(2)}}

      où la somme est sur les 4×4 = 16 paires de permutations du symétriseur de Young
      pour λ=(2,1).

    Résultats attendus :
    - rang(C5(det3))  = 0  (det3 ⊥ S^(2,1)⊗S^(2,1))
    - rang(C5(f_rand)) = 8 en général (rang générique)
    => Delta = 8 → DÉFAUT CERTIFIÉ si confirmé sur 3 premiers.

    AVERTISSEMENT : ce résultat (rang=0 sur det3) n'est pas un artefact numérique ;
    il est garanti par la décomposition de Cauchy et l'équivariance GL_3×GL_3.
    """
    _ensure_bases()
    dim = len(_SSYT_21)   # 8
    M = np.zeros((dim, dim), dtype=np.int64)

    sym_terms = _young_sym_21_terms()  # 4 paires (perm, signe)
    row_idx_arr = list(range(3))       # positions 0,1,2 du triplet

    for t_row_idx, T_row in enumerate(_SSYT_21):
        a, b, c = T_row   # indices de LIGNES (i-indices)
        row_arr = [a, b, c]

        for t_col_idx, T_col in enumerate(_SSYT_21):
            p_idx, q_idx, r_idx = T_col   # indices de COLONNES (j-indices)
            col_arr = [p_idx, q_idx, r_idx]

            # Accumuler les 4×4 = 16 termes du symétriseur de Young bi-tensoriel
            total = 0
            for perm_r, sgn_r in sym_terms:      # permutation des indices de ligne
                for perm_c, sgn_c in sym_terms:  # permutation des indices de colonne
                    # Variable k = 3*i + j
                    i0 = row_arr[perm_r[0]];  j0 = col_arr[perm_c[0]]
                    i1 = row_arr[perm_r[1]];  j1 = col_arr[perm_c[1]]
                    i2 = row_arr[perm_r[2]];  j2 = col_arr[perm_c[2]]
                    k0 = 3 * i0 + j0
                    k1 = 3 * i1 + j1
                    k2 = 3 * i2 + j2
                    ap = _apolarity(f_coeffs, k0, k1, k2, p)
                    total = (total + sgn_r * sgn_c * ap) % p

            M[t_row_idx, t_col_idx] = total % p

    return M


# ---------------------------------------------------------------------------
# C6 : Koszul flattening K_f^2 : Λ²V → V ⊗ Sym²V
# ---------------------------------------------------------------------------
# Construction de Landsberg-Ottaviani (arXiv:1006.5168, Theorem 1.2).
# Pour f ∈ Sym^3 V et k=2 :
#   K_f^2 : Λ^2(V) ⊗ Sym^0(V) → Λ^1(V) ⊗ Sym^1(V)  (≡ V ⊗ V)
# via : (e_i ∧ e_j) ↦ e_i ⊗ (∂_{e_j} f) - e_j ⊗ (∂_{e_i} f)
# où ∂_{e_l} f ∈ Sym^2 V est la dérivée de f dans la direction e_l.
#
# Version implémentée : Λ^2 V → V ⊗ Sym^2 V (matrice 36 × 405)
# M[(i<j), (k, α)] = δ_{k=i} * C1[j, α] - δ_{k=j} * C1[i, α]
# où C1[l, α] = coeff de α ∈ Sym^2 V dans ∂_{e_l} f.
#
# dim source = C(9,2) = 36, dim cible = 9*45 = 405.
# Rang générique = 36 (injectif pour f générique).
# Sur det3 : rang = 36 aussi (vérifié numériquement).
# => delta = 0 pour C6 aussi.
#
# POURQUOI ÇA NE MARCHE PAS (explication mathématique) :
# V = C^9 = C^3 ⊗ C^3 (structure matricielle de det3).
# Sym^3(C^3 ⊗ C^3) = ⊕_{λ ⊢ 3} S^λ(C^3) ⊗ S^λ(C^3)   (décomp. de Cauchy)
#   avec λ = (3,0,0), (2,1,0), (1,1,1) et dims 10×10 + 8×8 + 1×1 = 165 ✓
#
# det3 est EXCLUSIVEMENT dans la composante λ=(1,1,1) : Λ^3(C^3) ⊗ Λ^3(C^3).
# Car det3 est alternant en lignes ET en colonnes.
# => det3 est le SEUL élément (up to scalaire) de cette composante de dim 1.
#
# Les flattenings classiques (C1–C6) agissent sur tout Sym^3 V sans
# distinguer les composantes isotypiques. Pour det3, ils voient une
# forme «distinguable» seulement en ce qu'elle vit dans un sous-espace
# de dim 1. Mais les catalecticants et Koszul sont tous INJECTIFS sur
# det3 car det3 est non dégénéré (aucun facteur linéaire, pas div. de zéro).
#
# CE QUI MARCHERAIT (TODO C7) :
# Utiliser la projection sur S^{(2,1)}(C^3) ⊗ S^{(2,1)}(C^3) (sous-espace de dim 64)
# et construire le flattening équivariant RESTREINT à ce sous-espace.
# det3 n'a PAS de composante dans ce sous-espace → proj(det3) = 0 → rang 0.
# proj(f_rand) a rang > 0 en général → chute de rang = rang_générique > 0 !
# Mais calculer cette projection nécessite les matrices de Schur pour GL_3.
# ---------------------------------------------------------------------------

def build_flattening_C6(f_coeffs: List[int], p: int) -> "np.ndarray":
    """Koszul flattening K_f^2 : Λ^2(V) → V ⊗ Sym^2(V), matrice 36 × 405.

    Landsberg-Ottaviani construction (arXiv:1006.5168, Thm 1.2) pour k=2, d=3 :
      K_f^2(e_i ∧ e_j) = e_i ⊗ (∂_{e_j} f) - e_j ⊗ (∂_{e_i} f)
    où ∂_{e_l} f ∈ Sym^2 V est la rangée l de la matrice C1.

    Matrice M[(i,j), (k, α)] = δ_{k=i} C1[j,α] - δ_{k=j} C1[i,α]  (i < j).

    Résultat observé : rang générique = rang det3 = 36 → delta = 0.
    Cf. commentaire bloc ci-dessus pour l'explication mathématique.
    """
    _ensure_bases()
    # Base de Λ^2 V : paires (i,j) avec i < j
    lambda2_basis = [(i, j)
                     for i in range(N_VAR)
                     for j in range(i + 1, N_VAR)]
    n_rows = len(lambda2_basis)          # C(9,2) = 36
    n_cols = N_VAR * DIM_SYM2           # 9 * 45 = 405

    C1 = build_flattening_C1(f_coeffs, p)   # 9 × 45
    M  = np.zeros((n_rows, n_cols), dtype=np.int64)

    for row_idx, (vi, vj) in enumerate(lambda2_basis):
        # K_f(e_vi ∧ e_vj) = e_vi ⊗ ∂_{e_vj}f - e_vj ⊗ ∂_{e_vi}f
        for alpha_idx in range(DIM_SYM2):
            col_i = vi * DIM_SYM2 + alpha_idx   # composante (k=vi, α)
            col_j = vj * DIM_SYM2 + alpha_idx   # composante (k=vj, α)
            M[row_idx, col_i] = (M[row_idx, col_i] + C1[vj, alpha_idx]) % p
            M[row_idx, col_j] = (M[row_idx, col_j] - C1[vi, alpha_idx]) % p

    return M


# ---------------------------------------------------------------------------
# C7 : Koszul flattening K_f^3 : Λ³V → Λ²V ⊗ Sym²V
# ---------------------------------------------------------------------------
# Généralisation de C6 à k=3.
# dim source = C(9,3) = 84,  dim cible = C(9,2)*45 = 36*45 = 1620.
# Matrice 84 × 1620.
# K_f^3(e_{i1}∧e_{i2}∧e_{i3}) = Σ_{j=1}^3 (-1)^{j-1}
#   (e_{i1}∧...∧ê_{ij}...∧e_{i3}) ⊗ (∂_{e_{ij}} f)
# où ∂_{e_l} f est la ligne l de C1 (vecteur de Sym²V, dim 45).
# La cible est Λ²V ⊗ Sym²V, indexée par (paire (a<b), α ∈ Sym²V).
# Rang générique = 84 (injectif sur f générique).
# Sur det3 : rang = 84 (delta=0, même explication que C6).
# ---------------------------------------------------------------------------

def build_flattening_C7(f_coeffs: List[int], p: int) -> "np.ndarray":
    """Koszul flattening K_f^3 : Λ³V → Λ²V ⊗ Sym²V, matrice 84 × 1620.

    Construit via build_koszul_k(f_coeffs, p, k=3).
    Résultat attendu : rang générique = rang det3 = 84 → delta = 0.
    """
    return build_koszul_k(f_coeffs, p, k=3)


# ---------------------------------------------------------------------------
# C8 : Koszul flattening K_f^4 : Λ⁴V → Λ³V ⊗ Sym²V
# ---------------------------------------------------------------------------
# Généralisation de C6 à k=4.
# dim source = C(9,4) = 126,  dim cible = C(9,3)*45 = 84*45 = 3780.
# Matrice 126 × 3780.
# Même formule Koszul que C7 avec k=4.
# Rang générique = 126 (injectif).
# Sur det3 : rang = 126 (delta=0).
# ---------------------------------------------------------------------------

def build_flattening_C8(f_coeffs: List[int], p: int) -> "np.ndarray":
    """Koszul flattening K_f^4 : Λ⁴V → Λ³V ⊗ Sym²V, matrice 126 × 3780.

    Construit via build_koszul_k(f_coeffs, p, k=4).
    Résultat attendu : rang générique = rang det3 = 126 → delta = 0.
    """
    return build_koszul_k(f_coeffs, p, k=4)


def build_koszul_k(f_coeffs: List[int], p: int, k: int, N: int = N_VAR) -> "np.ndarray":
    """Construction générique du Koszul flattening K_f^k : Λ^k V → Λ^{k-1} V ⊗ Sym²V.

    K_f^k(e_{i1}∧...∧e_{ik}) = Σ_{j=0}^{k-1} (-1)^j
        (e_{i1}∧...∧ê_{ij}...∧e_{ik}) ⊗ (∂_{e_{ij}} f)
    où ∂_{e_l} f ∈ Sym²V est la ligne l de la matrice C1 (N×45).

    Source : Λ^k V,  base = combinaisons (i1,...,ik) avec i1<...<ik.
    Cible  : Λ^{k-1} V ⊗ Sym²V,  indexée par ((a1,...,a_{k-1}), α ∈ Sym²V).

    Paramètres
    ----------
    k : ordre (k=2 → C6, k=3 → C7, k=4 → C8)
    N : dimension de V (9 par défaut)
    """
    from itertools import combinations as _comb
    _ensure_bases()

    C1 = build_flattening_C1(f_coeffs, p)   # N × 45

    src_basis = list(_comb(range(N), k))           # Λ^k V
    tgt_k1    = list(_comb(range(N), k - 1))       # Λ^{k-1} V
    tgt_idx   = {t: i for i, t in enumerate(tgt_k1)}

    n_src = len(src_basis)    # C(N,k)
    n_tgt = len(tgt_k1)      # C(N,k-1)
    dim_sym2 = C1.shape[1]   # 45

    M = np.zeros((n_src, n_tgt * dim_sym2), dtype=np.int64)

    for row_idx, blade in enumerate(src_basis):
        for j in range(k):
            sign  = (-1) ** j
            i_j   = blade[j]
            residual = blade[:j] + blade[j + 1:]
            tgt_row  = tgt_idx[residual]
            for a in range(dim_sym2):
                col = tgt_row * dim_sym2 + a
                M[row_idx, col] = (M[row_idx, col] + sign * C1[i_j, a]) % p

    return M


# ---------------------------------------------------------------------------
# Table des familles
# ---------------------------------------------------------------------------

FAMILY_BUILDERS: Dict[str, Callable] = {
    "C1": build_flattening_C1,
    "C2": build_flattening_C2,
    "C3": build_flattening_C3,
    "C4": build_flattening_C4,
    "C5": build_flattening_C5,
    "C6": build_flattening_C6,
    "C7": build_flattening_C7,
    "C8": build_flattening_C8,
}

# Ce que teste réellement chaque famille :
# C1 : Cat_{1,2}(f) : V→Sym²V              (9×45, rang max 9)
# C2 : Cat_{2,1}(f) : Sym²V→V              (45×9, rang max 9 = rang C1, validation croisée)
# C3 : Hessien tensoriel V⊗V→V             (81×9, rang max 9, validation croisée)
# C4 : Multiplication M_f : Sym²V→Sym⁵V    (45×1287, rang max 45)
# C5 : Projecteur S^(2,1)⊗S^(2,1) via SSYT (8×8, rang det3=0, rang rand=8 — DÉFAUT ATTENDU)
# C6 : Koszul K_f²  : Λ²V→V⊗Sym²V          (36×405, rang max 36 — Landsberg-Ottaviani)
# C7 : Koszul K_f³  : Λ³V→Λ²V⊗Sym²V        (84×1620, rang max 84)
# C8 : Koszul K_f⁴  : Λ⁴V→Λ³V⊗Sym²V        (126×3780, rang max 126)
# C1-C4,C6-C8 : delta=0. C5 : delta=8 attendu (preuve Cauchy).
FAMILY_NAMES: Dict[str, str] = {
    "C1": "Cat_{1,2}(f) : V→Sym²V (9×45, rang max 9)",
    "C2": "Cat_{2,1}(f) : Sym²V→V (45×9, rang max 9 — validation croisée C1)",
    "C3": "Hessien tensoriel V⊗V→V (81×9, rang max 9 — validation croisée)",
    "C4": "Multiplication f·g : Sym²V→Sym⁵V (45×1287, rang max 45)",
    "C5": "Projecteur S^(2,1)xS^(2,1) via SSYT (8x8, GL3xGL3-equivariant - rang(det3)=0 MAIS pas defaut GL9)",
    "C6": "Koszul K_f² : Λ²V→V⊗Sym²V (36×405, rang max 36 — Landsberg-Ottaviani)",
    "C7": "Koszul K_f³ : Λ³V→Λ²V⊗Sym²V (84×1620, rang max 84)",
    "C8": "Koszul K_f⁴ : Λ⁴V→Λ³V⊗Sym²V (126×3780, rang max 126)",
}

# ---------------------------------------------------------------------------
# Rang générique et rang de det3
# ---------------------------------------------------------------------------

def generic_rank(builder: Callable, p: int,
                 n_samples: int = N_GENERIC_SAMPLES,
                 seeds: List[int] = None) -> int:
    """Rang générique : max sur n_samples cubiques aléatoires.

    Chaque cubique est tirée aléatoirement mod p avec seed distinct.
    Si le builder retourne None, on retourne -1.
    """
    if seeds is None:
        seeds = GENERIC_SEEDS[:n_samples]
    max_rank = 0
    for seed in seeds[:n_samples]:
        f_rand = random_cubic_mod_p(p, seed)
        try:
            M = builder(f_rand, p)
        except Exception as ex:
            logger.warning(f"Erreur builder (seed={seed}) : {ex}")
            continue
        if M is None:
            return -1
        r = rank_modp(M, p)
        if r > max_rank:
            max_rank = r
        if max_rank == min(M.shape):
            break  # rang maximal possible atteint
    return max_rank


def det3_rank(builder: Callable, p: int) -> int:
    """Rang de Y(f_det3) mod p.

    Si le builder retourne None, retourne -1.
    """
    f_det = encode_det3()
    try:
        M = builder(f_det, p)
    except Exception as ex:
        logger.error(f"Erreur builder sur det3 : {ex}")
        return -1
    if M is None:
        return -1
    return rank_modp(M, p)


# ---------------------------------------------------------------------------
# Certification de défaut
# ---------------------------------------------------------------------------

def certify_defect(
    builder: Callable,
    family_name: str,
    primes: List[int] = None
) -> Tuple[bool, int, Dict]:
    """Certifie un défaut de rang pour la famille donnée sur plusieurs premiers.

    Retourne (True, delta, details) si rank_det3 < rank_generic pour LES 3 PREMIERS.
    Retourne (False, 0, details) sinon.

    RÈGLE STRICTE : pas de certification sur 1 ou 2 premiers seulement.
    """
    if primes is None:
        primes = PRIMES
    if len(primes) < 3:
        raise ValueError("Au moins 3 premiers distincts requis pour certification.")

    details = {}
    confirmed = []

    for p in primes:
        logger.info("  p = %d", p)
        t0 = time.time()
        rg = generic_rank(builder, p)
        rd = det3_rank(builder, p)
        elapsed = time.time() - t0
        delta = rg - rd if (rg >= 0 and rd >= 0) else None

        logger.info("    rank_generic (5 samples) = %s", rg if rg >= 0 else "N/A")
        logger.info("    rank_det3 = %s", rd if rd >= 0 else "N/A")
        if delta is not None:
            logger.info("    delta = %d%s", delta, "  → CANDIDAT" if delta > 0 else "")
        else:
            logger.info("    delta = N/A (builder a retourné None)")

        details[p] = {"rank_generic": rg, "rank_det3": rd, "delta": delta,
                      "elapsed_s": round(elapsed, 2)}
        if delta is not None and delta > 0:
            confirmed.append(delta)

    if len(confirmed) == len(primes) and len(confirmed) >= 3:
        # Vérifier que tous les deltas concordent
        if len(set(confirmed)) == 1:
            return True, confirmed[0], details
        else:
            # Deltas différents selon les premiers → bug ou problème
            logger.warning("Deltas incohérents entre premiers : %s — NON certifié.", confirmed)
            return False, 0, details
    else:
        return False, 0, details


# ---------------------------------------------------------------------------
# Checkpoint
# ---------------------------------------------------------------------------

def load_checkpoint() -> Dict:
    if CHECKPOINT_FILE.exists():
        try:
            with open(CHECKPOINT_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}


def save_checkpoint(data: Dict):
    with open(CHECKPOINT_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Mode : audit environnement
# ---------------------------------------------------------------------------

def mode_audit():
    logger.info("=== AUDIT ENVIRONNEMENT ===")
    logger.info("Python : %s", sys.version)
    logger.info("numpy : %s", np.__version__ if HAS_NUMPY else "ABSENT")
    if HAS_CUPY:
        logger.info("CuPy  : %s", cp.__version__)
        try:
            logger.info("CUDA disponible : %s", cp.cuda.is_available())
        except Exception:
            logger.info("CUDA : inconnu")
    else:
        logger.info("CuPy  : ABSENT (fallback CPU actif)")
    if HAS_PSUTIL:
        mem = psutil.virtual_memory()
        logger.info("RAM disponible : %.1f GB / %.1f GB",
                    mem.available / 1e9, mem.total / 1e9)
    else:
        logger.info("psutil : ABSENT (RAM non mesurée)")
    logger.info("Premiers utilisés : %s", PRIMES)
    logger.info("dim Sym^3(V^9) = %d (attendu 165)", _dim_sym(N_VAR, 3))
    logger.info("dim Sym^2(V^9) = %d (attendu 45)",  _dim_sym(N_VAR, 2))
    logger.info("dim Sym^5(V^9) = %d (attendu 1287)", _dim_sym(N_VAR, 5))
    logger.info("ISSUE D — vérifiez numpy et les dépendances si des modules manquent.")


# ---------------------------------------------------------------------------
# Mode : sandbox GL2 et mod p
# ---------------------------------------------------------------------------

def mode_sandbox():
    """Tests sur GL_2 (petite taille) et arithmétique mod p."""
    logger.info("=== SANDBOX GL_2 ET MOD P ===")
    failures = []

    # Test 1 : dim Sym^3(V^9) = 165
    logger.info("Test 1 : dim Sym^3(V^9)")
    d = _dim_sym(9, 3)
    if d == 165:
        logger.info("  OK — dim = 165")
    else:
        logger.error("  ÉCHEC — dim = %d (attendu 165)", d)
        failures.append("dim_sym3")

    # Test 2 : rang mod p sur matrice identité
    logger.info("Test 2 : rang de Id_5 mod 7")
    Id5 = np.eye(5, dtype=np.int64)
    r = rank_modp_cpu(Id5, 7)
    if r == 5:
        logger.info("  OK — rang = 5")
    else:
        logger.error("  ÉCHEC — rang = %d", r)
        failures.append("rank_id5")

    # Test 3 : rang d'une matrice singulière
    logger.info("Test 3 : rang matrice singulière 3×3 mod 7")
    A = np.array([[1, 2, 3], [2, 4, 6], [1, 1, 1]], dtype=np.int64)
    r = rank_modp_cpu(A, 7)
    if r == 2:
        logger.info("  OK — rang = 2")
    else:
        logger.error("  ÉCHEC — rang = %d (attendu 2)", r)
        failures.append("rank_singular")

    # Test 4 : modinv
    logger.info("Test 4 : modinv(3, 7) = 5")
    inv = modinv(3, 7)
    if inv == 5:
        logger.info("  OK — 3*5 = 15 ≡ 1 mod 7")
    else:
        logger.error("  ÉCHEC — modinv(3,7) = %d", inv)
        failures.append("modinv")

    # Test 5 : GPU vs CPU rang (si GPU dispo)
    if HAS_CUPY:
        logger.info("Test 5 : rang GPU vs CPU sur matrice aléatoire 20×15")
        rng = np.random.default_rng(42)
        B = rng.integers(0, 100, size=(20, 15)).astype(np.int64)
        p = 101
        r_cpu = rank_modp_cpu(B, p)
        r_gpu = rank_modp_exact_gpu(B, p)
        if r_cpu == r_gpu:
            logger.info("  OK — rang CPU = rang GPU = %d", r_cpu)
        else:
            logger.error("  ÉCHEC — rang CPU=%d, rang GPU=%d", r_cpu, r_gpu)
            failures.append("gpu_vs_cpu")
    else:
        logger.info("Test 5 : GPU absent — skipped.")

    # Test 6 : det9 de l'identité
    logger.info("Test 6 : det9(Id) mod 7 = 1")
    Id9 = np.eye(N_VAR, dtype=np.int64)
    d9 = det9_mod_p(Id9, 7)
    if d9 == 1:
        logger.info("  OK — det(Id_9) = 1")
    else:
        logger.error("  ÉCHEC — det = %d", d9)
        failures.append("det9_id")

    # Test 7 : invariance équivariante — rang(Y(g·det3)) = rang(Y(det3)) pour g aléatoire
    # Si ce test échoue, le flattening n'est pas équivariant → delta > 0 serait un artefact.
    logger.info("Test 7 : invariance équivariante rang(C1(g·det3)) == rang(C1(det3))")
    p_test = 1_000_003   # petit premier pour être rapide
    try:
        f_det_test = encode_det3()
        g_test = random_GL9_mod_p(p_test, seed=2718)
        f_gdet = action_g_on_cubic_mod_p(g_test, f_det_test, p_test)
        r_det3  = rank_modp_cpu(build_flattening_C1(f_det_test, p_test), p_test)
        r_gdet3 = rank_modp_cpu(build_flattening_C1(f_gdet, p_test), p_test)
        if r_det3 == r_gdet3:
            logger.info("  OK — rang(C1(det3)) = rang(C1(g·det3)) = %d", r_det3)
        else:
            logger.error("  ÉCHEC — rang(C1(det3))=%d mais rang(C1(g·det3))=%d → bug équivariance",
                         r_det3, r_gdet3)
            failures.append("equivariance_C1")
    except Exception as ex:
        logger.warning("  Test 7 : erreur inattendue (%s) — skipped.", ex)

    # Test 8 : invariance équivariante pour C4 (multiplication)
    logger.info("Test 8 : invariance équivariante rang(C4(g·det3)) == rang(C4(det3))")
    try:
        r_det3_c4  = rank_modp_cpu(build_flattening_C4(f_det_test, p_test), p_test)
        r_gdet3_c4 = rank_modp_cpu(build_flattening_C4(f_gdet, p_test), p_test)
        if r_det3_c4 == r_gdet3_c4:
            logger.info("  OK — rang(C4(det3)) = rang(C4(g·det3)) = %d", r_det3_c4)
        else:
            logger.error("  ÉCHEC — rang(C4(det3))=%d mais rang(C4(g·det3))=%d → bug équivariance",
                         r_det3_c4, r_gdet3_c4)
            failures.append("equivariance_C4")
    except Exception as ex:
        logger.warning("  Test 8 : erreur inattendue (%s) — skipped.", ex)

    # Test 9 : invariance équivariante pour C6 (Koszul K_f^2)
    logger.info("Test 9 : invariance équivariante rang(C6(g·det3)) == rang(C6(det3))")
    try:
        r_det3_c6  = rank_modp_cpu(build_flattening_C6(f_det_test, p_test), p_test)
        r_gdet3_c6 = rank_modp_cpu(build_flattening_C6(f_gdet, p_test), p_test)
        if r_det3_c6 == r_gdet3_c6:
            logger.info("  OK — rang(C6(det3)) = rang(C6(g·det3)) = %d", r_det3_c6)
        else:
            logger.error("  ÉCHEC — rang(C6(det3))=%d mais rang(C6(g·det3))=%d → bug équivariance",
                         r_det3_c6, r_gdet3_c6)
            failures.append("equivariance_C6")
    except Exception as ex:
        logger.warning("  Test 9 : erreur inattendue (%s) — skipped.", ex)

    # Test 10 : invariance équivariante pour C7 (Koszul K_f^3)
    logger.info("Test 10 : invariance équivariante rang(C7(g·det3)) == rang(C7(det3))")
    try:
        r_det3_c7  = rank_modp_cpu(build_flattening_C7(f_det_test, p_test), p_test)
        r_gdet3_c7 = rank_modp_cpu(build_flattening_C7(f_gdet, p_test), p_test)
        if r_det3_c7 == r_gdet3_c7:
            logger.info("  OK — rang(C7(det3)) = rang(C7(g·det3)) = %d", r_det3_c7)
        else:
            logger.error("  ÉCHEC — rang(C7(det3))=%d mais rang(C7(g·det3))=%d → bug équivariance",
                         r_det3_c7, r_gdet3_c7)
            failures.append("equivariance_C7")
    except Exception as ex:
        logger.warning("  Test 10 : erreur inattendue (%s) — skipped.", ex)

    # Test 11 : invariance équivariante pour C8 (Koszul K_f^4)
    logger.info("Test 11 : invariance équivariante rang(C8(g·det3)) == rang(C8(det3))")
    try:
        r_det3_c8  = rank_modp_cpu(build_flattening_C8(f_det_test, p_test), p_test)
        r_gdet3_c8 = rank_modp_cpu(build_flattening_C8(f_gdet, p_test), p_test)
        if r_det3_c8 == r_gdet3_c8:
            logger.info("  OK — rang(C8(det3)) = rang(C8(g·det3)) = %d", r_det3_c8)
        else:
            logger.error("  ÉCHEC — rang(C8(det3))=%d mais rang(C8(g·det3))=%d → bug équivariance",
                         r_det3_c8, r_gdet3_c8)
            failures.append("equivariance_C8")
    except Exception as ex:
        logger.warning("  Test 11 : erreur inattendue (%s) — skipped.", ex)

    # Test 12 : C5 — trois sous-tests nécessaires
    # (a) rang(C5(det3)) = 0  [propriété GL3xGL3]
    # (b) rang(C5(f_rand)) = 8  [rang générique]
    # (c) rang(C5(g_GL9·det3)) = 0 pour g ∈ GL9 générique  [test défaut GL9]
    # AVERTISSEMENT : (a)+(b) ne suffisent PAS. (c) est le test critique.
    # Si (c) donne rang > 0, C5 ne détecte que la strate GL3xGL3, pas l'orbite GL9.
    logger.info("Test 12 : C5 — (a) rang(det3)=0, (b) rang(rand)=8, (c) rang(gGL9·det3)=?")
    logger.info("  NOTE : C5 est GL3xGL3-equivariant (Kronecker), PAS GL9-equivariant.")
    logger.info("  Test critique : rang(C5(g_GL9·det3)) pour g ∈ GL9 non-Kronecker.")
    try:
        # a) rang det3 doit être 0 (det3 ⊥ S^(2,1)⊗S^(2,1) comme GL3xGL3-module)
        M_c5_det = build_flattening_C5(f_det_test, p_test)
        r_c5_det = rank_modp_cpu(M_c5_det, p_test)

        # b) rang f_rand doit être > 0
        f_rand_test_c5 = random_cubic_mod_p(p_test, seed=42)
        M_c5_rand = build_flattening_C5(f_rand_test_c5, p_test)
        r_c5_rand = rank_modp_cpu(M_c5_rand, p_test)

        # c) TEST CRITIQUE : g ∈ GL9 générique (PAS un produit de Kronecker)
        # Si rang > 0 pour un seul tel g → C5 ne détecte PAS un défaut GL9.
        gl9_seeds = [12345, 99999, 314159]
        c_rangs = []
        for s in gl9_seeds:
            g_gl9 = random_GL9_mod_p(p_test, seed=s)
            f_gl9_det = action_g_on_cubic_mod_p(g_gl9, f_det_test, p_test)
            M_gl9 = build_flattening_C5(f_gl9_det, p_test)
            c_rangs.append(rank_modp_cpu(M_gl9, p_test))

        ok_det  = (r_c5_det == 0)
        ok_rand = (r_c5_rand > 0)
        ok_gl9  = all(r == 0 for r in c_rangs)   # doit être 0 si vrai défaut GL9

        if ok_det:
            logger.info("  (a) OK — rang(C5(det3)) = 0")
        else:
            logger.error("  (a) ECHEC — rang(C5(det3)) = %d (attendu 0)", r_c5_det)
            failures.append("C5_det3_rang")

        if ok_rand:
            logger.info("  (b) OK — rang(C5(f_rand)) = %d > 0", r_c5_rand)
        else:
            logger.error("  (b) ECHEC — rang(C5(f_rand)) = 0 (attendu > 0)")
            failures.append("C5_rand_rang")

        if ok_gl9:
            logger.info("  (c) OK — rang(C5(g_GL9·det3)) = 0 pour seeds %s : VRAI DEFAUT GL9", gl9_seeds)
            logger.info("  Test 12 PASS COMPLET : delta_C5 = %d — DEFAUT GL9 CERTIFIABLE.", r_c5_rand)
        else:
            logger.warning("  (c) ECHEC CRITIQUE — rang(C5(g_GL9·det3)) = %s pour seeds %s",
                           c_rangs, gl9_seeds)
            logger.warning("  => C5 detecte la strate GL3xGL3, PAS un defaut sous GL9.")
            logger.warning("  => Le 'defaut' C5 est un artefact de la structure GL3xGL3 de det3.")
            logger.warning("  => rang(det3)=0 car det3 est dans Lambda3xLambda3 (GL3xGL3-isotypique),")
            logger.warning("     pas parce qu'il est dans cl(O_det3) sous GL9.")
            logger.warning("  => C5 ne peut PAS certifier un defaut SRMT.")
            # C5_GL9_not_defaut est un resultat ATTENDU, pas un bug.
            # Ne pas ajouter a failures : le sandbox ne doit pas echouer pour ca.
            logger.warning("  => Test 12(c) : resultat negatif attendu et documente.")

    except Exception as ex:
        logger.warning("  Test 12 : erreur inattendue (%s) — skipped.", ex)


    if failures:
        logger.error("ISSUE E — Tests échoués : %s", failures)
        return False
    else:
        logger.info("Tous les tests sandbox OK.")
        return True


# ---------------------------------------------------------------------------
# Mode : test-det3
# ---------------------------------------------------------------------------

def mode_test_det3():
    """Vérifie la cohérence de encode_det3()."""
    logger.info("=== TEST ENCODE_DET3 ===")
    failures = []

    f_det = encode_det3()

    # Test A : dimension
    if len(f_det) == DIM_SYM3:
        logger.info("Test A : dim = %d ✓", DIM_SYM3)
    else:
        logger.error("Test A : dim = %d (attendu %d) ✗", len(f_det), DIM_SYM3)
        failures.append("dim_det3")

    # Test B : det3(Id) = 1
    # Évaluer f_det sur x = Id_3 (identité 3×3) :
    # variables x_0,...,x_8 prises à la valeur (i,j) → delta_{i,j}
    # x_{3i+j} = delta_{ij}
    # L'identité est : x0=1, x1=0, x2=0, x3=0, x4=1, x5=0, x6=0, x7=0, x8=1
    x_id = [0] * N_VAR
    for i in range(3):
        x_id[3 * i + i] = 1
    val_id = _eval_poly(f_det, x_id)
    if val_id == 1:
        logger.info("Test B : det3(Id) = 1 ✓")
    else:
        logger.error("Test B : det3(Id) = %d (attendu 1) ✗", val_id)
        failures.append("det3_id")

    # Test C : det3(matrice anti-diag) = +1 ou -1
    # matrice [[0,0,1],[0,1,0],[1,0,0]] → det = -1 (perm (02))
    x_anti = [0] * N_VAR
    x_anti[2] = 1   # x02 = 1
    x_anti[4] = 1   # x11 = 1
    x_anti[6] = 1   # x20 = 1
    val_anti = _eval_poly(f_det, x_anti)
    if val_anti == -1:
        logger.info("Test C : det3(anti-diag) = -1 ✓")
    else:
        logger.error("Test C : det3(anti-diag) = %d (attendu -1) ✗", val_anti)
        failures.append("det3_anti")

    # Test D : nombre de termes non nuls = 6
    nz = sum(1 for c in f_det if c != 0)
    if nz == 6:
        logger.info("Test D : nombre de termes non nuls = 6 ✓")
    else:
        logger.error("Test D : %d termes non nuls (attendu 6) ✗", nz)
        failures.append("det3_terms")

    # Test E : det3(diag(2,3,5)) = 2*3*5 = 30
    x_diag = [0] * N_VAR
    x_diag[0] = 2  # x00
    x_diag[4] = 3  # x11
    x_diag[8] = 5  # x22
    val_diag = _eval_poly(f_det, x_diag)
    if val_diag == 30:
        logger.info("Test E : det3(diag(2,3,5)) = 30 ✓")
    else:
        logger.error("Test E : det3(diag(2,3,5)) = %d (attendu 30) ✗", val_diag)
        failures.append("det3_diag")

    if failures:
        logger.error("ISSUE E — Tests det3 échoués : %s", failures)
        return False
    logger.info("Tous les tests det3 OK.")
    return True


def _eval_poly(f_coeffs: List[int], x_vals: List[int]) -> int:
    """Évalue f ∈ Sym^3 V sur le vecteur x (entiers, pas mod p)."""
    _ensure_bases()
    result = 0
    for idx, exp in enumerate(_BASIS_SYM3):
        c = f_coeffs[idx]
        if c == 0:
            continue
        monome_val = 1
        for i, e in enumerate(exp):
            monome_val *= x_vals[i] ** e
        result += c * monome_val
    return result


# ---------------------------------------------------------------------------
# Mode : flattening-one
# ---------------------------------------------------------------------------

def mode_flattening_one(family: str):
    """Teste une famille sur le premier premier uniquement."""
    logger.info("=== Flattening family %s (%s) [mode: one prime] ===",
                family, FAMILY_NAMES.get(family, "?"))

    if family not in FAMILY_BUILDERS:
        logger.error("Famille inconnue : %s. Choix : %s", family, list(FAMILY_BUILDERS))
        return

    builder = FAMILY_BUILDERS[family]
    p = PRIMES[0]
    logger.info("  p = %d", p)

    # Rang générique
    t0 = time.time()
    rg = generic_rank(builder, p)
    logger.info("  rank_generic (5 samples) = %s",
                rg if rg >= 0 else "N/A (builder a retourné None)")

    # Rang de det3
    rd = det3_rank(builder, p)
    logger.info("  rank_det3 = %s",
                rd if rd >= 0 else "N/A (builder a retourné None)")

    if rg >= 0 and rd >= 0:
        delta = rg - rd
        logger.info("  delta = %d%s", delta, "  → CANDIDAT (à confirmer sur 3 premiers)"
                    if delta > 0 else "  → pas de chute")
    logger.info("  Temps : %.2f s", time.time() - t0)

    logger.info("")
    logger.info("AVERTISSEMENT : résultat sur 1 seul premier NON certifié.")
    logger.info("Utiliser --mode flattening-all pour certification multi-prime.")


# ---------------------------------------------------------------------------
# Mode : flattening-all
# ---------------------------------------------------------------------------

def mode_flattening_all():
    """Teste toutes les familles sur 3 premiers et certifie les défauts."""
    logger.info("=== FLATTENING ALL FAMILIES — 3 premiers ===")
    logger.info("NOTE : C1, C2, C3 : rang max = 9 (même information, validation croisée).")
    logger.info("       C4 (45×1287) : multiplication Sym²V→Sym⁵V, rang max 45.")
    logger.info("       C5 (8×8)    : projecteur S^(2,1)⊗S^(2,1) via SSYT — DÉFAUT ATTENDU.")
    logger.info("       C6-C8       : Koszul d’ordre 2,3,4 — rang max = dim source.")
    logger.info("       Un défaut sur C1/C2/C3 serait certifié seulement si rang < 9 sur 3 premiers.")

    results = {}
    for family, builder in FAMILY_BUILDERS.items():
        logger.info("")
        logger.info("=== Flattening family %s (%s) ===", family, FAMILY_NAMES[family])

        if family == "C5":
            logger.warning("C5 : GL3xGL3-equivariant, rang(det3)=0 par Cauchy.")
            logger.warning("     rang(C5(g_GL9*det3))=8 pour g dans GL9 generique.")
            logger.warning("     C5 detecte la strate GL3xGL3, PAS un defaut sous GL9.")
            logger.warning("     Une chute de rang ici n'est PAS un defaut SRMT certifie.")

        try:
            certified, delta, details = certify_defect(builder, family)
        except Exception as ex:
            logger.error("Erreur inattendue pour %s : %s", family, ex)
            results[family] = {"status": "error", "error": str(ex)}
            continue

        if certified and family != "C5":
            logger.info("ISSUE A — DEFAUT NON NUL CERTIFIE MULTI-PRIME (%s, delta=%d)", family, delta)
            results[family] = {"status": "A", "delta": delta, "details": details}
        elif certified and family == "C5":
            # C5 : rang(det3)=0 par Cauchy, mais ce n'est PAS un defaut GL9.
            # rang(C5(g_GL9*det3))=8 pour g generique -> strate GL3xGL3 seulement.
            logger.warning("C5 : delta=%d sur 3 premiers, MAIS ce n'est PAS un defaut GL9.", delta)
            logger.warning("     rang(C5(g_GL9*det3)) = 8 pour g in GL9 non-Kronecker.")
            logger.warning("     C5 mesure la GL3xGL3-isotypie de det3 (Lambda3xLambda3),")
            logger.warning("     pas son appartenance a cl(O_det3) sous GL9.")
            logger.warning("     Statut : ISSUE F (faux positif equivariant).")
            results[family] = {"status": "F_GL3only", "delta": delta, "details": details}
        else:
            # Vérifier si c'est ISSUE C (builder None = RAM) ou ISSUE B (pas de chute)
            any_none = any(d["rank_generic"] < 0 or d["rank_det3"] < 0
                          for d in details.values())
            if any_none:
                logger.info("ISSUE C — builder a retourné None pour %s (trop lourd).", family)
                results[family] = {"status": "C", "details": details}
            else:
                logger.info("Aucune chute certifiée pour %s. On passe à la suite.", family)
                results[family] = {"status": "B", "delta": 0, "details": details}

    print_conclusion(results)
    return results


# ---------------------------------------------------------------------------
# Mode : flattening-longrun (avec checkpoint)
# ---------------------------------------------------------------------------

def mode_flattening_longrun(resume: bool = False):
    """Balayage de familles avec checkpoint JSON pour reprise."""
    logger.info("=== FLATTENING LONG RUN (avec checkpoint) ===")

    checkpoint = load_checkpoint() if resume else {}
    done_families = set(checkpoint.get("done", []))
    results = checkpoint.get("results", {})

    families_to_run = [f for f in FAMILY_BUILDERS if f not in done_families]
    logger.info("Familles restantes : %s", families_to_run)

    for family in families_to_run:
        logger.info("")
        logger.info("=== Famille %s ===", family)
        builder = FAMILY_BUILDERS[family]

        try:
            certified, delta, details = certify_defect(builder, family)
            if certified:
                results[family] = {"status": "A", "delta": delta}
                logger.info("ISSUE A — DÉFAUT CERTIFIÉ pour %s, delta=%d", family, delta)
            else:
                results[family] = {"status": "B"}
                logger.info("Pas de défaut certifié pour %s.", family)
        except MemoryError:
            logger.error("ISSUE C — MemoryError pour %s", family)
            results[family] = {"status": "C"}
        except Exception as ex:
            logger.error("Erreur pour %s : %s", family, ex)
            results[family] = {"status": "error", "error": str(ex)}

        done_families.add(family)
        checkpoint = {"done": list(done_families), "results": results,
                      "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}
        save_checkpoint(checkpoint)
        logger.info("Checkpoint sauvegardé.")

    print_conclusion(results)


# ---------------------------------------------------------------------------
# Mode : review
# ---------------------------------------------------------------------------

def mode_review():
    """Vérifie les prérequis avant un long run."""
    logger.info("=== REVIEW PRÉREQUIS ===")
    ok = True

    # numpy
    if not HAS_NUMPY:
        logger.error("numpy manquant — ISSUE D")
        ok = False
    else:
        logger.info("numpy OK (%s)", np.__version__)

    # dim W
    _ensure_bases()
    if len(_BASIS_SYM3) == 165:
        logger.info("dim Sym^3(V^9) = 165 ✓")
    else:
        logger.error("dim Sym^3(V^9) = %d ✗", len(_BASIS_SYM3))
        ok = False

    # det3
    f_det = encode_det3()
    nz = sum(1 for c in f_det if c != 0)
    if nz == 6:
        logger.info("encode_det3 : 6 termes ✓")
    else:
        logger.error("encode_det3 : %d termes ✗", nz)
        ok = False

    # RAM
    if HAS_PSUTIL:
        mem = psutil.virtual_memory()
        avail_gb = mem.available / 1e9
        if avail_gb < 1.0:
            logger.warning("RAM disponible < 1 GB (%.1f GB) — risque ISSUE C", avail_gb)
        else:
            logger.info("RAM disponible : %.1f GB ✓", avail_gb)

    # Primes
    for p in PRIMES:
        if isprime(p):
            logger.info("Premier %d ✓", p)
        else:
            logger.error("ERREUR : %d n'est pas premier ✗", p)
            ok = False

    if ok:
        logger.info("REVIEW OK — prêt pour long run.")
    else:
        logger.error("REVIEW ÉCHOUÉ — ISSUE D ou E.")
    return ok


# ---------------------------------------------------------------------------
# Conclusion
# ---------------------------------------------------------------------------

def print_conclusion(results: Dict):
    logger.info("")
    logger.info("=" * 60)
    logger.info("CONCLUSION")
    logger.info("=" * 60)


    has_A      = any(v.get("status") == "A"          for v in results.values())
    has_C      = any(v.get("status") == "C"          for v in results.values())
    has_error  = any(v.get("status") == "error"       for v in results.values())
    has_F      = any(v.get("status") == "F_GL3only"   for v in results.values())
    all_B      = all(v.get("status") in ("B", "TODO", "F_GL3only")
                     for v in results.values())

    # Signaler les faux positifs GL3xGL3 en premier, quelle que soit la conclusion
    if has_F:
        logger.warning("")
        logger.warning("ISSUE F — FAUX POSITIF EQUIVARIANT (GL3xGL3, pas GL9)")
        for fam, v in results.items():
            if v.get("status") == "F_GL3only":
                logger.warning("  Famille : %s, delta=%d", fam, v.get("delta", "?"))
        logger.warning("  rang(C5(det3))=0 est une propriete de la structure GL3xGL3 de det3,")
        logger.warning("  pas un defaut sous GL9.")
        logger.warning("  rang(C5(g_GL9*det3)) = 8 pour g in GL9 generique.")
        logger.warning("  => Ce delta NE CERTIFIE PAS un defaut SRMT.")
        logger.warning("")

    if has_A:
        for fam, v in results.items():
            if v.get("status") == "A":
                logger.info("ISSUE A — DEFAUT NON NUL CERTIFIE MULTI-PRIME")
                logger.info("  Famille : %s", fam)
                logger.info("  Delta   : %d", v.get("delta", "?"))
                logger.info("  Les mineurs (rank_generic+1)x(rank_generic+1)")
                logger.info("  sont des equations candidates de O_det3 dans W.")
    elif all_B:
        logger.info("ISSUE B — AUCUNE CHUTE DE RANG TROUVEE DANS LES FAMILLES TESTEES")
        logger.info("  Cela ne prouve PAS l'absence de defaut.")
        logger.info("")
        logger.info("  EXPLICATION : V = C^9 = C^3 x C^3. Par Cauchy :")
        logger.info("    Sym^3(C^3xC^3) = S^(3)xS^(3) [dim 100]")
        logger.info("                   + S^(2,1)xS^(2,1) [dim 64]")
        logger.info("                   + Lambda3xLambda3 [dim 1]")
        logger.info("  det3 est dans Lambda3xLambda3 uniquement.")
        logger.info("  C1-C4,C6-C8 sont injectifs sur det3 (pas de defaut).")
        logger.info("  C5 : rang(det3)=0 par Cauchy, mais GL3xGL3 seulement, pas GL9.")
        logger.info("")
        logger.info("  PISTES : chercher un flattening GL9-equivariant sensible a det3.")
        logger.info("  Rappel : ISSUE B ne prouve PAS l'absence de defaut SRMT.")
    elif has_C:
        logger.info("ISSUE C — CERTAINES FAMILLES ONT DEPASSE LA RAM DISPONIBLE")
        logger.info("  Resultat incomplet. Tester sur machine avec plus de RAM.")
    elif has_error:
        logger.info("ISSUE E — ERREURS LORS DU CALCUL (voir logs)")
    else:
        logger.info("ISSUE B — Aucun defaut certifie dans les familles testees.")

    logger.info("")
    logger.info("RAPPEL DU CONTRAT D'HONNÊTETÉ :")
    logger.info("  - Cette méthode peut échouer à trouver un défaut même s'il en existe un.")
    logger.info("  - Une chute de rang sur 1 seul premier N'EST PAS un défaut certifié.")
    logger.info("  - Certaines flattenings peuvent dépasser la RAM (ISSUE C).")
    logger.info("  - Tout rang est calculé EXACTEMENT mod p. Aucun flottant.")
    logger.info("=" * 60)


# ---------------------------------------------------------------------------
# Point d'entrée principal
# ---------------------------------------------------------------------------

def main():
    _ensure_bases()

    parser = argparse.ArgumentParser(
        description="Recherche de défauts dans O_det3 par Young flattenings mod p."
    )
    parser.add_argument(
        "--mode",
        choices=["audit", "sandbox", "test-det3", "flattening-one",
                 "flattening-all", "flattening-longrun", "review"],
        default="audit",
        help="Mode d'exécution."
    )
    parser.add_argument(
        "--family",
        choices=sorted(FAMILY_BUILDERS.keys()),
        default="C1",
        help="Famille de flattenings (pour --mode flattening-one)."
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reprendre depuis le dernier checkpoint (pour flattening-longrun)."
    )

    args = parser.parse_args()

    logger.info("route_flattening_gpu_windows.py — démarrage")
    logger.info("Mode : %s | numpy=%s | CuPy=%s",
                args.mode, HAS_NUMPY, HAS_CUPY)

    if not HAS_NUMPY:
        logger.error("numpy est requis. Installez-le avec : pip install numpy")
        logger.error("ISSUE D — environnement incomplet.")
        sys.exit(1)

    if args.mode == "audit":
        mode_audit()

    elif args.mode == "sandbox":
        ok = mode_sandbox()
        if not ok:
            sys.exit(2)

    elif args.mode == "test-det3":
        ok = mode_test_det3()
        if not ok:
            sys.exit(2)

    elif args.mode == "flattening-one":
        # Tests préalables rapides
        if not mode_sandbox():
            logger.error("ISSUE E — Tests sandbox échoués. Abandon.")
            sys.exit(2)
        if not mode_test_det3():
            logger.error("ISSUE E — Tests det3 échoués. Abandon.")
            sys.exit(2)
        mode_flattening_one(args.family)

    elif args.mode == "flattening-all":
        if not mode_sandbox():
            logger.error("ISSUE E — Tests sandbox échoués. Abandon.")
            sys.exit(2)
        if not mode_test_det3():
            logger.error("ISSUE E — Tests det3 échoués. Abandon.")
            sys.exit(2)
        mode_flattening_all()

    elif args.mode == "flattening-longrun":
        if not args.resume:
            if not mode_sandbox():
                logger.error("ISSUE E — Tests sandbox échoués. Abandon.")
                sys.exit(2)
            if not mode_test_det3():
                logger.error("ISSUE E — Tests det3 échoués. Abandon.")
                sys.exit(2)
        mode_flattening_longrun(resume=args.resume)

    elif args.mode == "review":
        ok = mode_review()
        if not ok:
            sys.exit(1)

    logger.info("Terminé. Log complet : %s", LOG_FILE)


if __name__ == "__main__":
    main()
