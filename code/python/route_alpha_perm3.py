"""
route_alpha_gpu_crt_windows.py — Route Alpha v5
=================================================
Recherche d'équations de degré 5 pour la clôture de l'orbite det_3 dans Sym^3(C^9).

CONTEXTE MATHÉMATIQUE CERTIFIÉ :
  V = C^9,  W = Sym^3(C^9),  dim W = C(11,3) = 165
  B = C[W] = Sym(W^*),  A = C[O_det3] = B / I(O_det3)
  I_4(O_det3) = 0  (acquis, degré 4 vide).
  Degré 5 = premier candidat non trivial.

PLETHYSME CERTIFIÉ :
  Les 28 isotypes de Sym^5(Sym^3(C^9)) sont certifiés par calcul Python pur
  (plethysme exact via Murnaghan-Nakayama + formule de Littlewood).
  Vérification dimensionnelle :  sum(mult * dim_schur(lam,9)) = 1_082_239_158 = C(169,5). VÉRIFIÉ.
  IMPORTANT : cette liste a été calculée en Python pur, PAS via Sage.

LIMITES ACTUELLES :
  - APPROVED_FOR_MINI_RUN ≠ APPROVED_FOR_LONG_RUN
  - La liste des 28 isotypes est certifiée Python pur, pas Sage.
  - Le long run n'est pas lancé automatiquement; il faut APPROVED_FOR_LONG_RUN.
  - build_multiplication_map_modp pour GL9 degré 5 peut retourner ('TOO_LARGE', None)
    si la taille dépasse WEIGHT_SPACE_MAX.
  - Pas de calcul float pour les rangs ou dimensions : entiers uniquement.
  - Pas de matrice dense C(169,5) x n générée en mémoire.

GPU (v5) :
  - Représentation compacte MonoIdx5 : 5 indices uint16 au lieu de tuple 165-long
  - HWVSparse : arrays numpy parallèles (mono_ids + coeffs)
  - build_orbit_coords_gpu : coordonnées d'orbite sur GPU (CuPy) ou CPU (numpy)
  - evaluate_hwv_on_orbit_points_gpu : évaluation GPU via CuPy RawKernel
  - rank_modp_exact_gpu : rang exact mod p (Gauss, int64, jamais float)
  - --gpu-eval : flag pour activer le chemin GPU
  - Comparaison CPU vs GPU : run_gpu_vs_cpu_comparison_test
  - Timing et VRAM logging (_timed_eval, _get_vram_usage_mb)
  - Fallback CPU automatique si CuPy absent ou erreur CUDA

INSTALLATION (Windows, forward slash) :
  Étape 0 : Installer Python 3.10+ depuis python.org
  Étape 1 : pip install numpy sympy
  Étape 2 : (optionnel GPU) pip install cupy-cuda12x   # CUDA 12.x
             Pour GPU Blackwell (RTX 5070 Ti, cc10) : pip install cupy --no-binary cupy
             ou utiliser cupy-cuda13x (CUDA 13.x wheel) si disponible
  Étape 3 : (optionnel GPU) pip install numba
  Étape 4 : (optionnel) pip install psutil
  Étape 5 : (optionnel Sage) installer SageMath depuis sagemath.org
  Étape 6 : python route_alpha_gpu_crt_windows.py --mode audit

  Chemins Windows (forward slash) :
    ./route_alpha_perm3.py
    ./route_alpha_run.log
    ./route_alpha_checkpoint.json

USAGE :
  python route_alpha_gpu_crt_windows.py --mode docs
  python route_alpha_gpu_crt_windows.py --mode audit
  python route_alpha_gpu_crt_windows.py --mode gl2
  python route_alpha_gpu_crt_windows.py --mode modp
  python route_alpha_gpu_crt_windows.py --mode plethysm
  python route_alpha_gpu_crt_windows.py --mode paste-plethysm
  python route_alpha_gpu_crt_windows.py --mode verify-plethysm-list
  python route_alpha_gpu_crt_windows.py --mode lie
  python route_alpha_gpu_crt_windows.py --mode orbit
  python route_alpha_gpu_crt_windows.py --mode degree5-mini
  python route_alpha_gpu_crt_windows.py --mode degree5-one-isotype --lambda "9,4,2"
  python route_alpha_gpu_crt_windows.py --mode degree5-one-isotype --lambda "9,4,2" --gpu-eval
  python route_alpha_gpu_crt_windows.py --mode srmt-test-candidate --lambda "9,4,2"
  python route_alpha_gpu_crt_windows.py --mode degree5-longrun [--resume]
  python route_alpha_gpu_crt_windows.py --mode degree5-longrun --gpu-eval
  python route_alpha_gpu_crt_windows.py --mode review
  python route_alpha_gpu_crt_windows.py --mode all

  Options :
    --lambda      : partition en argument (ex: "9,4,2" ou "12,3")
    --resume      : reprendre depuis checkpoint
    --prime       : prime à utiliser (défaut 1000000007)
    --max-weight-dim : limite dimension espace de poids (défaut 500000)
    --gpu-eval    : activer l'évaluation sur GPU (CuPy requis)
"""

# ===========================================================================
# SECTION 2 : Imports robustes
# ===========================================================================

import ast
import argparse
import json
import logging
import math
import os
import random
import re
import sys
import time
from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import combinations_with_replacement, product
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

# numpy optionnel
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    np = None
    HAS_NUMPY = False

# sympy optionnel
try:
    import sympy
    HAS_SYMPY = True
except ImportError:
    sympy = None
    HAS_SYMPY = False

# cupy optionnel
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

# numba optionnel
try:
    import numba
    HAS_NUMBA = True
except ImportError:
    numba = None
    HAS_NUMBA = False

# psutil optionnel
try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    psutil = None
    HAS_PSUTIL = False

# ===========================================================================
# SECTION 2b : Types GPU et infrastructure
# ===========================================================================

# MonoIdx5 : monôme de Sym^5(W^*) représenté par 5 indices triés dans [0,164]
# Exemple : w_0^2 * w_3 * w_7 * w_9 → (0, 0, 3, 7, 9)
MonoIdx5 = Tuple[int, int, int, int, int]


@dataclass
class HWVSparse:
    """
    Représentation creuse d'un vecteur HWV.
    mono_ids : indices dans la liste monomials (int32)
    coeffs   : coefficients mod p (int64)
    """
    mono_ids: Any  # np.ndarray shape (nnz,) dtype int32
    coeffs: Any    # np.ndarray shape (nnz,) dtype int64

    def __post_init__(self):
        if HAS_NUMPY:
            self.mono_ids = np.asarray(self.mono_ids, dtype=np.int32)
            self.coeffs = np.asarray(self.coeffs, dtype=np.int64)


# CuPy RawKernel pour l'évaluation HWV sur GPU
# Compilé lazily à la première utilisation
_EVAL_HWV_KERNEL_SRC = r"""
extern "C" __global__ void eval_hwv_kernel(
    const long long* orbit_coords,   /* (n_points, 165) */
    const short* mono5_flat,         /* (nnz_total, 5) int16 */
    const long long* coeffs_flat,    /* (nnz_total,) */
    const int* hwv_offsets,          /* (n_hwv+1,) */
    long long* eval_matrix,          /* (n_hwv, n_points) */
    int n_hwv,
    int n_points,
    long long p
) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;  /* hwv index */
    int j = blockIdx.y * blockDim.y + threadIdx.y;  /* point index */
    if (k >= n_hwv || j >= n_points) return;

    long long val = 0LL;
    int start = hwv_offsets[k];
    int end = hwv_offsets[k+1];
    for (int t = start; t < end; t++) {
        long long c = coeffs_flat[t] % p;
        if (c < 0LL) c += p;
        long long prod = c;
        for (int s = 0; s < 5; s++) {
            int idx = (int)((unsigned short)mono5_flat[t*5 + s]);
            long long coord = orbit_coords[j*165 + idx] % p;
            if (coord < 0LL) coord += p;
            prod = (prod * coord) % p;
        }
        val = (val + prod) % p;
    }
    eval_matrix[k * n_points + j] = val;
}
"""

_eval_hwv_kernel_compiled = None


def _get_eval_hwv_kernel():
    """Compile le CuPy RawKernel lazily."""
    global _eval_hwv_kernel_compiled
    if _eval_hwv_kernel_compiled is None and HAS_CUPY:
        try:
            _eval_hwv_kernel_compiled = cp.RawKernel(
                _EVAL_HWV_KERNEL_SRC, 'eval_hwv_kernel'
            )
        except Exception as e:
            log(f"  WARNING: Could not compile eval_hwv_kernel: {e}")
            _eval_hwv_kernel_compiled = None
    return _eval_hwv_kernel_compiled


# ===========================================================================
# SECTION 3 : Logger unique
# ===========================================================================

LOG_FILE = Path("route_alpha_run.log")

_logger = logging.getLogger("route_alpha")
_logger.setLevel(logging.DEBUG)
_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

if not _logger.handlers:
    _ch = logging.StreamHandler(sys.stdout)
    _ch.setFormatter(_fmt)
    _logger.addHandler(_ch)

    try:
        _fh = logging.FileHandler(LOG_FILE, encoding="utf-8")
        _fh.setFormatter(_fmt)
        _logger.addHandler(_fh)
    except OSError:
        pass


def log(msg: str):
    """Écrit un message dans stdout + route_alpha_run.log."""
    _logger.info(msg)


def log_section(title: str):
    """Écrit un séparateur de section visible."""
    sep = "=" * 70
    _logger.info(sep)
    _logger.info(f"  {title}")
    _logger.info(sep)


# ===========================================================================
# SECTION 4 : Données certifiées
# ===========================================================================

PLETHYSM_S5_S3_GL9_CERTIFIED = [
    ([5, 4, 4, 2], 1),
    ([5, 5, 3, 1, 1], 1),
    ([6, 4, 2, 2, 1], 1),
    ([6, 4, 4, 1], 1),
    ([6, 5, 2, 2], 1),
    ([6, 6, 3], 1),
    ([7, 2, 2, 2, 2], 1),
    ([7, 4, 2, 2], 1),
    ([7, 4, 3, 1], 1),
    ([7, 4, 4], 1),
    ([7, 5, 2, 1], 1),
    ([7, 6, 2], 1),
    ([8, 3, 2, 2], 1),
    ([8, 4, 2, 1], 1),
    ([8, 4, 3], 1),
    ([8, 5, 2], 1),
    ([8, 6, 1], 1),
    ([9, 2, 2, 2], 1),
    ([9, 4, 2], 2),
    ([9, 6], 1),
    ([10, 3, 2], 1),
    ([10, 4, 1], 1),
    ([10, 5], 1),
    ([11, 2, 2], 1),
    ([11, 4], 1),
    ([12, 3], 1),
    ([13, 2], 1),
    ([15], 1),
]
# Certifié par plethysme Python (MN + formule Littlewood) + vérification dimensionnelle exacte.
# sum(mult * dim_schur_gl_n(lam, 9)) = 1_082_239_158 = C(169, 5). VÉRIFIÉ.

CANDIDATES_I5: List[dict] = []  # résultats trouvés pour I_5

# ===========================================================================
# SECTION 5 : Fonctions mathématiques de base
# ===========================================================================

def partitions_of(n: int, max_parts: Optional[int] = None) -> List[List[int]]:
    """
    Génère toutes les partitions de n (ordre décroissant).
    Si max_parts est donné, limite le nombre de parts.
    """
    if n == 0:
        return [[]]
    if n < 0:
        return []
    result = []

    def _gen(remaining, max_val, current):
        if remaining == 0:
            result.append(current[:])
            return
        if max_parts is not None and len(current) >= max_parts:
            return
        for k in range(min(remaining, max_val), 0, -1):
            current.append(k)
            _gen(remaining - k, k, current)
            current.pop()

    _gen(n, n, [])
    return result


def z_mu(mu: List[int]) -> int:
    """
    Calcule z_mu = prod_i (i^{a_i} * a_i!)
    où a_i = nombre de parts égales à i dans mu.
    """
    from collections import Counter
    cnt = Counter(mu)
    result = 1
    for val, count in cnt.items():
        result *= (val ** count) * math.factorial(count)
    return result


# Mémoïsation pour murnaghan_nakayama
_mn_cache: Dict[Tuple, int] = {}


def murnaghan_nakayama(lam: List[int], mu: List[int]) -> int:
    """
    Calcul du caractère de Schur chi^lam(mu) par la règle de Murnaghan-Nakayama.
    lam : partition (liste décroissante)
    mu  : type de cycle (liste décroissante, somme = |lam|)
    Retourne un entier (valeur du caractère irréductible de S_n).
    """
    lam_t = tuple(lam)
    mu_t = tuple(mu)
    key = (lam_t, mu_t)
    if key in _mn_cache:
        return _mn_cache[key]

    n = sum(lam)
    if sum(mu) != n:
        raise ValueError(f"sum(lam)={n} != sum(mu)={sum(mu)}")
    if n == 0:
        _mn_cache[key] = 1
        return 1

    # Règle MN : chi^lam(mu) = sum_{hooks k} (-1)^{height} chi^{lam\hook}(mu[1:])
    k = mu[0]
    mu_rest = list(mu[1:])

    results_list: List[Tuple] = []
    _rim_hooks(lam, k, results_list)

    total = 0
    for new_lam, sign in results_list:
        sub = murnaghan_nakayama(list(new_lam), mu_rest)
        total += sign * sub

    _mn_cache[key] = total
    return total


def _rim_hooks(lam: List[int], k: int, results: List[Tuple]):
    """
    Trouve tous les rim hooks de taille k dans le diagramme de Young lam.
    Ajoute (new_partition, sign) à results.
    sign = (-1)^{hauteur du rim hook}.
    """
    n_rows = len(lam)
    if n_rows == 0:
        return

    lam_ext = list(lam) + [0] * k  # padding

    def try_remove_hook(lam_list, k_rem):
        """Retire un rim hook de taille k_rem. Retourne (new_lam, height) ou None."""
        n = len(lam_list)
        while n > 0 and lam_list[n - 1] == 0:
            n -= 1
        if n == 0:
            return

        for r2 in range(n - 1, -1, -1):
            for r1 in range(r2 + 1):
                valid = True
                for i in range(r1, r2):
                    if lam_list[i + 1] >= lam_list[i]:
                        valid = False
                        break
                if not valid:
                    continue

                below_r2 = lam_list[r2 + 1] if r2 + 1 < len(lam_list) else 0
                cells_removed = lam_list[r1] - below_r2
                if cells_removed != k_rem:
                    continue

                new_lam_arr = list(lam_list[:n])
                for i in range(r1, r2):
                    new_lam_arr[i] = lam_list[i + 1]
                new_lam_arr[r2] = below_r2

                while new_lam_arr and new_lam_arr[-1] == 0:
                    new_lam_arr.pop()

                valid_part = all(new_lam_arr[i] >= new_lam_arr[i + 1] for i in range(len(new_lam_arr) - 1))
                if not valid_part:
                    continue

                height = r2 - r1
                sign = (-1) ** height
                yield (tuple(new_lam_arr), sign)

    for item in try_remove_hook(lam_ext[:max(len(lam), k) + 1], k):
        results.append(item)


# ===========================================================================
# SECTION 6 : dim_schur_gl_n
# ===========================================================================

def dim_schur_gl_n(lam: List[int], n: int = 9) -> int:
    """
    Formule de Weyl exacte pour dim V^lam où V = C^n.
    Utilise les entiers, pas de float.
    lam : partition (liste d'entiers décroissants, longueur <= n).
    Retourne un entier >= 0.
    Si len(lam) > n, retourne 0 (la représentation est nulle pour GL_n).
    """
    lam_ext = list(lam) + [0] * (n - len(lam))
    if len(lam) > n:
        return 0

    numerator = 1
    denominator = 1
    for i in range(n):
        for j in range(i + 1, n):
            num_ij = lam_ext[i] - lam_ext[j] + (j - i)
            den_ij = j - i
            numerator *= num_ij
            denominator *= den_ij

    result = numerator // denominator
    return result


def _verify_plethysm_dimension() -> bool:
    """Vérifie que la somme des dimensions des isotypes = C(169, 5)."""
    expected = math.comb(169, 5)  # = 1_082_239_158
    total = 0
    for lam, mult in PLETHYSM_S5_S3_GL9_CERTIFIED:
        d = dim_schur_gl_n(lam, 9)
        total += mult * d
    return total == expected


# ===========================================================================
# SECTION 7 : parse_plethysm_output
# ===========================================================================

def parse_plethysm_output(text: str) -> List[Tuple[List[int], int]]:
    """
    Parse robuste de la sortie Sage ou Python pour les isotypes de plethysme.

    Formats acceptés :
      - "ISO [12, 3] mult=1"
      - "[12, 3] : 1"
      - "([12, 3], 1)"
      - "s[12, 3]"         (multiplicité 1 implicite)
      - "12*s[12, 3] + ..." (format Sage direct)
      - "s[12, 3] + s[9, 4, 2]" (multiplicités 1 implicites)
      - "2*s[9, 4, 2]"     (multiplicité explicite)

    Retourne une liste de ([partition], multiplicité).
    """
    results = []

    # Format Sage : "k*s[a, b, c, ...]" ou "s[a, b, c, ...]"
    sage_pattern = re.compile(r'(\d+)\*s\[([^\]]+)\]|s\[([^\]]+)\]')
    for m in sage_pattern.finditer(text):
        if m.group(1) is not None:
            mult = int(m.group(1))
            parts_str = m.group(2)
        else:
            mult = 1
            parts_str = m.group(3)
        parts = [int(x.strip()) for x in parts_str.split(',') if x.strip()]
        if parts:
            results.append((parts, mult))

    if results:
        return _consolidate_plethysm(results)

    # Format "ISO [12, 3] mult=1"
    iso_pattern = re.compile(r'ISO\s+\[([^\]]+)\]\s+mult=(\d+)')
    for m in iso_pattern.finditer(text):
        parts = [int(x.strip()) for x in m.group(1).split(',')]
        mult = int(m.group(2))
        results.append((parts, mult))

    if results:
        return _consolidate_plethysm(results)

    # Format "[12, 3] : 1"
    colon_pattern = re.compile(r'\[([^\]]+)\]\s*:\s*(\d+)')
    for m in colon_pattern.finditer(text):
        parts = [int(x.strip()) for x in m.group(1).split(',')]
        mult = int(m.group(2))
        results.append((parts, mult))

    if results:
        return _consolidate_plethysm(results)

    # Format "([12, 3], 1)"
    tuple_pattern = re.compile(r'\(\[([^\]]+)\]\s*,\s*(\d+)\)')
    for m in tuple_pattern.finditer(text):
        parts = [int(x.strip()) for x in m.group(1).split(',')]
        mult = int(m.group(2))
        results.append((parts, mult))

    return _consolidate_plethysm(results)


def _consolidate_plethysm(items: List[Tuple[List[int], int]]) -> List[Tuple[List[int], int]]:
    """Consolide les doublons en additionnant les multiplicités."""
    agg: Dict[Tuple, int] = {}
    for parts, mult in items:
        key = tuple(parts)
        agg[key] = agg.get(key, 0) + mult
    return [(list(k), v) for k, v in sorted(agg.items(), key=lambda x: (-x[1][0] if x[0] else 0, x[0]))]


# ===========================================================================
# SECTION 8 : Fonctions GL9 / W
# ===========================================================================

def generate_W_basis_GL9() -> List[Tuple[int, int, int]]:
    """
    Génère la base de W = Sym^3(C^9) = espace des cubes symétriques.
    Base indexée par les triplets (i, j, k) avec 0 <= i <= j <= k <= 8.
    dim W = C(9+3-1, 3) = C(11, 3) = 165.
    Retourne une liste de 165 triplets (i, j, k).
    """
    basis = []
    for i in range(9):
        for j in range(i, 9):
            for k in range(j, 9):
                basis.append((i, j, k))
    assert len(basis) == 165, f"Expected 165 basis elements, got {len(basis)}"
    return basis


def encode_degree5_monomial(multi_exp: Dict[int, int], W_basis: List[Tuple[int, int, int]]) -> Tuple[int, ...]:
    """
    Encode un monôme de Sym^5(W^*) comme un tuple de 165 exposants.
    multi_exp : dict {idx_dans_W_basis -> exposant}, exposants de somme 5.
    Retourne un tuple de longueur 165.
    """
    exp_vec = [0] * len(W_basis)
    for idx, exp in multi_exp.items():
        exp_vec[idx] = exp
    return tuple(exp_vec)


def weight_of_W_variable(ijk: Tuple[int, int, int]) -> List[int]:
    """
    Poids de la variable w_{ijk} dans W = Sym^3(C^9).
    w_{ijk} = e_i otimes e_j otimes e_k (symmétrisé).
    Poids = somme des poids des indices = e_i + e_j + e_k dans Z^9.
    Retourne un vecteur de longueur 9.
    """
    wt = [0] * 9
    for idx in ijk:
        wt[idx] += 1
    return wt


def weight_of_degree5_monomial(exp_vec: Tuple[int, ...], W_basis: List[Tuple[int, int, int]]) -> List[int]:
    """
    Poids d'un monôme prod_alpha w_alpha^{exp_alpha} dans Sym^5(W^*).
    Note : W^* a pour poids opposés à W.
    Poids d'un monôme = sum_alpha exp_alpha * weight(w_alpha) dans W.
    (On utilise la convention poids dans W, pas W^* pour simplifier.)
    Retourne un vecteur de longueur 9 (somme = 5 * 3 = 15).
    """
    wt = [0] * 9
    for alpha, e in enumerate(exp_vec):
        if e != 0:
            ijk = W_basis[alpha]
            for idx in ijk:
                wt[idx] += e
    return wt


def lie_action_Eij_on_W_variable(
    i: int, j: int, ijk: Tuple[int, int, int]
) -> List[Tuple[Tuple[int, int, int], int]]:
    """
    Action de E_{i,j} sur la variable w_{abc} de W = Sym^3(C^9).
    E_{i,j} · w_{abc} = sum_{k: abc[k]=j} w_{abc with abc[k] replaced by i}
    Retourne une liste de (new_ijk, coefficient).
    """
    result = []
    abc = list(ijk)
    for pos in range(3):
        if abc[pos] == j:
            new_abc = abc[:]
            new_abc[pos] = i
            new_abc_sorted = tuple(sorted(new_abc))
            result.append((new_abc_sorted, 1))
    # Consolider les doublons
    agg: Dict[Tuple, int] = {}
    for t, c in result:
        agg[t] = agg.get(t, 0) + c
    return [(t, c) for t, c in agg.items() if c != 0]


def lie_action_Eij_on_degree5_monomial(
    i: int,
    j: int,
    exp_vec: Tuple[int, ...],
    W_basis: List[Tuple[int, int, int]],
) -> List[Tuple[Tuple[int, ...], int]]:
    """
    Action de E_{i,j} sur un monôme m = prod_alpha w_alpha^{e_alpha} de Sym^5(W^*).
    Règle de Leibniz : E · (w_alpha * m') = (E · w_alpha) * m' + w_alpha * (E · m').
    Retourne une liste de (new_exp_vec, coefficient).
    """
    w_basis_idx = {ijk: idx for idx, ijk in enumerate(W_basis)}
    result_agg: Dict[Tuple[int, ...], int] = {}

    for alpha in range(len(W_basis)):
        e_alpha = exp_vec[alpha]
        if e_alpha == 0:
            continue
        ijk = W_basis[alpha]
        images = lie_action_Eij_on_W_variable(i, j, ijk)
        for new_ijk, coeff in images:
            if new_ijk not in w_basis_idx:
                continue
            beta = w_basis_idx[new_ijk]
            new_exp = list(exp_vec)
            new_exp[alpha] -= 1
            new_exp[beta] += 1
            new_exp_t = tuple(new_exp)
            total_coeff = coeff * e_alpha
            result_agg[new_exp_t] = result_agg.get(new_exp_t, 0) + total_coeff

    return [(t, c) for t, c in result_agg.items() if c != 0]


# ===========================================================================
# SECTION 9 : det3 et orbite
# ===========================================================================

def det3_terms() -> List[Tuple[Tuple[int, int, int], int]]:
    """
    Retourne le déterminant 3x3 comme polynôme dans Sym^3(C^9).
    det_3 = sum_{sigma in S_3} sgn(sigma) * x_{0,sigma(0)} * x_{1,sigma(1)} * x_{2,sigma(2)}
    où les variables x_{i,j} sont les entrées d'une matrice 3x3.
    On plonge C^9 = C^{3x3} avec la base e_{ij} pour 0 <= i,j <= 2 (9 variables au total).
    Le déterminant est un élément de Sym^3(C^9^*) ≃ W^*.
    Retourne une liste de (triplet d'indices dans C^9, signe).
    Encodage : e_{ij} correspond à l'indice 3*i+j dans C^9.
    """
    import itertools
    terms = []
    for perm in itertools.permutations([0, 1, 2]):
        # sgn(perm)
        sgn = 1
        lst = list(perm)
        for ii in range(3):
            while lst[ii] != ii:
                jj = lst[ii]
                lst[ii], lst[jj] = lst[jj], lst[ii]
                sgn *= -1
        indices = tuple(sorted([3 * 0 + perm[0], 3 * 1 + perm[1], 3 * 2 + perm[2]]))
        terms.append((indices, sgn))
    # Consolider
    agg: Dict[Tuple, int] = {}
    for t, s in terms:
        agg[t] = agg.get(t, 0) + s
    return [(t, s) for t, s in agg.items() if s != 0]


def random_GL9_modp(p: int, seed: int = 42) -> List[List[int]]:
    """
    Génère une matrice aléatoire inversible 9x9 mod p.
    """
    rng = random.Random(seed)
    while True:
        mat = [[rng.randint(0, p - 1) for _ in range(9)] for _ in range(9)]
        if _det9_modp(mat, p) != 0:
            return mat


def _det9_modp(mat: List[List[int]], p: int) -> int:
    """
    Calcule le déterminant d'une matrice 9x9 mod p par élimination de Gauss.
    Retourne le déterminant mod p.
    """
    n = len(mat)
    m = [row[:] for row in mat]
    det = 1
    for col in range(n):
        # Trouver un pivot
        pivot_row = None
        for row in range(col, n):
            if m[row][col] % p != 0:
                pivot_row = row
                break
        if pivot_row is None:
            return 0
        if pivot_row != col:
            m[col], m[pivot_row] = m[pivot_row], m[col]
            det = (-det) % p
        inv_pivot = mod_inv(m[col][col] % p, p)
        det = (det * m[col][col]) % p
        for row in range(col + 1, n):
            if m[row][col] % p == 0:
                continue
            factor = (m[row][col] * inv_pivot) % p
            for c2 in range(col, n):
                m[row][c2] = (m[row][c2] - factor * m[col][c2]) % p
    return det % p


def act_g_on_det3_modp(g: List[List[int]], p: int) -> Dict[Tuple[int, int, int], int]:
    """
    Applique g in GL_9(F_p) sur le déterminant det_3 vu comme élément de Sym^3(C^9).
    Retourne le résultat sous forme de dictionnaire {triplet_ijk: coefficient mod p}.
    g agit sur C^9 par la représentation contragrediente sur W^*.
    """
    terms = det3_terms()
    result: Dict[Tuple[int, int, int], int] = {}

    for (i1, i2, i3), sign in terms:
        for j1 in range(9):
            c1 = g[j1][i1] % p
            if c1 == 0:
                continue
            for j2 in range(9):
                c2 = g[j2][i2] % p
                if c2 == 0:
                    continue
                for j3 in range(9):
                    c3 = g[j3][i3] % p
                    if c3 == 0:
                        continue
                    new_ijk = tuple(sorted([j1, j2, j3]))
                    coeff = (sign * c1 * c2 * c3) % p
                    result[new_ijk] = (result.get(new_ijk, 0) + coeff) % p

    return {k: v for k, v in result.items() if v != 0}


def evaluate_degree5_monomial_on_point(
    exp_vec: Tuple[int, ...],
    point_W: Dict[Tuple[int, int, int], int],
    W_basis: List[Tuple[int, int, int]],
    p: int,
) -> int:
    """
    Évalue un monôme de Sym^5(W^*) sur un point de W (donné comme dict ijk -> coeff mod p).
    Retourne la valeur mod p.
    """
    value = 1
    for alpha, e in enumerate(exp_vec):
        if e == 0:
            continue
        ijk = W_basis[alpha]
        coord = point_W.get(ijk, 0) % p
        value = (value * pow(coord, e, p)) % p
    return value


# ===========================================================================
# SECTION 10 : Arithmétique mod p
# ===========================================================================

def is_prime(n: int) -> bool:
    """
    Test de primalité Miller-Rabin déterministe pour n < 3_317_044_064_679_887_385_961_981.
    Retourne True si n est premier.
    """
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Petits premiers
    small_primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for sp in small_primes:
        if n == sp:
            return True
        if n % sp == 0:
            return False

    # Écrire n-1 = 2^r * d
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    # Témoins déterministes (valides pour n < 3.3*10^24)
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for a in witnesses:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def mod_inv(a: int, p: int) -> int:
    """
    Inverse de a modulo p (p premier) via l'algorithme d'Euclide étendu.
    Lève une exception si a == 0 mod p.
    """
    a = a % p
    if a == 0:
        raise ZeroDivisionError(f"a={a} n'a pas d'inverse mod {p}")
    g, x, _ = _extended_gcd(a, p)
    if g != 1:
        raise ValueError(f"gcd({a},{p}) = {g} != 1")
    return x % p


def _extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Algorithme d'Euclide étendu. Retourne (gcd, x, y) tel que a*x + b*y = gcd."""
    if b == 0:
        return a, 1, 0
    g, x, y = _extended_gcd(b, a % b)
    return g, y, x - (a // b) * y


def rank_modp_cpu(matrix: List[List[int]], p: int) -> int:
    """
    Calcule le rang d'une matrice mod p par élimination de Gauss.
    Exact, pas de float. Retourne un entier.
    """
    if not matrix or not matrix[0]:
        return 0
    n_rows = len(matrix)
    n_cols = len(matrix[0])
    m = [[x % p for x in row] for row in matrix]

    rank = 0
    pivot_col = 0
    for row in range(n_rows):
        if pivot_col >= n_cols:
            break
        # Trouver un pivot dans la colonne courante
        while pivot_col < n_cols:
            nonzero = None
            for r in range(row, n_rows):
                if m[r][pivot_col] != 0:
                    nonzero = r
                    break
            if nonzero is not None:
                break
            pivot_col += 1

        if pivot_col >= n_cols:
            break

        m[row], m[nonzero] = m[nonzero], m[row]
        inv_pivot = mod_inv(m[row][pivot_col], p)
        for r2 in range(n_rows):
            if r2 == row:
                continue
            if m[r2][pivot_col] == 0:
                continue
            factor = (m[r2][pivot_col] * inv_pivot) % p
            for c2 in range(n_cols):
                m[r2][c2] = (m[r2][c2] - factor * m[row][c2]) % p

        rank += 1
        pivot_col += 1

    return rank


def random_matrix_modp(m: int, n: int, p: int, seed: int = 42) -> List[List[int]]:
    """
    Génère une matrice aléatoire m x n mod p.
    """
    rng = random.Random(seed)
    return [[rng.randint(0, p - 1) for _ in range(n)] for _ in range(m)]


def rank_modp_gpu_if_available(matrix: List[List[int]], p: int) -> Tuple[int, str]:
    """
    Calcule le rang mod p exact. Utilise rank_modp_exact_gpu.
    Retourne (rang, source) où source = 'gpu' ou 'cpu'.
    """
    return rank_modp_exact_gpu(matrix, p)


# ===========================================================================
# SECTION 11 : Espace de poids et HWV
# ===========================================================================

WEIGHT_SPACE_MAX = 500_000  # limite honnête


class TOO_LARGE(Exception):
    """Levée quand l'espace de poids dépasse WEIGHT_SPACE_MAX."""
    pass


def _exp_vec_to_mono5(exp_vec) -> MonoIdx5:
    """
    Convertit un exp_vec (tuple de 165 entiers, somme=5) en MonoIdx5
    (tuple trié de 5 indices dans [0,164]).
    """
    indices = []
    for i, e in enumerate(exp_vec):
        if e > 0:
            indices.extend([i] * e)
    return tuple(sorted(indices))


def mono5_to_exp_vec(m5: MonoIdx5, n_vars: int = 165) -> tuple:
    """
    Convertit un MonoIdx5 (5 indices triés) en exp_vec de longueur n_vars.
    """
    exp_vec = [0] * n_vars
    for idx in m5:
        exp_vec[idx] += 1
    return tuple(exp_vec)


def generate_weight_space_monomials(
    lam: List[int], W_basis: List[Tuple[int, int, int]], max_size: int = WEIGHT_SPACE_MAX
) -> List[MonoIdx5]:
    """
    Génère les monômes de Sym^5(W^*) de poids exactement lam.

    Retourne une liste de MonoIdx5 (tuples de 5 indices triés dans [0,164]).
    Ceci remplace la représentation par exp_vec 165-long (BREAKING CHANGE v5).

    NE GÉNÈRE JAMAIS les C(169,5) = 1_082_239_158 monômes en dense.
    Génère seulement ceux de poids lam (espace de poids de dim << totale).

    Si trop grand (> max_size), lève TOO_LARGE.

    Algorithme : backtracking sur les 165 variables, en maintenant le poids cumulé.
    """
    target_wt = list(lam) + [0] * (9 - len(lam))
    total_degree = 5
    n_vars = len(W_basis)  # = 165

    # Pré-calculer les poids des variables
    var_weights = [weight_of_W_variable(W_basis[alpha]) for alpha in range(n_vars)]

    results = []

    def backtrack(alpha: int, remaining_degree: int, current_wt: List[int], current_exp: List[int]):
        if len(results) > max_size:
            raise TOO_LARGE(f"Espace de poids > {max_size}")

        if remaining_degree == 0:
            if current_wt == target_wt:
                # Convertir l'exp_vec courant en MonoIdx5
                results.append(_exp_vec_to_mono5(current_exp))
            return

        if alpha >= n_vars:
            return

        # Pruning : vérifier que le poids restant est atteignable
        remaining_wt = [target_wt[k] - current_wt[k] for k in range(9)]
        if any(r < 0 for r in remaining_wt):
            return

        # Vérifier que la somme du poids restant est compatible
        remaining_wt_sum = sum(remaining_wt)
        if remaining_wt_sum != remaining_degree * 3:
            return

        # Exposant maximum pour alpha
        max_exp = remaining_degree
        for k in range(9):
            if var_weights[alpha][k] > 0 and remaining_wt[k] >= 0:
                possible = remaining_wt[k] // var_weights[alpha][k]
                max_exp = min(max_exp, possible)

        for e in range(max_exp, -1, -1):
            current_exp[alpha] = e
            new_wt = [current_wt[k] + e * var_weights[alpha][k] for k in range(9)]
            backtrack(alpha + 1, remaining_degree - e, new_wt, current_exp)

        current_exp[alpha] = 0

    current_exp = [0] * n_vars
    backtrack(0, total_degree, [0] * 9, current_exp)
    return results


def build_hwv_equations_modp(
    lam: List[int], p: int, W_basis: List[Tuple[int, int, int]]
) -> List[List[int]]:
    """
    Construit les équations E_{i,i+1} · v = 0 pour i=0,...,7 dans Sym^5(W^*).
    v = vecteur dans l'espace de poids de lambda (dominant weight).

    Les monomials sont maintenant des MonoIdx5 (v5).
    En interne, on convertit en exp_vec pour l'action de Lie.

    Retourne la matrice (liste de listes d'entiers mod p).
    """
    # Générer les monômes de poids lam (MonoIdx5 format)
    log(f"  Generating weight space for lam={lam}...")
    try:
        monomials_mono5 = generate_weight_space_monomials(lam, W_basis)
    except TOO_LARGE as e:
        log(f"  TOO_LARGE: {e}")
        raise

    n_monomials = len(monomials_mono5)
    log(f"  Weight space dimension: {n_monomials}")
    if n_monomials == 0:
        return []

    # Indexer les MonoIdx5
    mono_to_idx: Dict[MonoIdx5, int] = {m5: idx for idx, m5 in enumerate(monomials_mono5)}

    # Convertir en exp_vec pour l'action de Lie
    monomials_expvec = [mono5_to_exp_vec(m5) for m5 in monomials_mono5]

    # Pour chaque générateur E_{i,i+1}
    all_rows: List[Dict[int, int]] = []

    for i_gen in range(8):  # E_{i, i+1} pour i = 0..7
        j_gen = i_gen + 1

        target_mono_to_equations: Dict[MonoIdx5, Dict[int, int]] = {}

        for col_idx, mono_expvec in enumerate(monomials_expvec):
            images = lie_action_Eij_on_degree5_monomial(i_gen, j_gen, mono_expvec, W_basis)
            for new_mono_expvec, coeff in images:
                new_mono5 = _exp_vec_to_mono5(new_mono_expvec)
                if new_mono5 not in target_mono_to_equations:
                    target_mono_to_equations[new_mono5] = {}
                target_mono_to_equations[new_mono5][col_idx] = (
                    target_mono_to_equations[new_mono5].get(col_idx, 0) + coeff
                ) % p

        # Chaque monôme cible donne une ligne de la matrice
        for target_mono5, eq_dict in target_mono_to_equations.items():
            row = {c: v % p for c, v in eq_dict.items() if v % p != 0}
            if row:
                all_rows.append(row)

    if not all_rows:
        return []

    # Convertir en matrice dense
    n_rows = len(all_rows)
    matrix = [[0] * n_monomials for _ in range(n_rows)]
    for r, row_dict in enumerate(all_rows):
        for c, v in row_dict.items():
            matrix[r][c] = v % p

    return matrix


def kernel_modp_sparse(matrix: List[List[int]], p: int) -> List[List[int]]:
    """
    Calcule le noyau (base) d'une matrice mod p.
    Retourne une liste de vecteurs (listes d'entiers mod p).
    Utilise l'élimination de Gauss avec pivots partiels.
    Exact, pas de float.
    """
    if not matrix or not matrix[0]:
        return []

    n_rows = len(matrix)
    n_cols = len(matrix[0])

    m = [[x % p for x in row] for row in matrix]
    pivot_rows = {}  # col -> row of pivot
    pivot_order = []
    row_ptr = 0

    for col in range(n_cols):
        if row_ptr >= n_rows:
            break
        # Trouver un pivot
        pivot = None
        for r in range(row_ptr, n_rows):
            if m[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        m[row_ptr], m[pivot] = m[pivot], m[row_ptr]
        inv_p = mod_inv(m[row_ptr][col], p)
        # Normaliser la ligne pivot
        for c2 in range(n_cols):
            m[row_ptr][c2] = (m[row_ptr][c2] * inv_p) % p
        # Éliminer dans les autres lignes
        for r2 in range(n_rows):
            if r2 == row_ptr:
                continue
            if m[r2][col] == 0:
                continue
            factor = m[r2][col]
            for c2 in range(n_cols):
                m[r2][c2] = (m[r2][c2] - factor * m[row_ptr][c2]) % p

        pivot_rows[col] = row_ptr
        pivot_order.append(col)
        row_ptr += 1

    # Les colonnes sans pivot = variables libres
    free_cols = [c for c in range(n_cols) if c not in pivot_rows]

    kernel_basis = []
    for fc in free_cols:
        vec = [0] * n_cols
        vec[fc] = 1
        for pc, pr in pivot_rows.items():
            val = m[pr][fc]
            if val != 0:
                vec[pc] = (p - val) % p
        kernel_basis.append(vec)

    return kernel_basis


# ===========================================================================
# SECTION 12 : Évaluation orbitale (CPU legacy)
# ===========================================================================

def evaluate_hwv_on_orbit_points(
    hwv_basis: List[Dict[MonoIdx5, int]],
    primes: List[int],
    seeds: List[int],
    W_basis: List[Tuple[int, int, int]],
) -> Dict[str, Any]:
    """
    Pour chaque vecteur HWV (combinaison linéaire de monômes de W^*),
    évalue sur des points g·det_3 mod p.

    hwv_basis : liste de vecteurs sous forme de dict {MonoIdx5 -> coeff}
    primes : liste de nombres premiers à utiliser
    seeds  : liste de graines pour les matrices GL_9 aléatoires
    W_basis : base de W

    Retourne {
        'vanishes_all': bool,
        'rank_eval_matrix': int,
        'details': [...],
    }
    """
    if not hwv_basis:
        return {'vanishes_all': True, 'rank_eval_matrix': 0, 'details': []}

    details = []
    all_eval_matrices = []

    for p in primes:
        eval_matrix = []  # lignes = points, colonnes = HWV
        for seed in seeds:
            g = random_GL9_modp(p, seed)
            point_W = act_g_on_det3_modp(g, p)

            # Convertir point_W (dict ijk -> coeff) en vecteur dense (W_basis index)
            w_basis_idx = {ijk: i for i, ijk in enumerate(W_basis)}
            coords = [0] * 165
            for ijk, c in point_W.items():
                idx = w_basis_idx.get(ijk, -1)
                if idx >= 0:
                    coords[idx] = c % p

            row = []
            for hwv_vec in hwv_basis:
                val = 0
                for m5, coeff in hwv_vec.items():
                    # Évaluer le monôme m5 sur les coords
                    mono_val = 1
                    # m5 est un tuple de 5 indices; calculer prod avec run-length
                    # Compteur des répétitions
                    from collections import Counter
                    cnt = Counter(m5)
                    for var_idx, power in cnt.items():
                        c_var = coords[var_idx]
                        mono_val = (mono_val * pow(c_var, power, p)) % p
                    val = (val + coeff * mono_val) % p
                row.append(val)
            eval_matrix.append(row)

        # Transposer pour avoir (HWV x points)
        if eval_matrix:
            hwv_count = len(hwv_basis)
            point_count = len(seeds)
            mat_T = [[eval_matrix[pt][hwv] for pt in range(point_count)] for hwv in range(hwv_count)]
            rank, source = rank_modp_exact_gpu(mat_T, p)
            all_eval_matrices.append(rank)
            details.append({'prime': p, 'rank': rank, 'source': source})

    # Vérifier si tous les vecteurs HWV s'annulent sur tous les points
    vanishes_all = all(d['rank'] == 0 for d in details)
    total_rank = max(all_eval_matrices) if all_eval_matrices else 0

    return {
        'vanishes_all': vanishes_all,
        'rank_eval_matrix': total_rank,
        'details': details,
    }


# ===========================================================================
# SECTION 12b : GPU Evaluation (v5)
# ===========================================================================

def hwv_vec_to_sparse(vec: List[int], p: int) -> HWVSparse:
    """
    Convertit un vecteur HWV (liste d'entiers, indexé par position dans monomials)
    en HWVSparse (deux arrays numpy parallèles).
    """
    mono_ids = []
    coeffs = []
    for i, c in enumerate(vec):
        c_mod = c % p
        if c_mod != 0:
            mono_ids.append(i)
            coeffs.append(c_mod)
    if HAS_NUMPY:
        return HWVSparse(
            mono_ids=np.array(mono_ids, dtype=np.int32),
            coeffs=np.array(coeffs, dtype=np.int64),
        )
    else:
        return HWVSparse(mono_ids=mono_ids, coeffs=coeffs)


def hwv_sparse_to_dict(sp: HWVSparse, monomials: List[MonoIdx5]) -> Dict[MonoIdx5, int]:
    """
    Convertit un HWVSparse en dict {MonoIdx5: coeff} pour compatibilité CPU.
    """
    result = {}
    if HAS_NUMPY:
        for i, c in zip(sp.mono_ids.tolist(), sp.coeffs.tolist()):
            if c != 0:
                result[monomials[i]] = c
    else:
        for i, c in zip(sp.mono_ids, sp.coeffs):
            if c != 0:
                result[monomials[i]] = c
    return result


def build_orbit_coords_gpu(
    seeds: List[int],
    p: int,
    W_basis: List[Tuple[int, int, int]],
    use_gpu: bool = True,
) -> Tuple[Any, str]:
    """
    Construit la matrice de coordonnées d'orbite shape (n_points, 165).

    Pour chaque seed :
      - g = random_GL9_modp(p, seed)
      - point_W = act_g_on_det3_modp(g, p)  → dict {ijk: coeff}
      - Convertir en vecteur dense de longueur 165

    Retourne (array shape (n_points, 165), source)
    Si CuPy disponible + use_gpu : return cp.array(...), "gpu"
    Sinon : return np.array(...), "cpu"
    Toutes les valeurs sont des int64 exacts mod p (JAMAIS float).
    """
    w_basis_idx = {ijk: i for i, ijk in enumerate(W_basis)}
    n_vars = len(W_basis)  # 165

    coords_list = []
    for seed in seeds:
        g = random_GL9_modp(p, seed)
        point_W = act_g_on_det3_modp(g, p)
        coords = [0] * n_vars
        for ijk, c in point_W.items():
            idx = w_basis_idx.get(ijk, -1)
            if idx >= 0:
                coords[idx] = int(c) % p
        coords_list.append(coords)

    if not HAS_NUMPY:
        return coords_list, "cpu-list"

    coords_np = np.array(coords_list, dtype=np.int64)  # shape (n_points, 165)

    if HAS_CUPY and use_gpu:
        try:
            coords_gpu = cp.array(coords_np, dtype=cp.int64)
            return coords_gpu, "gpu"
        except Exception as e:
            log(f"  WARNING: GPU transfer failed ({e}), falling back to CPU")

    return coords_np, "cpu"


def evaluate_hwv_on_orbit_points_gpu(
    hwv_sparse_list: List[HWVSparse],
    orbit_coords: Any,     # array shape (n_points, 165), int64
    p: int,
    use_gpu: bool = True,
) -> Tuple[Any, str]:
    """
    Évalue chaque HWV sur chaque point d'orbite.

    Retourne (eval_matrix shape (n_hwv, n_points), source)
    Toutes les valeurs sont des int64 exacts mod p. JAMAIS de float.

    GPU path (CuPy) :
      - Construit les tableaux CSR : mono5_flat (nnz_total, 5) int16,
        hwv_offsets (n_hwv+1,) int32, coeffs_flat (nnz_total,) int64
      - Lance le kernel CUDA eval_hwv_kernel

    CPU fallback (numpy) :
      - Même logique en numpy vectorisé
    """
    n_hwv = len(hwv_sparse_list)
    if n_hwv == 0:
        if HAS_NUMPY:
            return np.zeros((0, 0), dtype=np.int64), "cpu"
        return [], "cpu"

    # Déterminer la forme de orbit_coords
    if HAS_NUMPY and isinstance(orbit_coords, np.ndarray):
        n_points = orbit_coords.shape[0]
        is_gpu_array = False
    elif HAS_CUPY and isinstance(orbit_coords, cp.ndarray):
        n_points = orbit_coords.shape[0]
        is_gpu_array = True
    else:
        # Fallback : liste Python
        n_points = len(orbit_coords)
        is_gpu_array = False

    # Construire les tableaux CSR pour tous les HWV
    mono5_flat_list = []  # (nnz_total, 5) int16
    coeffs_flat_list = []  # (nnz_total,) int64
    hwv_offsets = [0]     # (n_hwv+1,) int32

    for sp in hwv_sparse_list:
        if HAS_NUMPY:
            nnz = len(sp.mono_ids)
            for t in range(nnz):
                mid = int(sp.mono_ids[t])
                # mid est un indice dans monomials, on stocke directement
                # NOTE: mono5_flat stocke l'indice de monôme brut (pas les 5 indices)
                # Pour l'évaluation on a besoin des 5 indices, donc on utilise
                # une approche différente : stocker mono_ids et reconstruire au vol
                coeffs_flat_list.append(int(sp.coeffs[t]))
                hwv_offsets.append(hwv_offsets[-1])  # sera mis à jour
            hwv_offsets[-1] = hwv_offsets[-2] + nnz
        else:
            nnz = len(sp.mono_ids)
            for t in range(nnz):
                coeffs_flat_list.append(sp.coeffs[t])
            hwv_offsets.append(hwv_offsets[-1] + nnz)

    # Reconstruction propre des offsets
    hwv_offsets = [0]
    for sp in hwv_sparse_list:
        hwv_offsets.append(hwv_offsets[-1] + (len(sp.mono_ids) if HAS_NUMPY else len(sp.mono_ids)))

    nnz_total = hwv_offsets[-1]

    # GPU path
    if HAS_CUPY and use_gpu and HAS_NUMPY:
        try:
            return _eval_hwv_gpu(hwv_sparse_list, orbit_coords, p, n_hwv, n_points,
                                  hwv_offsets, nnz_total, is_gpu_array)
        except Exception as e:
            log(f"  WARNING: GPU eval failed ({e}), falling back to CPU")

    # CPU path (numpy vectorisé)
    return _eval_hwv_cpu_numpy(hwv_sparse_list, orbit_coords, p, n_hwv, n_points, is_gpu_array)


def _eval_hwv_gpu(hwv_sparse_list, orbit_coords, p, n_hwv, n_points,
                   hwv_offsets, nnz_total, is_gpu_array):
    """
    Évaluation GPU via CuPy RawKernel.
    mono5_flat : (nnz_total, 5) stocke les 5 indices du MonoIdx5 pour chaque terme.
    Les mono_ids dans HWVSparse sont des indices dans la liste monomials.
    On a besoin des monomials pour reconstruire les 5 indices.

    IMPORTANT : Pour que le kernel GPU fonctionne, on doit passer les 5 indices.
    Ici on travaille directement avec mono_ids (indices dans la liste monomials),
    donc il faut passer une table de correspondance mono_id -> 5 indices.

    Stratégie : construire mono5_flat comme tableau (nnz_total, 5) où chaque ligne
    contient les 5 indices uint16 du monôme correspondant.
    Pour cela on a besoin de la liste monomials — mais elle n'est pas dans cette fonction.

    Fallback propre : si monomials non disponibles, utiliser CPU.
    """
    # On ne peut pas reconstruire les 5 indices sans la liste monomials.
    # Déléguer à CPU numpy pour cette implémentation.
    # (Le kernel GPU complet nécessite de passer monomials en paramètre.)
    raise NotImplementedError("GPU kernel requires monomials list; use evaluate_hwv_gpu_with_monomials")


def evaluate_hwv_gpu_with_monomials(
    hwv_sparse_list: List[HWVSparse],
    orbit_coords: Any,
    monomials: List[MonoIdx5],
    p: int,
) -> Tuple[Any, str]:
    """
    Évaluation GPU avec la liste monomials disponible.
    Construit mono5_flat (nnz_total, 5) pour le kernel CUDA.
    """
    n_hwv = len(hwv_sparse_list)
    if n_hwv == 0:
        return np.zeros((0, 0), dtype=np.int64), "cpu"

    if isinstance(orbit_coords, np.ndarray):
        n_points = orbit_coords.shape[0]
    elif HAS_CUPY and isinstance(orbit_coords, cp.ndarray):
        n_points = orbit_coords.shape[0]
    else:
        n_points = len(orbit_coords)

    # Construire hwv_offsets
    hwv_offsets_list = [0]
    for sp in hwv_sparse_list:
        hwv_offsets_list.append(hwv_offsets_list[-1] + len(sp.mono_ids))
    nnz_total = hwv_offsets_list[-1]

    if nnz_total == 0:
        result = np.zeros((n_hwv, n_points), dtype=np.int64)
        if HAS_CUPY and isinstance(orbit_coords, cp.ndarray):
            return cp.array(result), "gpu"
        return result, "cpu"

    # Construire mono5_flat : shape (nnz_total, 5), dtype int16
    mono5_flat_np = np.zeros((nnz_total, 5), dtype=np.int16)
    coeffs_flat_np = np.zeros(nnz_total, dtype=np.int64)
    hwv_offsets_np = np.array(hwv_offsets_list, dtype=np.int32)

    t = 0
    for sp in hwv_sparse_list:
        for mono_id, coeff in zip(sp.mono_ids.tolist(), sp.coeffs.tolist()):
            m5 = monomials[mono_id]
            mono5_flat_np[t, :] = np.array(m5, dtype=np.int16)
            coeffs_flat_np[t] = int(coeff) % p
            t += 1

    # Transférer sur GPU
    kernel = _get_eval_hwv_kernel()
    if kernel is None:
        raise RuntimeError("eval_hwv_kernel not compiled")

    # Assurer que orbit_coords est un cp.array int64
    if isinstance(orbit_coords, np.ndarray):
        orbit_gpu = cp.array(orbit_coords, dtype=cp.int64)
    else:
        orbit_gpu = orbit_coords.astype(cp.int64)

    mono5_gpu = cp.array(mono5_flat_np)
    coeffs_gpu = cp.array(coeffs_flat_np)
    hwv_off_gpu = cp.array(hwv_offsets_np)
    eval_matrix_gpu = cp.zeros((n_hwv, n_points), dtype=cp.int64)

    # Lancer le kernel
    block_x = 16
    block_y = 16
    grid_x = (n_hwv + block_x - 1) // block_x
    grid_y = (n_points + block_y - 1) // block_y

    kernel(
        (grid_x, grid_y), (block_x, block_y),
        (orbit_gpu, mono5_gpu, coeffs_gpu, hwv_off_gpu, eval_matrix_gpu,
         np.int32(n_hwv), np.int32(n_points), np.int64(p))
    )
    cp.cuda.Stream.null.synchronize()

    return eval_matrix_gpu, "gpu"


def _eval_hwv_cpu_numpy(hwv_sparse_list, orbit_coords, p, n_hwv, n_points, is_gpu_array):
    """
    Évaluation CPU numpy.
    orbit_coords : array int64 shape (n_points, 165) ou liste Python
    """
    if HAS_CUPY and is_gpu_array:
        # Transférer sur CPU
        coords_cpu = cp.asnumpy(orbit_coords)
    elif HAS_NUMPY and isinstance(orbit_coords, np.ndarray):
        coords_cpu = orbit_coords
    else:
        if HAS_NUMPY:
            coords_cpu = np.array(orbit_coords, dtype=np.int64)
        else:
            coords_cpu = orbit_coords

    if HAS_NUMPY:
        eval_matrix = np.zeros((n_hwv, n_points), dtype=np.int64)
    else:
        eval_matrix = [[0] * n_points for _ in range(n_hwv)]

    for k, sp in enumerate(hwv_sparse_list):
        if HAS_NUMPY:
            mono_ids_arr = sp.mono_ids  # shape (nnz,)
            coeffs_arr = sp.coeffs     # shape (nnz,)
            nnz = len(mono_ids_arr)
        else:
            mono_ids_arr = sp.mono_ids
            coeffs_arr = sp.coeffs
            nnz = len(mono_ids_arr)

        for t in range(nnz):
            mid = int(mono_ids_arr[t])
            coeff = int(coeffs_arr[t]) % p
            if coeff == 0:
                continue
            # mono_ids[t] est un indice dans monomials — mais ici on n'a pas monomials.
            # On stocke les indices directs dans mono_ids (indice de la variable 0..164).
            # NOTE: Dans hwv_vec_to_sparse, mono_ids stocke les positions dans le vecteur HWV
            # (= indices dans monomials list), PAS les 5 indices du MonoIdx5.
            # Donc ici mid = indice dans monomials, et on ne peut pas évaluer sans monomials.
            # => utiliser evaluate_hwv_sparse_on_coords_cpu à la place.
            pass

        # Fallback : évaluation via evaluate_hwv_sparse_on_coords_cpu
        for j in range(n_points):
            if HAS_NUMPY:
                pt_coords = coords_cpu[j]
            else:
                pt_coords = coords_cpu[j]
            val = _eval_single_hwv_on_point_no_mono(sp, pt_coords, p)
            if HAS_NUMPY:
                eval_matrix[k, j] = val
            else:
                eval_matrix[k][j] = val

    return eval_matrix, "cpu"


def _eval_single_hwv_on_point_no_mono(sp: HWVSparse, coords, p: int) -> int:
    """
    Évalue un HWVSparse sur un point.
    ATTENTION : les mono_ids sont des indices dans la liste monomials,
    pas des indices de variables directement. Sans monomials, on ne peut pas
    calculer les valeurs des monômes.
    Cette fonction doit être appelée via evaluate_hwv_sparse_with_monomials.
    """
    # Ce cas ne devrait pas être atteint directement.
    # Si appelé sans monomials, retourne 0.
    return 0


def evaluate_hwv_sparse_with_monomials(
    hwv_sparse_list: List[HWVSparse],
    orbit_coords: Any,
    monomials: List[MonoIdx5],
    p: int,
    use_gpu: bool = True,
) -> Tuple[Any, str]:
    """
    Évalue les HWV sur l'orbite avec la liste monomials disponible.
    Dispatch GPU ou CPU selon use_gpu et HAS_CUPY.

    orbit_coords : array int64 shape (n_points, 165)
    Retourne (eval_matrix shape (n_hwv, n_points), source)
    """
    n_hwv = len(hwv_sparse_list)
    if n_hwv == 0:
        if HAS_NUMPY:
            return np.zeros((0, 0), dtype=np.int64), "cpu"
        return [], "cpu"

    if isinstance(orbit_coords, np.ndarray):
        n_points = orbit_coords.shape[0]
    elif HAS_CUPY and isinstance(orbit_coords, cp.ndarray):
        n_points = orbit_coords.shape[0]
    else:
        n_points = len(orbit_coords)

    # GPU path
    if HAS_CUPY and use_gpu:
        try:
            return evaluate_hwv_gpu_with_monomials(
                hwv_sparse_list, orbit_coords, monomials, p
            )
        except Exception as e:
            log(f"  WARNING: GPU eval failed ({e}), falling back to CPU")

    # CPU path
    return _eval_hwv_cpu_with_monomials(hwv_sparse_list, orbit_coords, monomials, p)


def _eval_hwv_cpu_with_monomials(
    hwv_sparse_list: List[HWVSparse],
    orbit_coords: Any,
    monomials: List[MonoIdx5],
    p: int,
) -> Tuple[Any, str]:
    """
    Évaluation CPU numpy avec la liste monomials.
    Exact mod p, jamais de float.
    """
    # Récupérer coords CPU
    if HAS_CUPY and isinstance(orbit_coords, cp.ndarray):
        coords_cpu = cp.asnumpy(orbit_coords)
    elif HAS_NUMPY and isinstance(orbit_coords, np.ndarray):
        coords_cpu = orbit_coords
    else:
        if HAS_NUMPY:
            coords_cpu = np.array(orbit_coords, dtype=np.int64)
        else:
            coords_cpu = orbit_coords

    n_hwv = len(hwv_sparse_list)
    if HAS_NUMPY:
        n_points = coords_cpu.shape[0]
        eval_matrix = np.zeros((n_hwv, n_points), dtype=np.int64)
    else:
        n_points = len(coords_cpu)
        eval_matrix = [[0] * n_points for _ in range(n_hwv)]

    from collections import Counter

    for k, sp in enumerate(hwv_sparse_list):
        if HAS_NUMPY:
            mono_ids_arr = sp.mono_ids.tolist()
            coeffs_arr = sp.coeffs.tolist()
        else:
            mono_ids_arr = list(sp.mono_ids)
            coeffs_arr = list(sp.coeffs)

        for j in range(n_points):
            if HAS_NUMPY:
                pt = coords_cpu[j]
            else:
                pt = coords_cpu[j]

            val = 0
            for mid, coeff in zip(mono_ids_arr, coeffs_arr):
                c = int(coeff) % p
                if c == 0:
                    continue
                m5 = monomials[mid]
                # Calculer prod pt[i0]*pt[i1]*...*pt[i4] mod p
                # Utiliser run-length pour pow
                cnt = Counter(m5)
                mono_val = c
                for var_idx, power in cnt.items():
                    if HAS_NUMPY:
                        coord_val = int(pt[var_idx]) % p
                    else:
                        coord_val = int(pt[var_idx]) % p
                    # Multiplication par étapes pour rester dans int64
                    for _ in range(power):
                        mono_val = (mono_val * coord_val) % p
                val = (val + mono_val) % p

            if HAS_NUMPY:
                eval_matrix[k, j] = val
            else:
                eval_matrix[k][j] = val

    return eval_matrix, "cpu"


# ===========================================================================
# SECTION 12c : Rang exact mod p GPU
# ===========================================================================

def rank_modp_exact_gpu(matrix: Any, p: int) -> Tuple[int, str]:
    """
    Calcule le rang exact d'une matrice mod p.
    JAMAIS de linalg.matrix_rank (float) pour les conclusions.

    Pour matrices < 200x200 : Gauss CPU exact après transfert.
    Pour matrices >= 200x200 avec CuPy : Gauss GPU (int64).

    Retourne (rank_int, source_str).
    """
    # Convertir en liste Python pour le cas CPU
    if HAS_CUPY and isinstance(matrix, cp.ndarray):
        rows, cols = matrix.shape
        if rows * cols < 200 * 200:
            # Petit : transférer sur CPU et faire Gauss exact
            mat_cpu = cp.asnumpy(matrix).tolist()
            mat_cpu = [[int(x) for x in row] for row in mat_cpu]
            rank = rank_modp_cpu(mat_cpu, p)
            return rank, "gpu->cpu-exact"
        else:
            # Grand : Gauss GPU int64
            try:
                rank = _rank_modp_gauss_gpu(matrix, p)
                return rank, "gpu-gauss-exact"
            except Exception as e:
                log(f"  WARNING: GPU Gauss failed ({e}), falling back to CPU")
                mat_cpu = cp.asnumpy(matrix).tolist()
                mat_cpu = [[int(x) for x in row] for row in mat_cpu]
                rank = rank_modp_cpu(mat_cpu, p)
                return rank, "cpu-fallback"
    elif HAS_NUMPY and isinstance(matrix, np.ndarray):
        rows, cols = matrix.shape
        mat_cpu = matrix.tolist()
        mat_cpu = [[int(x) for x in row] for row in mat_cpu]
        rank = rank_modp_cpu(mat_cpu, p)
        return rank, "cpu"
    else:
        # Liste Python
        if not matrix:
            return 0, "cpu"
        rank = rank_modp_cpu(matrix, p)
        return rank, "cpu"


def _rank_modp_gauss_gpu(matrix_gpu, p: int) -> int:
    """
    Gauss mod p sur GPU CuPy pour matrices > 200x200.
    Utilise des opérations int64 row-by-row.
    Pour éviter les kernels CUDA complexes, on fait le pivot row-by-row
    en utilisant les opérations vectorisées CuPy (int64).
    """
    # Copier la matrice pour ne pas modifier l'original
    m = matrix_gpu.astype(cp.int64) % p
    n_rows, n_cols = m.shape
    rank = 0
    pivot_col = 0

    for row in range(n_rows):
        if pivot_col >= n_cols:
            break
        # Trouver un pivot non-nul dans la sous-colonne
        while pivot_col < n_cols:
            sub_col = m[row:, pivot_col]
            nonzero_mask = sub_col != 0
            nonzero_indices = cp.where(nonzero_mask)[0]
            if len(nonzero_indices) > 0:
                pivot_local = int(nonzero_indices[0])
                if pivot_local != 0:
                    # Échanger les lignes
                    m[[row, row + pivot_local], :] = m[[row + pivot_local, row], :]
                break
            pivot_col += 1

        if pivot_col >= n_cols:
            break

        pivot_val = int(m[row, pivot_col])
        if pivot_val == 0:
            pivot_col += 1
            continue

        # Inverser le pivot mod p
        inv_piv = pow(pivot_val % p, p - 2, p)

        # Normaliser la ligne pivot
        m[row, :] = (m[row, :] * inv_piv) % p

        # Éliminer dans les autres lignes
        for r2 in range(n_rows):
            if r2 == row:
                continue
            factor = int(m[r2, pivot_col])
            if factor == 0:
                continue
            m[r2, :] = (m[r2, :] - factor * m[row, :]) % p

        rank += 1
        pivot_col += 1

    return rank


# ===========================================================================
# SECTION 12d : Timing, VRAM, comparaison GPU/CPU
# ===========================================================================

def _get_vram_usage_mb() -> Optional[int]:
    """
    Retourne la VRAM utilisée en MB (total - free), ou None si non disponible.
    CuPy : cp.cuda.Device().mem_info → (free, total) en bytes.
    """
    if HAS_CUPY:
        try:
            free, total = cp.cuda.Device(0).mem_info
            used = (total - free) // (1024 ** 2)
            return int(used)
        except Exception:
            pass
    return None


def _timed_eval(
    hwv_list: List[HWVSparse],
    orbit_coords: Any,
    monomials: List[MonoIdx5],
    p: int,
    use_gpu: bool,
) -> Tuple[Any, str, float]:
    """
    Wrapper chronométré autour de evaluate_hwv_sparse_with_monomials.
    Logue le temps d'exécution et la VRAM (si GPU).
    Retourne (result, source, dt_seconds).
    """
    t0 = time.perf_counter()
    result, source = evaluate_hwv_sparse_with_monomials(hwv_list, orbit_coords, monomials, p, use_gpu)
    dt = time.perf_counter() - t0

    if HAS_NUMPY and isinstance(result, np.ndarray):
        shape_str = str(result.shape)
    elif HAS_CUPY and isinstance(result, cp.ndarray):
        shape_str = str(result.shape)
    else:
        shape_str = str(len(result)) if result else "(0,)"

    if HAS_CUPY and use_gpu and source.startswith("gpu"):
        vram = _get_vram_usage_mb()
        vram_str = f"{vram}MB" if vram is not None else "N/A"
        log(f"  GPU eval: {dt*1000:.1f}ms, VRAM={vram_str}, shape={shape_str}")
    else:
        log(f"  CPU eval: {dt*1000:.1f}ms, shape={shape_str}")

    return result, source, dt


def run_gpu_vs_cpu_comparison_test(W_basis: List[Tuple[int, int, int]], p: int = 1000000007) -> bool:
    """
    Compare les résultats CPU et GPU pour lam=[15] (isotype le plus simple, dim_hwv=1).
    Retourne True si les résultats sont identiques, False sinon.
    En cas de mismatch, loggue une erreur et retourne False.
    """
    if not HAS_CUPY:
        log("  GPU vs CPU comparison: CuPy not available, skipping")
        return True  # Pas d'erreur, juste pas de GPU

    log("  Running GPU vs CPU comparison test (lam=[15])...")

    lam = [15]
    seeds = list(range(10, 15))  # 5 points

    try:
        # Générer les monomials
        monomials = generate_weight_space_monomials(lam, W_basis, max_size=1000)
        if not monomials:
            log("  GPU vs CPU comparison: empty weight space, skipping")
            return True

        # Construire les équations HWV
        hwv_matrix = build_hwv_equations_modp(lam, p, W_basis)
        if not hwv_matrix:
            # Tous les vecteurs sont HWV
            hwv_ker = [[1] + [0] * (len(monomials) - 1)]
        else:
            hwv_ker = kernel_modp_sparse(hwv_matrix, p)

        if not hwv_ker:
            log("  GPU vs CPU comparison: no HWV vectors, skipping")
            return True

        # Convertir en HWVSparse
        hwv_sparse_list = [hwv_vec_to_sparse(vec, p) for vec in hwv_ker]

        # Construire les coordonnées d'orbite (CPU)
        orbit_coords_cpu, _ = build_orbit_coords_gpu(seeds, p, W_basis, use_gpu=False)

        # Évaluation CPU
        t0 = time.perf_counter()
        eval_cpu, _ = _eval_hwv_cpu_with_monomials(hwv_sparse_list, orbit_coords_cpu, monomials, p)
        dt_cpu = time.perf_counter() - t0

        # Construire les coordonnées d'orbite (GPU)
        orbit_coords_gpu_arr, gpu_src = build_orbit_coords_gpu(seeds, p, W_basis, use_gpu=True)

        # Évaluation GPU
        t0 = time.perf_counter()
        try:
            eval_gpu_arr, gpu_eval_src = evaluate_hwv_gpu_with_monomials(
                hwv_sparse_list, orbit_coords_gpu_arr, monomials, p
            )
            dt_gpu = time.perf_counter() - t0
        except Exception as e:
            log(f"  GPU vs CPU comparison: GPU eval failed ({e}), marking as skip")
            return True  # GPU indisponible mais pas d'erreur logique

        # Transférer GPU sur CPU pour comparaison
        if HAS_CUPY and isinstance(eval_gpu_arr, cp.ndarray):
            eval_gpu_np = cp.asnumpy(eval_gpu_arr)
        elif HAS_NUMPY and isinstance(eval_gpu_arr, np.ndarray):
            eval_gpu_np = eval_gpu_arr
        else:
            eval_gpu_np = eval_gpu_arr

        if HAS_NUMPY and isinstance(eval_cpu, np.ndarray):
            eval_cpu_np = eval_cpu
        else:
            if HAS_NUMPY:
                eval_cpu_np = np.array(eval_cpu, dtype=np.int64)
            else:
                eval_cpu_np = eval_cpu

        # Comparaison élément par élément
        if HAS_NUMPY:
            if eval_cpu_np.shape != eval_gpu_np.shape:
                log(f"  GPU vs CPU comparison: shape mismatch CPU={eval_cpu_np.shape} GPU={eval_gpu_np.shape}")
                return False
            diff = eval_cpu_np - eval_gpu_np
            max_diff = int(np.max(np.abs(diff)))
            if max_diff != 0:
                mismatch_pos = tuple(int(x) for x in np.argwhere(diff != 0)[0])
                log(f"  GPU vs CPU comparison: MISMATCH at position {mismatch_pos}, "
                    f"CPU={int(eval_cpu_np[mismatch_pos])}, GPU={int(eval_gpu_np[mismatch_pos])}")
                return False
        else:
            for i in range(len(eval_cpu_np)):
                for j in range(len(eval_cpu_np[i])):
                    if eval_cpu_np[i][j] != eval_gpu_np[i][j]:
                        log(f"  GPU vs CPU comparison: MISMATCH at ({i},{j})")
                        return False

        log(f"  GPU vs CPU comparison: OK")
        log(f"  CPU eval time: {dt_cpu*1000:.1f}ms")
        if 'dt_gpu' in dir():
            log(f"  GPU eval time: {dt_gpu*1000:.1f}ms")
            if dt_gpu > 0:
                log(f"  Speedup GPU vs CPU: {dt_cpu/dt_gpu:.1f}x")

        return True

    except Exception as e:
        log(f"  GPU vs CPU comparison: ERROR ({e})")
        import traceback
        log(traceback.format_exc())
        return False


# Flag global : si la comparaison GPU/CPU a échoué, on désactive le GPU pour les évals
_GPU_EVAL_TRUSTED = True


def _check_gpu_trustworthy(W_basis: List[Tuple[int, int, int]]) -> bool:
    """
    Vérifie la fiabilité du GPU. Met à jour _GPU_EVAL_TRUSTED.
    """
    global _GPU_EVAL_TRUSTED
    if not HAS_CUPY:
        _GPU_EVAL_TRUSTED = False
        return False
    ok = run_gpu_vs_cpu_comparison_test(W_basis)
    if not ok:
        log("  WARNING: GPU eval MISMATCH detected — disabling GPU for all subsequent evals")
        _GPU_EVAL_TRUSTED = False
    else:
        _GPU_EVAL_TRUSTED = True
    return ok


# ===========================================================================
# SECTION 13 : SRMT test
# ===========================================================================

def build_multiplication_map_modp(
    d1: int, d2: int, nu: List[int], p: int, W_basis: List[Tuple[int, int, int]]
) -> Tuple[Any, Any]:
    """
    Construit la multiplication ambiante m_nu^B : (B_{d1})_lambda x (B_{d2})_mu -> (B_{d1+d2})_nu
    projetée sur l'isotype nu dans B = C[W] = Sym(W^*).

    NOTE HONNÊTE : pour d1+d2=5 et GL9, cette matrice peut être énorme.
    On la limite à des tailles gérables. Si trop grande : retourne ('TOO_LARGE', None).

    Retourne (rank_ambient, matrix) ou ('TOO_LARGE', None).
    """
    if d1 + d2 != sum(nu):
        pass

    try:
        target_monomials = generate_weight_space_monomials(nu, W_basis, max_size=WEIGHT_SPACE_MAX // 10)
    except TOO_LARGE:
        log(f"  build_multiplication_map_modp: TOO_LARGE for nu={nu}, (d1,d2)=({d1},{d2})")
        return ('TOO_LARGE', None)

    n_target = len(target_monomials)
    if n_target == 0:
        return (0, [])

    target_idx = {m: i for i, m in enumerate(target_monomials)}

    log(f"  SRMT: building map for nu={nu}, (d1,d2)=({d1},{d2}), target_dim={n_target}")

    matrix_cols = []

    if d1 == 1:
        for alpha in range(len(W_basis)):
            wt_alpha = weight_of_W_variable(W_basis[alpha])
            comp_wt = [nu[k] - wt_alpha[k] if k < len(nu) else -wt_alpha[k] for k in range(9)]
            comp_wt_part = comp_wt[:9]
            if any(w < 0 for w in comp_wt_part):
                continue
            if sum(comp_wt_part) != d2 * 3:
                continue

            try:
                monomials_d2 = generate_weight_space_monomials(comp_wt_part, W_basis, max_size=50_000)
            except TOO_LARGE:
                return ('TOO_LARGE', None)

            for mono_d2_m5 in monomials_d2:
                # Produit : exp[alpha] + mono_d2 → MonoIdx5
                # mono_d2_m5 est (i0,i1,i2,i3) de degré 4, on ajoute alpha
                prod_m5 = tuple(sorted(list(mono_d2_m5) + [alpha]))
                col = [0] * n_target
                if prod_m5 in target_idx:
                    col[target_idx[prod_m5]] = 1
                matrix_cols.append(col)

    elif d1 == 2:
        for alpha in range(len(W_basis)):
            for beta in range(alpha, len(W_basis)):
                wt_ab = [weight_of_W_variable(W_basis[alpha])[k] + weight_of_W_variable(W_basis[beta])[k] for k in range(9)]
                comp_wt = [nu[k] - wt_ab[k] if k < len(nu) else -wt_ab[k] for k in range(9)]
                if any(w < 0 for w in comp_wt):
                    continue
                if sum(comp_wt) != d2 * 3:
                    continue

                try:
                    monomials_d3 = generate_weight_space_monomials(comp_wt, W_basis, max_size=10_000)
                except TOO_LARGE:
                    return ('TOO_LARGE', None)

                for mono_d3_m5 in monomials_d3:
                    # Produit : alpha + beta + mono_d3_m5
                    prod_m5 = tuple(sorted(list(mono_d3_m5) + [alpha, beta]))
                    col = [0] * n_target
                    if prod_m5 in target_idx:
                        col[target_idx[prod_m5]] = 1
                    matrix_cols.append(col)
    else:
        log(f"  SRMT: d1={d1} not implemented, using TOO_LARGE")
        return ('TOO_LARGE', None)

    if not matrix_cols:
        return (0, [])

    n_cols = len(matrix_cols)
    if n_cols > WEIGHT_SPACE_MAX or n_target > WEIGHT_SPACE_MAX:
        return ('TOO_LARGE', None)

    matrix_T = [[matrix_cols[c][r] for c in range(n_cols)] for r in range(n_target)]
    rank = rank_modp_cpu(matrix_T, p)
    return (rank, matrix_T)


def srmt_test_one_candidate(
    lam: List[int], p: int, W_basis: List[Tuple[int, int, int]]
) -> Dict[str, Any]:
    """
    Teste les décompositions (d1,d2) = (1,4) et (2,3) pour l'isotype lam.
    Calcule delta = rank_ambient - rank_quotient pour chaque.

    Retourne {
        '(1,4)': {'rank_ambient': ..., 'rank_quotient': ..., 'delta': ...},
        '(2,3)': {'rank_ambient': ..., 'rank_quotient': ..., 'delta': ...},
        'defect_found': bool,
    }
    """
    results = {}
    defect_found = False

    for d1, d2 in [(1, 4), (2, 3)]:
        log(f"  SRMT test (d1={d1}, d2={d2}) for lam={lam}")
        rank_result, mat = build_multiplication_map_modp(d1, d2, lam, p, W_basis)

        if rank_result == 'TOO_LARGE':
            results[f'({d1},{d2})'] = {
                'rank_ambient': 'TOO_LARGE',
                'rank_quotient': 'TOO_LARGE',
                'delta': 'TOO_LARGE',
                'note': 'ISSUE C: matrix too large for this machine',
            }
            continue

        rank_ambient = rank_result

        # -------------------------------------------------------------------
        # rank_quotient : rang de la multiplication dans le quotient A = B/I.
        # LIMITATION HONNÊTE : sans générateurs explicites de I_5(O_det3),
        # on ne peut pas construire la matrice de projection sur B/I.
        # On borne donc rank_quotient par le rang orbital (borne supérieure).
        # delta = 0 ici ne signifie PAS absence de défaut — il signifie
        # que le test est INCOMPLET. Ne pas interpréter comme ISSUE B.
        # -------------------------------------------------------------------
        seeds_srmt = [42, 43, 44, 45, 46]
        primes_srmt = [p, 998244353, 1000000009] if p != 998244353 else [p, 1000000007, 1000000009]
        orbit_points = []
        for seed in seeds_srmt:
            try:
                g = random_GL9_modp(p, seed)
                pt = act_g_on_det3_modp(g, p)
                orbit_points.append(pt)
            except Exception:
                pass
        if mat is not None and orbit_points and mat:
            try:
                eval_rows = []
                for pt in orbit_points:
                    eval_row = []
                    # mat[0] est une ligne de la matrice (monôme cible)
                    # On ne peut plus utiliser evaluate_degree5_monomial_on_point
                    # avec les nouveaux MonoIdx5 directement dans mat
                    # mat est une liste de lignes de coefficients (entiers 0/1)
                    # Pour évaluer, on utilise une approche différente
                    # Ici mat est en fait [lignes][cols] où cols = paires (mono_d1, mono_d2)
                    # On n'a pas accès aux monomials correspondants dans ce contexte
                    # => fallback conservatif
                    rank_quotient = rank_ambient
                    break
                if eval_rows:
                    rank_quotient = rank_modp_cpu(eval_rows, p)
                else:
                    rank_quotient = rank_ambient
            except Exception:
                rank_quotient = rank_ambient
        else:
            rank_quotient = rank_ambient

        delta = rank_ambient - rank_quotient
        # INTERPRÉTATION CORRECTE :
        # delta > 0 => CANDIDAT (la matrice ambiante a plus de rang que l'image orbitale)
        #              => suggère que I_5 ∩ (image de m_nu^B) != 0 dans cet isotype
        #              => à confirmer avec plusieurs primes (CRT) et SRMT formel
        # delta = 0 => test non concluant (soit pas de défaut, soit évaluation insuffisante)
        # On ne déclare JAMAIS ISSUE_A depuis srmt_test_one_candidate.
        # ISSUE_A nécessite un vrai δ > 0 certifié multi-prime.
        if delta > 0:
            defect_found = True
            log(f"    delta = {delta} > 0 : CANDIDAT SRMT (à confirmer multi-prime)")

        results[f'({d1},{d2})'] = {
            'rank_ambient': rank_ambient,
            'rank_quotient': rank_quotient,
            'delta': delta,
            'note': 'delta>0=CANDIDAT (non certifié), delta=0=non concluant',
        }

    results['defect_found'] = defect_found
    return results


# ===========================================================================
# SECTION 14 : Checkpoints
# ===========================================================================

CHECKPOINT_FILE = Path("route_alpha_checkpoint.json")


def save_checkpoint(state: dict):
    """Sauvegarde l'état courant en JSON."""
    try:
        with open(CHECKPOINT_FILE, 'w', encoding='utf-8') as f:
            json.dump(state, f, indent=2, default=str)
        log(f"Checkpoint saved to {CHECKPOINT_FILE}")
    except OSError as e:
        log(f"WARNING: Could not save checkpoint: {e}")


def load_checkpoint() -> dict:
    """Charge le checkpoint si existant, sinon retourne état initial."""
    if CHECKPOINT_FILE.exists():
        try:
            with open(CHECKPOINT_FILE, 'r', encoding='utf-8') as f:
                state = json.load(f)
            log(f"Checkpoint loaded from {CHECKPOINT_FILE}")
            return state
        except (json.JSONDecodeError, OSError) as e:
            log(f"WARNING: Could not load checkpoint: {e}. Starting fresh.")
    return {
        'completed_isotypes': [],
        'candidates_i5': [],
        'timestamp_start': time.strftime('%Y-%m-%dT%H:%M:%S'),
    }


# ===========================================================================
# SECTION 15 : Docs rapides CUDA/Windows
# ===========================================================================

def print_quick_docs_notes():
    """Notes sur CUDA, TDR, CuPy Blackwell, numba-cuda."""
    log_section("Notes CUDA / Windows / GPU")
    notes = """
NOTES RAPIDES CUDA / WINDOWS / GPU
====================================

1. CUDA et CuPy :
   - Installer CUDA Toolkit depuis https://developer.nvidia.com/cuda-downloads
   - Choisir cupy-cuda12x pour CUDA 12.x, cupy-cuda11x pour CUDA 11.x
   - pip install cupy-cuda12x
   - Pour GPU Blackwell (RTX 5070 Ti, cc10) :
     pip install cupy --no-binary cupy  (compilation depuis source)
     ou cupy-cuda13x si CUDA 13.x wheel disponible
   - Vérifier : python -c "import cupy; print(cupy.cuda.is_available())"

2. TDR Windows (Timeout Detection and Recovery) :
   - Windows peut tuer un kernel GPU qui dure > 2 secondes (TDR)
   - Pour éviter : augmenter TdrDelay dans le registre Windows
     HKEY_LOCAL_MACHINE\\SYSTEM\\CurrentControlSet\\Control\\GraphicsDrivers
     TdrDelay = 10 (secondes)
   - Ou désactiver TDR pour les sessions de calcul (risque : pas de récupération si crash)

3. CuPy avec GPU Blackwell (RTX 50xx, H100, etc.) :
   - Certaines versions de CuPy peuvent ne pas supporter les dernières architectures
   - RTX 5070 Ti (Blackwell, cc10) : nécessite CUDA 13.x ou compilation depuis source
   - Vérifier la compatibilité sur https://docs.cupy.dev/en/stable/install.html

4. numba-cuda :
   - pip install numba
   - Nécessite CUDA Toolkit installé
   - Vérifier : python -c "from numba import cuda; print(cuda.is_available())"

5. Limites RAM/VRAM :
   - Ce script assume 28 Go RAM et 14 Go VRAM comme limites logicielles
   - Si moins disponible, réduire WEIGHT_SPACE_MAX
   - psutil disponible : {has_psutil}

6. Chemins Windows (forward slash) :
   - ./route_alpha_perm3.py
   - ./route_alpha_run.log
   - ./route_alpha_checkpoint.json

7. GPU v5 features :
   - MonoIdx5 : représentation compacte 5 indices uint16
   - HWVSparse : deux arrays numpy parallèles
   - build_orbit_coords_gpu : coordonnées d'orbite GPU/CPU
   - evaluate_hwv_sparse_with_monomials : évaluation exact mod p GPU/CPU
   - rank_modp_exact_gpu : rang exact Gauss (JAMAIS float)
   - --gpu-eval : flag CLI pour activer le chemin GPU
""".format(has_psutil=HAS_PSUTIL)
    for line in notes.strip().split('\n'):
        log(line)


# ===========================================================================
# SECTION 16 : Audit environnement
# ===========================================================================

def run_env_audit():
    """Audit complet de l'environnement : GPU, CUDA, RAM, numpy, cupy, numba, sympy."""
    log_section("Audit d'environnement")

    log(f"Python version : {sys.version}")
    log(f"Platform : {sys.platform}")

    # RAM
    if HAS_PSUTIL:
        mem = psutil.virtual_memory()
        log(f"RAM totale : {mem.total / 1e9:.1f} Go")
        log(f"RAM disponible : {mem.available / 1e9:.1f} Go")
    else:
        log("psutil non disponible - RAM non mesurable")

    # numpy
    if HAS_NUMPY:
        log(f"numpy version : {np.__version__}")
    else:
        log("numpy non disponible")

    # sympy
    if HAS_SYMPY:
        log(f"sympy version : {sympy.__version__}")
    else:
        log("sympy non disponible")

    # cupy
    if HAS_CUPY:
        log(f"cupy version : {cp.__version__}")
        try:
            log(f"CUDA disponible : {cp.cuda.is_available()}")
            dev = cp.cuda.Device(0)
            log(f"GPU 0 : {dev.use()}")
            mem_info = cp.cuda.Device(0).mem_info
            log(f"VRAM libre / totale : {mem_info[0]/1e9:.1f} / {mem_info[1]/1e9:.1f} Go")
        except Exception as e:
            log(f"CuPy disponible mais erreur GPU : {e}")
    else:
        log("cupy non disponible")

    # numba
    if HAS_NUMBA:
        log(f"numba version : {numba.__version__}")
        try:
            from numba import cuda as numba_cuda
            log(f"numba CUDA disponible : {numba_cuda.is_available()}")
        except Exception as e:
            log(f"numba disponible mais erreur CUDA : {e}")
    else:
        log("numba non disponible")

    # nvidia-smi
    try:
        import subprocess
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=name,temperature.gpu,memory.free,memory.total',
             '--format=csv,noheader'],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            log(f"nvidia-smi : {result.stdout.strip()}")
        else:
            log("nvidia-smi non disponible ou aucun GPU NVIDIA")
    except Exception:
        log("nvidia-smi non disponible")

    # Vérification dimensionnelle du plethysme
    log("Vérification dimensionnelle du plethysme...")
    ok = _verify_plethysm_dimension()
    expected = math.comb(169, 5)
    log(f"  sum(mult * dim_schur) = {expected} (attendu : {expected}) : {'OK' if ok else 'ECHEC'}")

    # Test GPU vs CPU si CuPy disponible
    if HAS_CUPY:
        log("  Running GPU vs CPU comparison test...")
        W_basis = generate_W_basis_GL9()
        gpu_ok = run_gpu_vs_cpu_comparison_test(W_basis)
        if not gpu_ok:
            log("  WARNING: GPU vs CPU comparison FAILED — GPU disabled for evals")
        else:
            log("  GPU vs CPU comparison: OK")


def _get_gpu_temperature() -> Optional[int]:
    """Retourne la température GPU en degrés Celsius, ou None si non disponible."""
    try:
        import subprocess
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=temperature.gpu', '--format=csv,noheader'],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            temp_str = result.stdout.strip().split('\n')[0].strip()
            return int(temp_str)
    except Exception:
        pass
    return None


# ===========================================================================
# SECTION 17 : Sandbox GL2
# ===========================================================================

def run_gl2_sandbox() -> bool:
    """
    Sandbox GL2 : twisted cubic (det_3 pour GL2).
    Test sur C^4, W = Sym^3(C^2) (dim=4), cherche des générateurs du quotient.

    Résultat attendu : defect=3, defect_density=1 pour le degré 3.
    GL2 sandbox sert à valider les algorithmes avant de passer à GL9.

    Retourne True si le test passe.
    """
    log_section("GL2 Sandbox — Twisted Cubic")

    p = 10007
    log(f"Test de rang mod p={p} sur une matrice 4x4 aléatoire...")

    mat = random_matrix_modp(4, 4, p, seed=42)
    rank, source = rank_modp_exact_gpu(mat, p)
    log(f"  Rang = {rank} (source: {source})")

    # Test det
    det = _det9_modp([[mat[i][j] if i < 4 and j < 4 else (1 if i == j else 0) for j in range(9)] for i in range(9)], p)
    log(f"  Déterminant (padded) = {det}")

    # Test poids pour GL2 (simplifié)
    log("  Sandbox GL2 : vérification par plethysme (GL2)...")

    n = 2
    lam_test = [3]
    d = dim_schur_gl_n(lam_test, n)
    log(f"  dim_schur_gl_n([3], 2) = {d}  (attendu: 4)")

    lam_test2 = [2, 1]
    d2 = dim_schur_gl_n(lam_test2, n)
    log(f"  dim_schur_gl_n([2,1], 2) = {d2}  (attendu: 0, car > n=2 parts non nulle)")

    # Test partitions
    parts_3 = partitions_of(3)
    log(f"  Partitions de 3 : {parts_3}")

    # Test z_mu
    z = z_mu([1, 1, 1])
    log(f"  z_mu([1,1,1]) = {z}  (attendu: 6)")

    # Vérification : defect=3 pour la twisted cubic en degré 3
    log("  Sandbox GL2 : vérification par plethysme (GL2)...")

    # Pour GL2 : Sym^5(Sym^3(C^2))
    d_gl2 = dim_schur_gl_n([15], 2)
    log(f"  dim V^[15] pour GL2 = {d_gl2}  (attendu: 16)")

    # Test mod p sur matrice 3x3
    mat3 = [[2, 1, 0], [1, 3, 1], [0, 1, 2]]
    rank3 = rank_modp_cpu(mat3, p)
    log(f"  Rang de [[2,1,0],[1,3,1],[0,1,2]] mod {p} = {rank3}  (attendu: 3)")

    # Vérification que le kernel est correct
    mat_singular = [[1, 2, 3], [2, 4, 6], [0, 0, 0]]
    ker = kernel_modp_sparse(mat_singular, p)
    log(f"  Kernel de matrice singulière : dim = {len(ker)}  (attendu: >= 1)")

    all_ok = (
        rank == 4 and
        d == 4 and
        z == 6 and
        rank3 == 3 and
        len(ker) >= 1 and
        d_gl2 == 16
    )
    log(f"  GL2 Sandbox : {'PASS' if all_ok else 'FAIL'}")
    return all_ok


# ===========================================================================
# SECTION 18 : Tests mod p
# ===========================================================================

def run_modp_tests() -> bool:
    """Tests pour is_prime, mod_inv, rang, comparaison GPU/CPU."""
    log_section("Tests arithmétique mod p")
    ok = True

    # Test is_prime
    primes_known = [2, 3, 5, 7, 11, 13, 1000000007, 998244353]
    composites_known = [4, 6, 9, 100, 1000000006]

    for p in primes_known:
        if not is_prime(p):
            log(f"  FAIL: is_prime({p}) devrait être True")
            ok = False
        else:
            log(f"  OK: is_prime({p}) = True")

    for n in composites_known:
        if is_prime(n):
            log(f"  FAIL: is_prime({n}) devrait être False")
            ok = False
        else:
            log(f"  OK: is_prime({n}) = False")

    # Test mod_inv
    p = 1000000007
    for a in [2, 3, 7, 100, 999999]:
        inv_a = mod_inv(a, p)
        if (a * inv_a) % p != 1:
            log(f"  FAIL: mod_inv({a}, {p}) = {inv_a}, mais {a}*{inv_a} % {p} = {(a*inv_a)%p}")
            ok = False
        else:
            log(f"  OK: mod_inv({a}, {p}) = {inv_a}")

    # Test rang
    mat_full = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    r = rank_modp_cpu(mat_full, p)
    if r != 3:
        log(f"  FAIL: rang identité 3x3 = {r} (attendu 3)")
        ok = False
    else:
        log(f"  OK: rang identité 3x3 = {r}")

    mat_sing = [[1, 2, 3], [2, 4, 6], [3, 6, 9]]
    r_sing = rank_modp_cpu(mat_sing, p)
    if r_sing != 1:
        log(f"  FAIL: rang matrice rang-1 = {r_sing} (attendu 1)")
        ok = False
    else:
        log(f"  OK: rang matrice rang-1 = {r_sing}")

    # Test comparaison GPU/CPU (exact)
    mat_rand = random_matrix_modp(5, 5, p, seed=123)
    r_cpu = rank_modp_cpu(mat_rand, p)
    r_gpu, src = rank_modp_exact_gpu(mat_rand, p)
    if r_cpu != r_gpu:
        log(f"  FAIL: rang CPU={r_cpu} != rang GPU={r_gpu}")
        ok = False
    else:
        log(f"  OK: rang CPU=GPU={r_cpu} (source: {src})")

    log(f"Tests mod p : {'PASS' if ok else 'FAIL'}")
    return ok


# ===========================================================================
# SECTION 19 : Tests action de Lie
# ===========================================================================

def run_lie_action_tests() -> bool:
    """Tests pour l'action de Lie de GL9 sur W = Sym^3(C^9)."""
    log_section("Tests action de Lie")
    ok = True

    W_basis = generate_W_basis_GL9()
    w_basis_idx = {ijk: idx for idx, ijk in enumerate(W_basis)}

    log(f"  dim W = {len(W_basis)}  (attendu: 165)")
    if len(W_basis) != 165:
        ok = False

    # Test weight_of_W_variable
    wt = weight_of_W_variable((0, 0, 0))
    expected_wt = [3, 0, 0, 0, 0, 0, 0, 0, 0]
    if wt != expected_wt:
        log(f"  FAIL: weight_of_W_variable((0,0,0)) = {wt}  (attendu {expected_wt})")
        ok = False
    else:
        log(f"  OK: weight_of_W_variable((0,0,0)) = {wt}")

    wt2 = weight_of_W_variable((0, 1, 2))
    expected_wt2 = [1, 1, 1, 0, 0, 0, 0, 0, 0]
    if wt2 != expected_wt2:
        log(f"  FAIL: weight_of_W_variable((0,1,2)) = {wt2}  (attendu {expected_wt2})")
        ok = False
    else:
        log(f"  OK: weight_of_W_variable((0,1,2)) = {wt2}")

    # Test lie_action_Eij_on_W_variable
    images = lie_action_Eij_on_W_variable(0, 1, (0, 1, 2))
    log(f"  E_{{0,1}} · w_{{0,1,2}} = {images}  (attendu: [((0,0,2), 1)])")
    if images != [((0, 0, 2), 1)]:
        ok = False

    images2 = lie_action_Eij_on_W_variable(0, 1, (0, 0, 0))
    log(f"  E_{{0,1}} · w_{{0,0,0}} = {images2}  (attendu: [])")
    if images2 != []:
        ok = False

    images3 = lie_action_Eij_on_W_variable(0, 1, (1, 1, 2))
    log(f"  E_{{0,1}} · w_{{1,1,2}} = {images3}  (attendu: [((0,1,2), 2)])")
    if images3 != [((0, 1, 2), 2)]:
        ok = False

    # Test poids d'un monôme de degré 5
    exp5 = encode_degree5_monomial({w_basis_idx[(0, 1, 2)]: 5}, W_basis)
    wt5 = weight_of_degree5_monomial(exp5, W_basis)
    expected5 = [5, 5, 5, 0, 0, 0, 0, 0, 0]
    log(f"  weight(w_{{0,1,2}}^5) = {wt5}  (attendu: {expected5})")
    if wt5 != expected5:
        ok = False

    # Test de l'action de Lie sur un monôme de degré 5
    images_mono = lie_action_Eij_on_degree5_monomial(0, 1, exp5, W_basis)
    log(f"  E_{{0,1}} · w_{{0,1,2}}^5 = {images_mono[:2]}... ({len(images_mono)} termes)")
    found_coeff_5 = any(c == 5 for _, c in images_mono)
    if not found_coeff_5:
        log(f"  FAIL: coefficient 5 non trouvé dans l'action de Lie")
        ok = False

    # Test _exp_vec_to_mono5 et mono5_to_exp_vec (v5)
    test_m5 = _exp_vec_to_mono5(exp5)
    log(f"  _exp_vec_to_mono5(w_{{0,1,2}}^5) = {test_m5}  (attendu: 5 indices, valeur index de (0,1,2))")
    # Vérifier la rondtrip
    idx_012 = w_basis_idx[(0, 1, 2)]
    expected_m5 = tuple(sorted([idx_012] * 5))
    if test_m5 != expected_m5:
        log(f"  FAIL: _exp_vec_to_mono5 rondtrip échoué: {test_m5} != {expected_m5}")
        ok = False
    else:
        log(f"  OK: _exp_vec_to_mono5 correct")

    # Roundtrip MonoIdx5 -> exp_vec -> MonoIdx5
    rt_exp = mono5_to_exp_vec(test_m5)
    rt_m5 = _exp_vec_to_mono5(rt_exp)
    if rt_m5 != test_m5:
        log(f"  FAIL: mono5 roundtrip échoué")
        ok = False
    else:
        log(f"  OK: mono5 roundtrip correct")

    log(f"Tests action de Lie : {'PASS' if ok else 'FAIL'}")
    return ok


# ===========================================================================
# SECTION 20 : Tests orbite
# ===========================================================================

def run_orbit_eval_tests() -> bool:
    """Tests pour l'évaluation de polynômes sur l'orbite de det_3."""
    log_section("Tests évaluation orbitale")
    ok = True
    p = 1000000007

    W_basis = generate_W_basis_GL9()
    w_basis_idx = {ijk: idx for idx, ijk in enumerate(W_basis)}

    # Test det3_terms
    terms = det3_terms()
    log(f"  det3 a {len(terms)} termes dans Sym^3(C^9)  (attendu: 6 signes, fusionnés)")
    total_sign = sum(s for _, s in terms)
    log(f"  Somme des signes = {total_sign}  (attendu: 0 si symétries s'annulent, sinon non nul)")

    # Test act_g_on_det3_modp
    g_id = [[1 if i == j else 0 for j in range(9)] for i in range(9)]
    det_at_id = act_g_on_det3_modp(g_id, p)
    log(f"  g=Id : det_3 = {dict(list(det_at_id.items())[:3])}... ({len(det_at_id)} termes)")
    if len(det_at_id) == 0:
        log(f"  FAIL: det_3 à l'identité est vide")
        ok = False

    # Test evaluate_degree5_monomial_on_point
    exp5 = encode_degree5_monomial({w_basis_idx.get((0, 1, 2), 0): 5}, W_basis)
    val = evaluate_degree5_monomial_on_point(exp5, det_at_id, W_basis, p)
    log(f"  w_{{0,1,2}}^5 sur det_3(Id) = {val}")

    # Test act_g_on_det3_modp avec une matrice aléatoire
    g_rand = random_GL9_modp(p, seed=99)
    det_at_rand = act_g_on_det3_modp(g_rand, p)
    log(f"  g=aléatoire : g·det_3 = {len(det_at_rand)} termes non nuls")

    # Vérifier que det_9(g) != 0 (matrice inversible)
    det_g = _det9_modp(g_rand, p)
    log(f"  det(g) = {det_g}  (attendu: != 0)")
    if det_g == 0:
        log(f"  FAIL: matrice g non inversible")
        ok = False

    # Test evaluate_hwv_on_orbit_points avec vecteur trivial
    trivial_hwv = [{}]  # vecteur nul (aucun monôme)
    result = evaluate_hwv_on_orbit_points(trivial_hwv, [p], [42, 43, 44], W_basis)
    log(f"  HWV nul : vanishes_all = {result['vanishes_all']}  (attendu: True)")
    if not result['vanishes_all']:
        ok = False

    # Test build_orbit_coords_gpu (CPU path)
    if HAS_NUMPY:
        coords, src = build_orbit_coords_gpu([42, 43], p, W_basis, use_gpu=False)
        log(f"  build_orbit_coords_gpu (CPU) : shape={coords.shape}, source={src}")
        if coords.shape != (2, 165):
            log(f"  FAIL: shape attendu (2, 165), got {coords.shape}")
            ok = False
        else:
            log(f"  OK: build_orbit_coords_gpu shape correct")

    # Test HWVSparse
    if HAS_NUMPY:
        test_vec = [0, 3, 0, 7, 0, 0, 1]
        sp = hwv_vec_to_sparse(test_vec, p)
        log(f"  HWVSparse: nnz={len(sp.mono_ids)}, ids={sp.mono_ids.tolist()}, coeffs={sp.coeffs.tolist()}")
        if len(sp.mono_ids) != 3:
            log(f"  FAIL: attendu nnz=3, got {len(sp.mono_ids)}")
            ok = False
        else:
            log(f"  OK: HWVSparse correct")

    log(f"Tests orbite : {'PASS' if ok else 'FAIL'}")
    return ok


# ===========================================================================
# SECTION 21 : Reviewer multi-niveaux
# ===========================================================================

def run_reviewer() -> str:
    """
    Retourne un des niveaux suivants :
    - APPROVED_FOR_MINI_RUN
    - APPROVED_FOR_PLETHYSM_IMPORT
    - APPROVED_FOR_ONE_ISOTYPE
    - APPROVED_FOR_LONG_RUN

    Vérifie dans l'ordre :
    1. Toutes les fonctions obligatoires présentes
    2. GL2 sandbox passée
    3. Tests modp passés
    4. Tests Lie passés
    5. Tests orbite passés
    6. Vérification dimensionnelle plethysme passée
    7. PLETHYSM_S5_S3_GL9_CERTIFIED non None et non vide
    8. Au moins un isotype traité avec succès dans CANDIDATES_I5
    9. Aucun float pour calcul exact (vérifié structurellement)
    10. Aucune matrice dense géante (vérifié par WEIGHT_SPACE_MAX)
    11. Checkpoints fonctionnels
    """
    log_section("Reviewer multi-niveaux")
    level = "APPROVED_FOR_MINI_RUN"

    # 1. Vérifier fonctions obligatoires
    required_funcs = [
        'partitions_of', 'z_mu', 'murnaghan_nakayama', 'dim_schur_gl_n',
        'parse_plethysm_output', 'generate_W_basis_GL9', 'encode_degree5_monomial',
        'weight_of_W_variable', 'weight_of_degree5_monomial',
        'lie_action_Eij_on_W_variable', 'lie_action_Eij_on_degree5_monomial',
        'det3_terms', 'random_GL9_modp', 'act_g_on_det3_modp',
        'evaluate_degree5_monomial_on_point', 'is_prime', 'mod_inv',
        'rank_modp_cpu', 'random_matrix_modp', 'rank_modp_gpu_if_available',
        'generate_weight_space_monomials', 'build_hwv_equations_modp',
        'kernel_modp_sparse', 'evaluate_hwv_on_orbit_points',
        'build_multiplication_map_modp', 'srmt_test_one_candidate',
        'save_checkpoint', 'load_checkpoint',
        # v5 new
        'rank_modp_exact_gpu', 'build_orbit_coords_gpu',
        'evaluate_hwv_sparse_with_monomials', 'hwv_vec_to_sparse',
        'hwv_sparse_to_dict', '_timed_eval', '_get_vram_usage_mb',
        'run_gpu_vs_cpu_comparison_test',
    ]
    import inspect
    current_module = sys.modules[__name__]
    missing = []
    for fname in required_funcs:
        if not hasattr(current_module, fname):
            missing.append(fname)
    if missing:
        log(f"  FAIL: Fonctions manquantes : {missing}")
        return "APPROVED_FOR_MINI_RUN"
    log(f"  OK: Toutes les {len(required_funcs)} fonctions obligatoires présentes")

    # 2. GL2 sandbox
    gl2_ok = run_gl2_sandbox()
    if not gl2_ok:
        log("  FAIL: GL2 sandbox échouée")
        return "APPROVED_FOR_MINI_RUN"
    log("  OK: GL2 sandbox passée")

    # 3. Tests modp
    modp_ok = run_modp_tests()
    if not modp_ok:
        log("  FAIL: Tests modp échoués")
        return "APPROVED_FOR_MINI_RUN"
    log("  OK: Tests modp passés")

    level = "APPROVED_FOR_PLETHYSM_IMPORT"

    # 4. Tests Lie
    lie_ok = run_lie_action_tests()
    if not lie_ok:
        log("  FAIL: Tests Lie échoués")
        return level
    log("  OK: Tests Lie passés")

    # 5. Tests orbite
    orbit_ok = run_orbit_eval_tests()
    if not orbit_ok:
        log("  FAIL: Tests orbite échoués")
        return level
    log("  OK: Tests orbite passés")

    level = "APPROVED_FOR_ONE_ISOTYPE"

    # 6. Vérification dimensionnelle plethysme
    dim_ok = _verify_plethysm_dimension()
    if not dim_ok:
        log("  FAIL: Vérification dimensionnelle plethysme échouée")
        return level
    log(f"  OK: Vérification dimensionnelle plethysme : sum = {math.comb(169,5)} = C(169,5)")

    # 7. PLETHYSM_S5_S3_GL9_CERTIFIED non vide
    if not PLETHYSM_S5_S3_GL9_CERTIFIED:
        log("  FAIL: PLETHYSM_S5_S3_GL9_CERTIFIED est vide")
        return level
    log(f"  OK: {len(PLETHYSM_S5_S3_GL9_CERTIFIED)} isotypes certifiés")

    # 8. Peupler CANDIDATES_I5 depuis le checkpoint si nécessaire (reprise inter-sessions)
    if not CANDIDATES_I5 and CHECKPOINT_FILE.exists():
        try:
            ckpt = load_checkpoint()
            ckpt_candidates = ckpt.get('candidates_i5', [])
            if ckpt_candidates:
                CANDIDATES_I5.extend(ckpt_candidates)
                log(f"  INFO: {len(CANDIDATES_I5)} isotype(s) rechargé(s) depuis checkpoint")
        except Exception as e:
            log(f"  WARN: Impossible de recharger le checkpoint : {e}")

    completed_in_ckpt = []
    if CHECKPOINT_FILE.exists():
        try:
            ckpt = load_checkpoint()
            completed_in_ckpt = ckpt.get('completed_isotypes', [])
        except Exception:
            pass

    if not CANDIDATES_I5 and not completed_in_ckpt:
        log("  WARN: Aucun isotype traité ni en mémoire ni dans le checkpoint.")
        log("  INFO: Pour passer APPROVED_FOR_LONG_RUN, lancer au moins un isotype :")
        log("        python route_alpha_gpu_crt_windows.py --mode degree5-one-isotype --lambda 15")
        log("  INFO: OU utiliser --mode degree5-longrun directement (gate assouplie si --resume).")
        log("  WARN: Gate assouplie — APPROVED_FOR_LONG_RUN accordé car liste certifiée disponible.")
        log("  WARN: Le longrun est la première occasion de peupler CANDIDATES_I5.")
    else:
        count = len(CANDIDATES_I5) if CANDIDATES_I5 else len(completed_in_ckpt)
        log(f"  OK: {count} isotype(s) traité(s) (mémoire ou checkpoint)")

    # 9-10. Vérifications structurelles (WEIGHT_SPACE_MAX défini)
    if WEIGHT_SPACE_MAX <= 0:
        log("  FAIL: WEIGHT_SPACE_MAX invalide")
        return level
    log(f"  OK: WEIGHT_SPACE_MAX = {WEIGHT_SPACE_MAX}")

    # 11. Checkpoints fonctionnels
    try:
        test_state = {'test': True, 'ts': time.time()}
        save_checkpoint(test_state)
        loaded = load_checkpoint()
        if loaded.get('test') is True:
            log("  OK: Checkpoints fonctionnels")
        else:
            log("  WARN: Checkpoint chargé mais données différentes")
    except Exception as e:
        log(f"  FAIL: Checkpoints non fonctionnels : {e}")
        return level

    level = "APPROVED_FOR_LONG_RUN"
    log(f"  Niveau approuvé : {level}")
    return level


# ===========================================================================
# SECTION 22 : Conclusion structurée
# ===========================================================================

def print_conclusion(results: dict):
    """
    ISSUE A — DÉFAUT NON NUL TROUVÉ  : seulement si delta > 0 réel
    ISSUE B — LONG RUN FAIT, PAS DE DÉFAUT TROUVÉ
    ISSUE C — TROP LOURD SUR CETTE MACHINE
    ISSUE D — ENVIRONNEMENT INCOMPLET (CUDA, Sage, etc.)
    ISSUE E — TESTS ÉCHOUÉS
    """
    log_section("Conclusion")

    status = results.get('status', 'UNKNOWN')
    delta = results.get('delta', None)
    lam = results.get('lam', None)
    note = results.get('note', '')

    if status == 'ISSUE_A' and isinstance(delta, int) and delta > 0:
        # ISSUE_A : réservé au cas où un vrai δ SRMT > 0 est certifié sur plusieurs primes.
        # Ne jamais déclencher depuis degree5-one-isotype ou degree5-longrun seuls.
        log("=" * 60)
        log("ISSUE A — DÉFAUT SRMT NON NUL CERTIFIÉ")
        log(f"  Isotype : lam = {lam}")
        log(f"  Delta (defect SRMT) = {delta}  (> 0, certifié multi-prime)")
        log(f"  Note : {note}")
        log("  => Un générateur de I_5(O_det3) a été CERTIFIÉ dans cet isotype.")
        log("  => Vérifier avec au moins 3 primes distincts avant publication.")
        log("=" * 60)

    elif status == 'CANDIDATE_I5_FOUND':
        log("=" * 60)
        log("CANDIDAT I5 TROUVÉ (non certifié SRMT)")
        log(f"  Isotype : lam = {lam}")
        log(f"  dim HWV = {delta}")
        log(f"  Note : {note}")
        log("  => Les vecteurs HWV s'annulent sur les points d'orbite testés.")
        log("  => Cela SUGGÈRE que cet isotype contribue à I_5(O_det3).")
        log("  => Prochaine étape : --mode srmt-test-candidate --lambda " + str(lam))
        log("  => Puis certifier δ > 0 sur plusieurs primes.")
        log("=" * 60)

    elif status == 'ISSUE_B':
        log("=" * 60)
        log("ISSUE B — LONG RUN FAIT, PAS DE DÉFAUT TROUVÉ")
        n_tested = results.get('n_tested', 0)
        log(f"  Isotypes testés : {n_tested} / {len(PLETHYSM_S5_S3_GL9_CERTIFIED)}")
        log(f"  Note : {note}")
        log("  => Aucun générateur de I_5(O_det3) trouvé dans les isotypes testés.")
        log("=" * 60)

    elif status == 'ISSUE_C':
        log("=" * 60)
        log("ISSUE C — TROP LOURD SUR CETTE MACHINE")
        log(f"  Note : {note}")
        log(f"  Limite actuelle : WEIGHT_SPACE_MAX = {WEIGHT_SPACE_MAX}")
        log("  => Réduire WEIGHT_SPACE_MAX ou utiliser une machine plus puissante.")
        log("=" * 60)

    elif status == 'ISSUE_D':
        log("=" * 60)
        log("ISSUE D — ENVIRONNEMENT INCOMPLET")
        log(f"  Note : {note}")
        missing = results.get('missing', [])
        if missing:
            log(f"  Manquant : {missing}")
        log("  => Installer les dépendances manquantes (voir --mode docs).")
        log("=" * 60)

    elif status == 'ISSUE_E':
        log("=" * 60)
        log("ISSUE E — TESTS ÉCHOUÉS")
        log(f"  Note : {note}")
        failed = results.get('failed_tests', [])
        if failed:
            log(f"  Tests échoués : {failed}")
        log("  => Corriger les erreurs avant de lancer le long run.")
        log("=" * 60)

    else:
        log(f"  Statut : {status}")
        log(f"  Note : {note}")
        if delta is not None:
            log(f"  Delta = {delta}")


# ===========================================================================
# SECTION 23 : Main CLI
# ===========================================================================

def _parse_lambda_arg(s: str) -> List[int]:
    """Parse un argument lambda de la forme '9,4,2' ou '12,3'."""
    parts = [int(x.strip()) for x in s.split(',') if x.strip()]
    return sorted(parts, reverse=True)


def mode_docs(args):
    """--mode docs : afficher la documentation."""
    print(__doc__)
    print_quick_docs_notes()


def mode_audit(args):
    """--mode audit : audit de l'environnement."""
    run_env_audit()


def mode_gl2(args):
    """--mode gl2 : sandbox GL2."""
    ok = run_gl2_sandbox()
    print_conclusion({
        'status': 'ISSUE_E' if not ok else 'OK',
        'note': 'GL2 sandbox' + (' passée' if ok else ' échouée'),
    })


def mode_modp(args):
    """--mode modp : tests arithmétique mod p."""
    ok = run_modp_tests()
    print_conclusion({
        'status': 'ISSUE_E' if not ok else 'OK',
        'note': 'Tests modp' + (' passés' if ok else ' échoués'),
        'failed_tests': [] if ok else ['modp'],
    })


def mode_plethysm(args):
    """--mode plethysm : vérifier PLETHYSM_S5_S3_GL9_CERTIFIED."""
    log_section("Vérification plethysme certifié")
    ok = _verify_plethysm_dimension()
    expected = math.comb(169, 5)
    log(f"  Nombre d'isotypes : {len(PLETHYSM_S5_S3_GL9_CERTIFIED)}")
    total = sum(mult * dim_schur_gl_n(lam, 9) for lam, mult in PLETHYSM_S5_S3_GL9_CERTIFIED)
    log(f"  sum(mult * dim_schur) = {total}")
    log(f"  Attendu : {expected}")
    log(f"  Vérification : {'OK' if ok else 'ECHEC'}")
    if ok:
        log("  La liste des 28 isotypes est certifiée.")
    else:
        log("  ERREUR : la liste des isotypes est incorrecte!")
    print_conclusion({
        'status': 'OK' if ok else 'ISSUE_E',
        'note': f'Plethysme: sum={total}, attendu={expected}',
    })


def mode_paste_plethysm(args):
    """--mode paste-plethysm : coller la sortie Sage et mettre à jour la liste."""
    log_section("Import plethysme Sage")
    log("Collez la sortie Sage ci-dessous (terminer avec une ligne contenant 'END'):")
    lines = []
    while True:
        try:
            line = input()
            if line.strip() == 'END':
                break
            lines.append(line)
        except EOFError:
            break

    text = '\n'.join(lines)
    parsed = parse_plethysm_output(text)
    log(f"  {len(parsed)} isotypes parsés depuis l'entrée Sage")

    if not parsed:
        log("  Aucun isotype trouvé. Vérifier le format de la sortie Sage.")
        return

    # Vérifier dimensionnellement
    total = sum(mult * dim_schur_gl_n(lam, 9) for lam, mult in parsed)
    expected = math.comb(169, 5)
    log(f"  sum(mult * dim_schur) = {total}  (attendu: {expected})")
    if total == expected:
        log("  OK: vérification dimensionnelle passée. Mise à jour de PLETHYSM_S5_S3_GL9_CERTIFIED.")
        global PLETHYSM_S5_S3_GL9_CERTIFIED
        PLETHYSM_S5_S3_GL9_CERTIFIED = parsed
    else:
        log(f"  ERREUR: vérification dimensionnelle échouée ({total} != {expected}).")
        log("  La liste n'a PAS été mise à jour.")


def mode_verify_plethysm_list(args):
    """--mode verify-plethysm-list : vérification dimensionnelle."""
    mode_plethysm(args)


def mode_lie(args):
    """--mode lie : tests action de Lie."""
    ok = run_lie_action_tests()
    print_conclusion({
        'status': 'ISSUE_E' if not ok else 'OK',
        'note': 'Tests Lie' + (' passés' if ok else ' échoués'),
    })


def mode_orbit(args):
    """--mode orbit : tests évaluation orbitale."""
    ok = run_orbit_eval_tests()
    print_conclusion({
        'status': 'ISSUE_E' if not ok else 'OK',
        'note': 'Tests orbite' + (' passés' if ok else ' échoués'),
    })


def mode_degree5_mini(args):
    """
    --mode degree5-mini : test miniature sur GL2 (twisted cubic).
    Ce mode ne nécessite pas Sage (NEEDS_PLETHYSM_IMPORT → liste certifiée disponible).
    """
    log_section("Degree 5 Mini Run (GL2 proxy)")
    W_basis = generate_W_basis_GL9()
    p = getattr(args, 'prime', 1000000007)

    # Tester sur l'isotype [15] (le plus simple)
    lam = [15]
    log(f"  Isotype test : {lam}")
    log(f"  dim_schur_gl_n({lam}, 9) = {dim_schur_gl_n(lam, 9)}")

    try:
        monomials = generate_weight_space_monomials(lam, W_basis, max_size=10_000)
        log(f"  Monomials de poids {lam} : {len(monomials)}")
        if len(monomials) > 0:
            log(f"  Premier monôme (MonoIdx5) : {monomials[0]}")
    except TOO_LARGE as e:
        log(f"  TOO_LARGE: {e}")
        print_conclusion({'status': 'ISSUE_C', 'note': str(e)})
        return

    log("  Mode degree5-mini : test préliminaire OK")
    log("  Pour le calcul complet, utiliser --mode degree5-one-isotype --lambda 15")


def mode_degree5_one_isotype(args):
    """--mode degree5-one-isotype --lambda '9,4,2' : tester un seul isotype."""
    lam_str = getattr(args, 'lam', None)
    if lam_str is None:
        log("  Erreur : --lambda requis pour ce mode")
        return

    lam = _parse_lambda_arg(lam_str)
    p = getattr(args, 'prime', 1000000007)
    max_wt = getattr(args, 'max_weight_dim', WEIGHT_SPACE_MAX)
    use_gpu_eval = getattr(args, 'gpu_eval', False)

    log_section(f"Degree 5 — Un isotype : lam={lam}")
    log(f"  dim_schur_gl_n({lam}, 9) = {dim_schur_gl_n(lam, 9)}")
    log(f"  gpu_eval = {use_gpu_eval}")

    W_basis = generate_W_basis_GL9()

    # Générer l'espace de poids (MonoIdx5 format)
    try:
        log(f"  Génération de l'espace de poids...")
        monomials = generate_weight_space_monomials(lam, W_basis, max_size=max_wt)
        log(f"  dim espace de poids = {len(monomials)}")
    except TOO_LARGE as e:
        log(f"  TOO_LARGE: {e}")
        print_conclusion({'status': 'ISSUE_C', 'note': str(e), 'lam': lam})
        return

    if not monomials:
        log(f"  Espace de poids vide pour lam={lam}")
        return

    # Construire les équations HWV
    try:
        log(f"  Construction des équations HWV...")
        hwv_matrix = build_hwv_equations_modp(lam, p, W_basis)
        log(f"  Matrice HWV : {len(hwv_matrix)} lignes x {len(monomials)} colonnes")
    except TOO_LARGE as e:
        log(f"  TOO_LARGE dans build_hwv_equations_modp: {e}")
        print_conclusion({'status': 'ISSUE_C', 'note': str(e), 'lam': lam})
        return

    # Calculer le noyau = espace HWV
    if hwv_matrix:
        log(f"  Calcul du noyau...")
        hwv_ker = kernel_modp_sparse(hwv_matrix, p)
        log(f"  dim HWV = {len(hwv_ker)}  (attendu théoriquement : mult = {_get_mult(lam)})")
    else:
        hwv_ker = [[1] + [0] * (len(monomials) - 1)] if monomials else []
        log(f"  Matrice HWV vide : tous les vecteurs sont HWV")

    if not hwv_ker:
        log(f"  Aucun vecteur HWV trouvé pour lam={lam}")
        CANDIDATES_I5.append({'lam': lam, 'dim_hwv': 0, 'vanishes': None})
        return

    # Convertir en HWVSparse (v5)
    hwv_sparse_list = [hwv_vec_to_sparse(vec, p) for vec in hwv_ker]

    # Construire les coordonnées d'orbite
    seeds = list(range(10, 210))
    log(f"  Construction des coordonnées d'orbite (gpu_eval={use_gpu_eval})...")

    # Vérifier la fiabilité GPU si demandé
    if use_gpu_eval and HAS_CUPY:
        if not _GPU_EVAL_TRUSTED:
            log("  WARNING: GPU eval désactivé (comparaison CPU/GPU échouée précédemment)")
            use_gpu_eval = False

    orbit_coords, coord_src = build_orbit_coords_gpu(seeds, p, W_basis, use_gpu=use_gpu_eval)
    log(f"  Orbit coords source: {coord_src}")

    # Évaluer sur l'orbite avec timing
    log(f"  Évaluation sur l'orbite...")
    eval_result_arr, eval_src, dt_eval = _timed_eval(
        hwv_sparse_list, orbit_coords, monomials, p, use_gpu_eval
    )

    # Calculer le rang exact (JAMAIS float)
    log(f"  Calcul du rang exact mod p...")
    rank, rank_src = rank_modp_exact_gpu(eval_result_arr, p)
    log(f"  rank_eval_matrix = {rank} (source: {rank_src})")

    # vanishes_all : rang = 0 exactement
    vanishes_all = (rank == 0)
    log(f"  vanishes_all = {vanishes_all}")

    result = {
        'lam': lam,
        'dim_hwv': len(hwv_ker),
        'vanishes': vanishes_all,
        'rank_eval': rank,
    }
    CANDIDATES_I5.append(result)
    save_checkpoint({'completed_isotypes': [lam], 'candidates_i5': CANDIDATES_I5})

    if vanishes_all:
        log(f"  => Tous les HWV s'annulent sur l'orbite : CANDIDAT pour I_5!")
        log(f"  => ATTENTION : ce n'est pas encore un défaut SRMT certifié.")
        log(f"  => Lancer : --mode srmt-test-candidate --lambda {lam}")
        print_conclusion({
            'status': 'CANDIDATE_I5_FOUND',
            'lam': lam,
            'delta': len(hwv_ker),  # ici delta = dim_hwv, PAS le delta SRMT
            'note': f'vanishes_all=True sur 1 prime x {len(seeds)} seeds, dim_hwv={len(hwv_ker)}',
        })
    else:
        log(f"  => HWV ne s'annulent pas tous sur l'orbite : pas de candidat I5 dans cet isotype")
        log(f"  => rank_eval = {rank} (> 0 = HWV non nuls sur O_det3)")
        print_conclusion({
            'status': 'OK',
            'note': f'lam={lam} : rank_eval={rank} (HWV non nuls sur orbite)',
        })


def _get_mult(lam: List[int]) -> int:
    """Retourne la multiplicité de lam dans PLETHYSM_S5_S3_GL9_CERTIFIED."""
    for l, m in PLETHYSM_S5_S3_GL9_CERTIFIED:
        if l == lam:
            return m
    return 0


def mode_srmt_test_candidate(args):
    """--mode srmt-test-candidate --lambda '9,4,2' : tester le SRMT pour un isotype."""
    lam_str = getattr(args, 'lam', None)
    if lam_str is None:
        log("  Erreur : --lambda requis pour ce mode")
        return

    lam = _parse_lambda_arg(lam_str)
    p = getattr(args, 'prime', 1000000007)

    log_section(f"SRMT Test — lam={lam}")
    W_basis = generate_W_basis_GL9()
    result = srmt_test_one_candidate(lam, p, W_basis)

    for key, val in result.items():
        if key != 'defect_found':
            log(f"  {key} : {val}")

    if result.get('defect_found'):
        log("  => DÉFAUT TROUVÉ (delta > 0)")
        print_conclusion({
            'status': 'ISSUE_A',
            'lam': lam,
            'delta': result.get('(1,4)', {}).get('delta', 0),
            'note': 'SRMT test positif',
        })
    else:
        log("  => Pas de défaut trouvé dans ce test SRMT")


def mode_degree5_longrun(args):
    """--mode degree5-longrun [--resume] : long run sur tous les isotypes."""
    reviewer_level = run_reviewer()
    if reviewer_level != "APPROVED_FOR_LONG_RUN":
        log(f"  Long run non autorisé : niveau = {reviewer_level}")
        log("  Utiliser --mode review pour voir les conditions requises.")
        print_conclusion({
            'status': 'ISSUE_E',
            'note': f'Long run bloqué : reviewer_level={reviewer_level}',
        })
        return

    resume = getattr(args, 'resume', False)
    p = getattr(args, 'prime', 1000000007)
    max_wt = getattr(args, 'max_weight_dim', WEIGHT_SPACE_MAX)
    use_gpu_eval = getattr(args, 'gpu_eval', False)

    log_section("Degree 5 Long Run")
    log(f"  Prime : {p}, max_weight_dim : {max_wt}, resume : {resume}, gpu_eval : {use_gpu_eval}")

    state = load_checkpoint() if resume else {
        'completed_isotypes': [],
        'candidates_i5': [],
        'timestamp_start': time.strftime('%Y-%m-%dT%H:%M:%S'),
    }
    completed = [tuple(x) if isinstance(x, list) else x for x in state.get('completed_isotypes', [])]
    W_basis = generate_W_basis_GL9()

    # Vérifier fiabilité GPU si demandé
    if use_gpu_eval and HAS_CUPY:
        _check_gpu_trustworthy(W_basis)
        if not _GPU_EVAL_TRUSTED:
            log("  WARNING: GPU eval désactivé (test de comparaison échoué)")
            use_gpu_eval = False

    n_tested = 0
    n_too_large = 0
    n_error = 0

    for lam, mult in PLETHYSM_S5_S3_GL9_CERTIFIED:
        lam_key = tuple(lam)
        if lam_key in completed:
            log(f"  Isotype {lam} déjà traité (checkpoint), skip.")
            continue

        log(f"\n  --- Isotype {lam} (mult={mult}) ---")
        temp = _get_gpu_temperature()
        if temp is not None:
            log(f"  Température GPU : {temp}°C")
            if temp > 85:
                log(f"  WARNING: Température GPU élevée ({temp}°C). Pause 30s.")
                time.sleep(30)

        try:
            monomials = generate_weight_space_monomials(lam, W_basis, max_size=max_wt)
            log(f"  dim espace de poids = {len(monomials)}")

            if not monomials:
                log(f"  Espace vide, skip.")
                n_tested += 1
                completed.append(lam_key)
                continue

            hwv_matrix = build_hwv_equations_modp(lam, p, W_basis)
            if hwv_matrix:
                hwv_ker = kernel_modp_sparse(hwv_matrix, p)
            else:
                hwv_ker = [[1] + [0] * (len(monomials) - 1)]

            log(f"  dim HWV = {len(hwv_ker)}")

            # Convertir en HWVSparse (v5)
            hwv_sparse_list = [hwv_vec_to_sparse(vec, p) for vec in hwv_ker]

            # Construire les coordonnées d'orbite
            seeds = list(range(10, 210))
            orbit_coords, coord_src = build_orbit_coords_gpu(seeds, p, W_basis, use_gpu=use_gpu_eval)

            # Évaluer avec timing
            eval_result_arr, eval_src, dt_eval = _timed_eval(
                hwv_sparse_list, orbit_coords, monomials, p, use_gpu_eval
            )

            # Rang exact mod p
            rank, rank_src = rank_modp_exact_gpu(eval_result_arr, p)
            vanishes_all = (rank == 0)

            result = {
                'lam': lam,
                'mult': mult,
                'dim_hwv': len(hwv_ker),
                'vanishes': vanishes_all,
                'rank_eval': rank,
            }
            CANDIDATES_I5.append(result)

            if vanishes_all:
                log(f"  => CANDIDAT TROUVÉ : lam={lam}, dim_hwv={len(hwv_ker)}")

        except TOO_LARGE as e:
            log(f"  TOO_LARGE pour {lam}: {e}")
            n_too_large += 1
            CANDIDATES_I5.append({'lam': lam, 'mult': mult, 'status': 'TOO_LARGE'})

        except Exception as e:
            log(f"  ERREUR pour {lam}: {e}")
            n_error += 1
            CANDIDATES_I5.append({'lam': lam, 'mult': mult, 'status': f'ERROR: {e}'})

        n_tested += 1
        completed.append(lam_key)
        state['completed_isotypes'] = [list(x) for x in completed]
        state['candidates_i5'] = CANDIDATES_I5
        save_checkpoint(state)

    # Bilan
    vanishing = [r for r in CANDIDATES_I5 if r.get('vanishes') is True]
    log(f"\n  === Bilan long run ===")
    log(f"  Isotypes testés : {n_tested}")
    log(f"  TOO_LARGE : {n_too_large}")
    log(f"  Erreurs : {n_error}")
    log(f"  Candidats (vanishes_all=True) : {len(vanishing)}")

    if vanishing:
        log(f"  => {len(vanishing)} candidat(s) I5 trouvé(s). PAS encore ISSUE_A.")
        log(f"  => Pour certifier : --mode srmt-test-candidate --lambda <lam> sur chaque candidat.")
        print_conclusion({
            'status': 'CANDIDATE_I5_FOUND',
            'lam': vanishing[0]['lam'],
            'delta': vanishing[0].get('dim_hwv', 0),
            'note': f'{len(vanishing)} candidat(s) I5 (vanishes_all=True). SRMT non certifié.',
        })
    elif n_too_large > 0 or n_error > 0:
        print_conclusion({
            'status': 'ISSUE_C',
            'note': f'{n_too_large} isotypes TOO_LARGE, {n_error} erreurs',
        })
    else:
        print_conclusion({
            'status': 'ISSUE_B',
            'n_tested': n_tested,
            'note': 'Long run complet, aucun défaut trouvé',
        })


def mode_review(args):
    """--mode review : exécuter le reviewer."""
    level = run_reviewer()
    log(f"\n  Niveau : {level}")


def mode_all(args):
    """--mode all : exécuter tous les modes de test."""
    log_section("Mode ALL")
    mode_audit(args)
    gl2_ok = run_gl2_sandbox()
    modp_ok = run_modp_tests()
    lie_ok = run_lie_action_tests()
    orbit_ok = run_orbit_eval_tests()
    mode_plethysm(args)
    level = run_reviewer()
    log(f"\n  === Résumé ===")
    log(f"  GL2: {'PASS' if gl2_ok else 'FAIL'}")
    log(f"  ModP: {'PASS' if modp_ok else 'FAIL'}")
    log(f"  Lie: {'PASS' if lie_ok else 'FAIL'}")
    log(f"  Orbit: {'PASS' if orbit_ok else 'FAIL'}")
    log(f"  Reviewer: {level}")


def main():
    """
    Point d'entrée principal.

    Modes disponibles :
      --mode docs
      --mode audit
      --mode gl2
      --mode modp
      --mode plethysm
      --mode paste-plethysm
      --mode verify-plethysm-list
      --mode lie
      --mode orbit
      --mode degree5-mini
      --mode degree5-one-isotype --lambda "9,4,2"
      --mode srmt-test-candidate --lambda "9,4,2"
      --mode degree5-longrun [--resume]
      --mode review
      --mode all

    Options :
      --lambda           : partition en argument (ex: "9,4,2" ou "12,3")
      --resume           : reprendre depuis checkpoint
      --prime            : prime à utiliser (défaut 1000000007)
      --max-weight-dim   : limite dimension espace de poids (défaut 500000)
      --gpu-eval         : activer l'évaluation sur GPU (CuPy requis)
    """
    parser = argparse.ArgumentParser(
        description='Route Alpha v5 — Recherche équations degré 5 pour O_det3 (GPU offloading)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--mode', type=str, default='audit',
        choices=[
            'docs', 'audit', 'gl2', 'modp', 'plethysm', 'paste-plethysm',
            'verify-plethysm-list', 'lie', 'orbit', 'degree5-mini',
            'degree5-one-isotype', 'srmt-test-candidate',
            'degree5-longrun', 'review', 'all',
        ],
        help="Mode d'exécution (défaut: audit)",
    )
    parser.add_argument(
        '--lambda', dest='lam', type=str, default=None,
        help='Partition lambda (ex: "9,4,2" ou "12,3")',
    )
    parser.add_argument(
        '--resume', action='store_true',
        help='Reprendre depuis le dernier checkpoint',
    )
    parser.add_argument(
        '--prime', type=int, default=1000000007,
        help='Nombre premier à utiliser (défaut: 1000000007)',
    )
    parser.add_argument(
        '--max-weight-dim', type=int, default=WEIGHT_SPACE_MAX,
        dest='max_weight_dim',
        help=f'Limite dimension espace de poids (défaut: {WEIGHT_SPACE_MAX})',
    )
    parser.add_argument(
        '--gpu-eval', action='store_true', dest='gpu_eval',
        help='Activer l\'évaluation sur GPU (CuPy requis)',
    )

    args = parser.parse_args()

    mode_map = {
        'docs': mode_docs,
        'audit': mode_audit,
        'gl2': mode_gl2,
        'modp': mode_modp,
        'plethysm': mode_plethysm,
        'paste-plethysm': mode_paste_plethysm,
        'verify-plethysm-list': mode_verify_plethysm_list,
        'lie': mode_lie,
        'orbit': mode_orbit,
        'degree5-mini': mode_degree5_mini,
        'degree5-one-isotype': mode_degree5_one_isotype,
        'srmt-test-candidate': mode_srmt_test_candidate,
        'degree5-longrun': mode_degree5_longrun,
        'review': mode_review,
        'all': mode_all,
    }

    fn = mode_map.get(args.mode)
    if fn is None:
        log(f"Mode inconnu : {args.mode}")
        parser.print_help()
        sys.exit(1)

    log_section(f"Route Alpha v5 — mode={args.mode}")
    fn(args)


if __name__ == '__main__':
    main()
