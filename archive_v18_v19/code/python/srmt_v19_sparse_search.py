#!/usr/bin/env python3
# =====================================================================
# srmt_v19_sparse_search.py
#
# Recherche sparse de candidats quartiques dans I_4(\overline{GL_9.det_3}),
# version Python autonome. Macaulay2 n'est utilise qu'a la fin pour
# certifier 1-5 candidats survivants via testJet(P, 4).
#
# Pipeline :
#   1. Lire mons3 (165 monomes degre 3 en y_(0,0)..y_(2,2)) et les 81 matrices
#      gl_9 sparse depuis code/m2_data.txt (dumpe par M2).
#   2. Choisir un carrier W = ensemble de quartiques c_i c_j c_k c_l
#      (i <= j <= k <= l), structure par distance Lie depuis support det3.
#   3. Decomposer W par poids torique (le tore de SL_3 x SL_3 agit
#      diagonalement ; on raffine en utilisant tout le tore diagonal de GL_9
#      = ZZ^9, ce qui donne une decomposition tres fine).
#   4. Pour chaque bloc de poids w :
#        a. Calculer la matrice de contraintes jet<=1 modulo p=32003.
#        b. Si dim ker > 0, raffiner avec jet<=2.
#        c. Si dim ker <= 5 et > 0 a l'ordre 2, exporter chaque candidat
#           individuel reconstruit sur QQ vers un fichier .m2.
#   5. Cote M2 : pour chaque .m2 candidat, lancer testJet(P, 4).
#
# Conventions :
#   * c_alpha = variable polynomiale correspondant au monome mons3#alpha de R = QQ[y_00..y_22].
#   * Quartique = combinaison lineaire de monomes c_i c_j c_k c_l (i<=j<=k<=l).
#     Stocke comme dict {tuple_trie: int} (coeffs entiers).
#   * Derivation D_ab sur S = QQ[c_0..c_164] : c_k -> sum_i M_{ab}[i,k] c_i.
#   * Pour un quartique m = c_i c_j c_k c_l, D_ab(m) = somme sur les 4 facteurs
#     remplaces, par regle de Leibniz.
#   * Evaluation en vDet3 : un monome contribue ssi ses 4 indices sont dans
#     SUPPORT_DET3 = {34, 37, 65, 72, 92, 96}.
#   * Pruning : si rDeriv := derivees restantes, et un monome a > rDeriv indices
#     hors SUPPORT_DET3, il finira a zero, on le supprime.
#
# Statuts (alignes editorial v18) :
#   PROVED-ALG   : verifie algebriquement (rang tangent, ...).
#   QQ-CERTIFIED : un candidat P passe testJet(P, 4) sur QQ exact.
#   Computational obs. : obtenu modulo p ou avec kmax < 4.
# =====================================================================

import sys
import os
from collections import defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement, product

DATA_PATH = "./m2_data.txt"
CANDIDATES_DIR = "./candidates"
PRIME = 32003   # corps fini pour filtrage rapide
N3 = 165        # dim S_1 = dim Sym^3 C^9 = 165

# =====================================================================
# 1. Chargement des donnees M2
# =====================================================================

def load_m2_data(path=DATA_PATH):
    """
    Retourne :
      mons3_exp : list[tuple[int]*9]  # exposants des 165 monomes y^alpha
      gl_mats   : dict[(a,b)] -> dict[k] -> dict[i] -> int
                  (action de E_{ab} sur S_1 : c_k -> sum_i M[i,k] c_i)
      vDet3     : list[int] de longueur 165 (coefficients det3 dans mons3)
      support_det3 : set des indices alpha avec vDet3[alpha] != 0
    """
    with open(path, "r") as f:
        lines = f.readlines()

    idx = 0
    n3 = int(lines[idx].strip())
    idx += 1
    assert n3 == N3, f"unexpected N3={n3}"

    # 165 lignes d'exposants
    mons3_exp = []
    for i in range(N3):
        e = tuple(int(x) for x in lines[idx].split())
        assert len(e) == 9 and sum(e) == 3
        mons3_exp.append(e)
        idx += 1

    gl_mats = {}
    while True:
        line = lines[idx].strip()
        if line == "ENDGL":
            idx += 1
            break
        a, b = map(int, line.split())
        idx += 1
        # entrees jusqu'a "ENDM"
        # M[(a,b)] : k -> {i: coef}  (action sur c_k = colonne k)
        col = defaultdict(dict)
        while True:
            l = lines[idx].strip()
            idx += 1
            if l == "ENDM":
                break
            i, k, v = l.split()
            i, k, v = int(i), int(k), int(v)
            col[k][i] = v
        gl_mats[(a, b)] = dict(col)

    vDet3 = [int(x) for x in lines[idx].split()]
    assert len(vDet3) == N3
    support_det3 = frozenset(i for i, c in enumerate(vDet3) if c != 0)
    return mons3_exp, gl_mats, vDet3, support_det3


# =====================================================================
# 2. Quartiques sparse + derivation
# =====================================================================
# Cle = tuple (i,j,k,l) trie, valeur = entier modulo p ou Fraction.
# On normalise toujours par sortie + suppression des 0.

def normalize(d):
    """Supprime les zeros et trie les cles."""
    return {k: v for k, v in d.items() if v != 0}


def add_inplace(d, k, v, mod=None):
    """d[k] += v, avec gestion du modulo et suppression des 0."""
    if v == 0:
        return
    nv = d.get(k, 0) + v
    if mod is not None:
        nv = nv % mod
    if nv == 0:
        if k in d:
            del d[k]
    else:
        d[k] = nv


def derive_quartic(quart, ab, gl_mats, mod=None):
    """
    Derivation D_ab d'un quartique sparse.
    D_ab(c_i c_j c_k c_l) = somme sur le facteur remplace :
      somme_x M[x,i] c_x c_j c_k c_l + ... (4 termes, un par facteur).
    Comme le tuple est trie, on retrie systematiquement.
    Retourne un nouveau dict, normalise.
    """
    M = gl_mats[ab]   # M[k] = {i: coef}, action sur c_k.
    out = {}
    for tup, coef in quart.items():
        # 4 positions (avec doublons gérés naturellement : on remplace
        # chaque occurrence indexee, ce qui donne le bon facteur multiplicatif
        # via le nombre d'occurrences).
        for pos in range(4):
            k_old = tup[pos]
            if k_old not in M:
                continue
            for i_new, m_ik in M[k_old].items():
                # Reconstruire le tuple en remplaçant la pos
                new_list = list(tup)
                new_list[pos] = i_new
                new_list.sort()
                new_tup = tuple(new_list)
                contrib = coef * m_ik
                add_inplace(out, new_tup, contrib, mod=mod)
    return out


# =====================================================================
# 3. Evaluation en vDet3 (support_det3, signes vDet3[i])
# =====================================================================

def eval_at_vDet3(quart, vDet3, support_det3, mod=None):
    """
    P(v_det3) = sum_{tup} coef * prod_{i in tup} vDet3[i].
    On somme uniquement sur les tuples entierement dans support_det3.
    """
    s = 0
    for tup, coef in quart.items():
        ok = True
        prod = 1
        for x in tup:
            if x not in support_det3:
                ok = False
                break
            prod *= vDet3[x]
        if ok:
            s += coef * prod
            if mod is not None:
                s %= mod
    if mod is not None:
        return s % mod
    return s


# =====================================================================
# 4. Pruning par budget de derivees restantes
# =====================================================================

def prune_quartic(quart, r_remaining, support_det3):
    """
    Apres r_remaining derivees encore a appliquer, un monome dont strictement
    plus de r_remaining indices sont hors support_det3 NE PEUT PAS etre
    ramene dans le support en r_remaining derivees (chaque derivee remplace
    UN indice, donc reduit le nombre d'indices hors support d'au plus 1
    a chaque pas).
    => on supprime ces monomes.

    NB : c'est un pruning sur (jet ordre k, P) -> evaluation finale en vDet3.
    Si r_remaining = 0, on garde uniquement les monomes a 4 indices in support.
    """
    out = {}
    for tup, coef in quart.items():
        n_out = sum(1 for x in tup if x not in support_det3)
        if n_out <= r_remaining:
            out[tup] = coef
    return out


# (apply_seq_with_pruning supprime : factorise dans constraint_value)


# =====================================================================
# 5. Decomposition par poids torique
# =====================================================================
# Le tore diagonal de GL_9 = ZZ^9 agit sur c_alpha avec le poids
# = exposant du monome y^alpha dans Sym^3 C^9.
# Un quartique c_i c_j c_k c_l a pour poids = sum des 4 exposants (vecteur ZZ^9).

def cubic_weight(alpha, mons3_exp):
    """Poids torique du monome cubique c_alpha = exposants de mons3_exp[alpha]."""
    return mons3_exp[alpha]


def quartic_weight(tup, mons3_exp):
    """Poids torique de c_i c_j c_k c_l (i,j,k,l = tup)."""
    w = [0] * 9
    for x in tup:
        e = mons3_exp[x]
        for j in range(9):
            w[j] += e[j]
    return tuple(w)


# Note : le poids d'un element de I_4(\overline{GL_9.det_3}) DOIT etre invariant
# sous le sous-tore qui stabilise det3, soit le tore diagonal restreint a la
# sous-algebre {diag(t,t,t,...)*det = t^3 det}, c'est-a-dire le sous-espace
# w = (3,3,3,3,3,3,3,3,3) - normalise. Plus precisement : det3 a poids
# (3,3,3) sous chaque facteur SL_3 du parabolique. Mais en travaillant avec
# tout GL_9, l'orbite GL_9.det3 contient des elements de TOUS les poids
# correspondant a des monomes y^alpha de degre 3.
#
# DONC : un polynome P de degre 4 dans I_4(\overline{GL_9.det_3}) n'a pas
# besoin d'etre torique-pur ! C'est seulement quand on cherche P qui annule
# det3 lui-meme (jet 0) que le poids de P doit egaler 4*(3,3,3,...,3)/9 ...
# Bref, restriction par poids = restriction au sous-espace torique-pur.
# On gagne en taille mais on perd des candidats.
#
# COMPROMIS PRATIQUE : decomposer le carrier W par poids et chercher le
# noyau bloc par bloc -- c'est CORRECT ssi les contraintes (D^seq P)(v_det3)
# sont torique-equivariantes, ce qui est le cas pour seq de poids fixe.
# Le poids de v_det3 est sum de 6 monomes de poids divers (chaque monome
# y^sigma a poids = perm(1,1,1,...) selon sigma).
# En fait, les 6 monomes y^sigma de det3 ont TOUS le meme poids
# (1,1,1,1,1,1,1,1,1)/3 ... Non : chaque y_{i,sigma(i)} a poids unitaire
# sur la coordonnee (i,sigma(i)) ; le produit donne un poids
# (1 sur 3 lignes, 1 sur 3 cols, 0 ailleurs) = somme = (1,1,1,1,1,1,1,1,1)
# ?? non, sigma permute 3 colonnes, donc le poids total est (1,1,1,1,1,1,1,1,1)
# vu sur 9 coords. Verifions par exemple pour sigma = id :
# y_{0,0} y_{1,1} y_{2,2} -> exposant 1 sur (0,0), (1,1), (2,2), 0 ailleurs.
# Donc poids = (1,0,0, 0,1,0, 0,0,1) en serialisant (i,j) -> 3i+j.
# Different de sigma=(01) : y_{0,1} y_{1,0} y_{2,2} -> (0,1,0, 1,0,0, 0,0,1).
# DIFFERENT POIDS. Donc les 6 monomes de det3 ont 6 poids DISTINCTS.
# v_det3 = somme des 6 poids -> N'EST PAS torique-pur.
#
# CONSEQUENCE : la decomposition par poids n'est PAS triviale. Si P a un
# poids w_P, alors (D^seq P)(v_det3) = somme sur les 6 sigma de
# coefficient[poids w_P + somme(seq) - poids(sigma) = 0].
# Donc cette evaluation est non-nulle UNIQUEMENT si w_P + somme(seq)
# matche un des 6 poids de det3.
# CELA DONNE UNE CONTRAINTE FORTE et reduit l'espace.

# Pour l'instant, on raffine W non par poids torique strict mais par
# structure combinatoire (W0, W1, W2 par profil de support).


# =====================================================================
# 6. Construction du carrier
# =====================================================================

def support_levels(mons3_exp, gl_mats, support_det3, max_level=2):
    """
    Calcule IDX0 = support_det3, IDX1 = indices atteignables par UN
    coup gl_9 depuis IDX0, IDX2 = atteignables en 2 coups, etc.
    Atteignable = i tel qu'il existe (a,b) et k in IDX_prev avec M_{ab}[i,k] != 0.
    Retourne dict level -> frozenset.
    """
    levels = {0: frozenset(support_det3)}
    seen = set(support_det3)
    for L in range(1, max_level + 1):
        new = set()
        for ab, M in gl_mats.items():
            for k in levels[L - 1]:
                if k in M:
                    for i in M[k].keys():
                        if i not in seen:
                            new.add(i)
        seen |= new
        levels[L] = frozenset(new)
    return levels


def build_carrier(idx_levels, profile):
    """
    profile : tuple (n0, n1, n2) avec n0+n1+n2 = 4.
    Retourne la liste des tuples (i,j,k,l) tries, avec
       n0 indices dans IDX0, n1 dans IDX1, n2 dans IDX2.
    """
    n0, n1, n2 = profile
    assert n0 + n1 + n2 == 4
    L0 = sorted(idx_levels[0])
    L1 = sorted(idx_levels[1]) if 1 in idx_levels else []
    L2 = sorted(idx_levels[2]) if 2 in idx_levels else []
    out = []
    for c0 in combinations_with_replacement(L0, n0):
        for c1 in combinations_with_replacement(L1, n1):
            for c2 in combinations_with_replacement(L2, n2):
                tup = tuple(sorted(list(c0) + list(c1) + list(c2)))
                out.append(tup)
    # dedup (puisqu'un meme tuple peut venir de partitions differentes
    # quand n_i a doublons d'un IDX a l'autre -- ici disjoints, donc OK)
    return sorted(set(out))


# =====================================================================
# 7. Sequences Lie
# =====================================================================

def all_lie_sequences(kmax):
    """
    Toutes les sequences (E_{a1,b1}, ..., E_{ak,bk}) avec k <= kmax.
    Format : list de tuples-de-tuples.
    Compte = sum_{k=0}^{kmax} 81^k.
    """
    pairs = [(a, b) for a in range(9) for b in range(9)]
    seqs = [()]
    cur = [()]
    for _ in range(kmax):
        cur = [s + (p,) for s in cur for p in pairs]
        seqs += cur
    return seqs


# =====================================================================
# 8. Construction de la matrice de contraintes (sparse, mod p)
# =====================================================================
# Pour W = liste de tuples quartiques c1, c2, ..., cn :
# contraintes lignes indexees par seq Lie de longueur <= kmax,
# colonnes indexees par les n carriers.
# C[seq, alpha] = ((D^seq c_alpha)(v_det3)) mod p.
# noyau lineaire de C = combinaisons des c_alpha qui sont dans I^{<=kmax}.

def constraint_value(carrier_tup, seq, gl_mats, vDet3, support_det3, mod):
    """
    Pour le QUARTIQUE c_{i1} c_{i2} c_{i3} c_{i4} (carrier_tup tri\u00e9),
    calcule (D^seq P)(v_det3) mod p, avec pruning a chaque etape.
    """
    # Q = dict{carrier_tup: 1}
    Q = {carrier_tup: 1}
    L = len(seq)
    for step, ab in enumerate(seq):
        if not Q:
            return 0
        Q = derive_quartic(Q, ab, gl_mats, mod=mod)
        r_rem = L - step - 1
        Q = prune_quartic(Q, r_rem, support_det3)
    if not Q:
        return 0
    return eval_at_vDet3(Q, vDet3, support_det3, mod=mod)


def build_constraints_matrix(W, seqs, gl_mats, vDet3, support_det3, mod):
    """
    Retourne une matrice numpy ou liste-de-listes mod p :
      C[s, j] = constraint_value(W[j], seqs[s], ...) mod p.
    Renvoie comme list de lignes (chaque ligne = list de longueur n).
    """
    n = len(W)
    rows = []
    for s_idx, seq in enumerate(seqs):
        row = [constraint_value(W[j], seq, gl_mats, vDet3, support_det3, mod)
               for j in range(n)]
        rows.append(row)
    return rows


# =====================================================================
# 9. Noyau modulo p (Gauss-Jordan integrer)
# =====================================================================

def modinv(a, p):
    return pow(a, -1, p)


def mat_kernel_mod_p(rows, mod):
    """
    rows : list de m lignes, chacune longueur n, entries dans Z/p.
    Retourne (rang, list de vecteurs base du noyau, longueur n chacun).
    """
    if not rows:
        return 0, []
    m = len(rows)
    n = len(rows[0])
    # Reduction sur des copies
    A = [list(r) for r in rows]
    pivots = [-1] * m   # ligne i pivote sur colonne pivots[i]
    col_pivot = {}      # col -> ligne qui la pivote
    for col in range(n):
        # Cherche un pivot dans cette colonne, parmi lignes non encore pivotees
        pivot_row = None
        for i in range(m):
            if pivots[i] == -1 and A[i][col] % mod != 0:
                pivot_row = i
                break
        if pivot_row is None:
            continue
        # Normalise
        inv = modinv(A[pivot_row][col] % mod, mod)
        A[pivot_row] = [(x * inv) % mod for x in A[pivot_row]]
        # Eliminer dans les autres lignes
        for i in range(m):
            if i == pivot_row:
                continue
            v = A[i][col] % mod
            if v != 0:
                A[i] = [(A[i][k] - v * A[pivot_row][k]) % mod for k in range(n)]
        pivots[pivot_row] = col
        col_pivot[col] = pivot_row
    rank = sum(1 for p_ in pivots if p_ >= 0)
    # Variables libres = colonnes non pivot
    free_cols = [c for c in range(n) if c not in col_pivot]
    kernel = []
    for fc in free_cols:
        v = [0] * n
        v[fc] = 1
        # Pour chaque colonne pivot, valeur = -A[pivot_row][fc]
        for pc, pr in col_pivot.items():
            v[pc] = (-A[pr][fc]) % mod
        kernel.append(v)
    return rank, kernel


# =====================================================================
# 10. Reconstruction QQ d'un candidat (apres filtrage mod p)
# =====================================================================
# Si on a un noyau de petite dimension, on relance le calcul de la matrice
# de contraintes sur QQ (Fraction), et on prend le noyau exact.

def constraint_value_qq(carrier_tup, seq, gl_mats, vDet3, support_det3):
    """Idem constraint_value, mais sans mod (ZZ)."""
    Q = {carrier_tup: 1}
    L = len(seq)
    for step, ab in enumerate(seq):
        if not Q:
            return 0
        Q = derive_quartic(Q, ab, gl_mats, mod=None)
        r_rem = L - step - 1
        Q = prune_quartic(Q, r_rem, support_det3)
    if not Q:
        return 0
    return eval_at_vDet3(Q, vDet3, support_det3, mod=None)


def kernel_qq(W, seqs, gl_mats, vDet3, support_det3):
    """
    Construit la matrice de contraintes sur QQ (entrees ZZ -> resoudre via
    Gauss). Retourne (rank, kernel_basis) avec kernel_basis = list de
    vecteurs de longueur n (Fraction).
    """
    n = len(W)
    rows = []
    for seq in seqs:
        row = [Fraction(constraint_value_qq(W[j], seq, gl_mats, vDet3, support_det3))
               for j in range(n)]
        rows.append(row)
    if not rows:
        return 0, []
    m = len(rows)
    A = [list(r) for r in rows]
    pivots = [-1] * m
    col_pivot = {}
    for col in range(n):
        pivot_row = None
        for i in range(m):
            if pivots[i] == -1 and A[i][col] != 0:
                pivot_row = i
                break
        if pivot_row is None:
            continue
        inv = Fraction(1, 1) / A[pivot_row][col]
        A[pivot_row] = [x * inv for x in A[pivot_row]]
        for i in range(m):
            if i == pivot_row:
                continue
            v = A[i][col]
            if v != 0:
                A[i] = [A[i][k] - v * A[pivot_row][k] for k in range(n)]
        pivots[pivot_row] = col
        col_pivot[col] = pivot_row
    rank = sum(1 for p_ in pivots if p_ >= 0)
    free_cols = [c for c in range(n) if c not in col_pivot]
    kernel = []
    for fc in free_cols:
        v = [Fraction(0)] * n
        v[fc] = Fraction(1)
        for pc, pr in col_pivot.items():
            v[pc] = -A[pr][fc]
        kernel.append(v)
    return rank, kernel


def reconstruct_polynomial_from_kernel_vec(kvec, W):
    """
    kvec : list de Fraction de longueur len(W).
    Retourne dict{tuple_trie: Fraction} = polynome P = sum_j kvec[j] * c_W[j].
    """
    P = {}
    for j, c in enumerate(kvec):
        if c != 0:
            tup = W[j]
            P[tup] = P.get(tup, Fraction(0)) + c
    return {k: v for k, v in P.items() if v != 0}


def normalize_to_integers(P):
    """Echelonne par PPCM des denominateurs et divise par PGCD des numerateurs."""
    from math import gcd
    if not P:
        return P
    # PPCM des denoms
    denom = 1
    for v in P.values():
        denom = denom * v.denominator // gcd(denom, v.denominator)
    Q = {k: int(v.numerator * (denom // v.denominator)) for k, v in P.items()}
    # PGCD des nums
    g = 0
    for v in Q.values():
        g = gcd(g, abs(v))
    if g > 1:
        Q = {k: v // g for k, v in Q.items()}
    # Normaliser le signe (premier coef > 0)
    keys = sorted(Q.keys())
    if Q[keys[0]] < 0:
        Q = {k: -v for k, v in Q.items()}
    return Q


# =====================================================================
# 11. Export d'un candidat vers .m2
# =====================================================================

def export_candidate_m2(P_dict, idx, profile, weight, kmax_filter, out_dir=CANDIDATES_DIR):
    """
    Ecrit un fichier candidate_<idx>.m2 contenant :
      - la definition de l'anneau S (165 c_i)
      - l'expression Pcandidate = ...
      - l'invocation de testJet (suppose load du harness avant).
    """
    os.makedirs(out_dir, exist_ok=True)
    fpath = os.path.join(out_dir, f"candidate_{idx:03d}.m2")
    with open(fpath, "w") as f:
        f.write(f"-- candidate_{idx:03d}.m2\n")
        f.write(f"-- profile = {profile}, weight = {weight}, "
                f"filtered at kmax = {kmax_filter} (mod p + QQ ker)\n")
        f.write(f"-- {len(P_dict)} monomials\n")
        f.write("Pcandidate = (\n")
        terms = []
        for tup, coef in sorted(P_dict.items()):
            cs = " * ".join(f"cVars#{x}" for x in tup)
            sign = "+" if coef > 0 else "-"
            terms.append(f"     {sign} {abs(coef)} * {cs}")
        # Premier terme : retirer le "+ "
        if terms[0].lstrip().startswith("+"):
            terms[0] = terms[0].replace("+ ", "  ", 1)
        f.write("\n".join(terms))
        f.write("\n     );\n")
        f.write('print("--- Candidate " | toString ' + str(idx) + ' | " ---");\n')
        f.write("print(toString Pcandidate);\n")
        f.write("jetResult := testJet(Pcandidate, 4);\n")
        f.write('print("kmax=4 result : " | toString jetResult);\n')
        f.write('<< "kmax=4 result : " << jetResult << endl << flush;\n')
    return fpath


# =====================================================================
# 12. Pipeline principal
# =====================================================================

def main():
    print("=== SRMT v19 sparse search ===", flush=True)
    print(f"  loading M2 data from {DATA_PATH}...", flush=True)
    mons3_exp, gl_mats, vDet3, support_det3 = load_m2_data()
    print(f"  N3 = {N3}, |support_det3| = {len(support_det3)}", flush=True)
    print(f"  support_det3 = {sorted(support_det3)}", flush=True)
    print(f"  vDet3 on support : "
          f"{[(i, vDet3[i]) for i in sorted(support_det3)]}", flush=True)

    # Levels
    print("  computing IDX levels (Lie distance from support)...", flush=True)
    idx_levels = support_levels(mons3_exp, gl_mats, support_det3, max_level=2)
    for L in sorted(idx_levels):
        print(f"    IDX{L} : {len(idx_levels[L])} indices", flush=True)

    # ATTENTION : IDX0, IDX1, IDX2 peuvent ne PAS etre disjoints
    # (en general L>=1 n'inclut pas L=0 car on stocke seulement les nouveaux).
    # Verification :
    for L in idx_levels:
        for L2 in idx_levels:
            if L < L2 and idx_levels[L] & idx_levels[L2]:
                print(f"  WARN : IDX{L} et IDX{L2} se chevauchent", flush=True)

    # Lie sequences
    seqs1 = all_lie_sequences(1)   # 82
    seqs2 = all_lie_sequences(2)   # 82 + 81*81 = 82 + 6561 = 6643
    print(f"  #seqs1 = {len(seqs1)}, #seqs2 = {len(seqs2)}", flush=True)

    # Carriers a explorer (ordre de complexite croissante).
    # On commence par les plus contraints : (4,0,0) puis on relache.
    profiles = [
        (4, 0, 0),  # tout dans IDX0 = support det3 (cas v18 actuel)
        (3, 1, 0),
        (2, 2, 0),
    ]

    candidate_idx = 0
    # Hard cap : on ne va pas exporter > 50 candidats. Au-dela, on rapporte
    # juste les statistiques.
    MAX_EXPORT = 50

    for profile in profiles:
        W = build_carrier(idx_levels, profile)
        print(f"\n=== profile {profile} : |W| = {len(W)} ===", flush=True)
        if len(W) == 0:
            continue

        # Decomposer W par poids torique pour reduire taille des blocs
        blocks = defaultdict(list)
        for tup in W:
            w = quartic_weight(tup, mons3_exp)
            blocks[w].append(tup)
        print(f"  {len(blocks)} weight blocks", flush=True)

        n_blocks_with_ker1 = 0
        n_blocks_with_ker2 = 0
        n_blocks_total = 0
        n_exported_this_profile = 0

        # Iterer par bloc, jet<=1 d'abord
        for w, Wblock in sorted(blocks.items()):
            n = len(Wblock)
            if n == 0:
                continue
            n_blocks_total += 1
            # jet<=1 mod p
            rank1, ker1 = mat_kernel_mod_p(
                build_constraints_matrix(Wblock, seqs1, gl_mats, vDet3,
                                         support_det3, PRIME),
                PRIME
            )
            dim_ker1 = len(ker1)
            if dim_ker1 == 0:
                continue   # rien a faire
            n_blocks_with_ker1 += 1

            # jet<=2 mod p (raffinement)
            rank2, ker2 = mat_kernel_mod_p(
                build_constraints_matrix(Wblock, seqs2, gl_mats, vDet3,
                                         support_det3, PRIME),
                PRIME
            )
            dim_ker2 = len(ker2)
            if dim_ker2 == 0:
                continue
            n_blocks_with_ker2 += 1

            if dim_ker2 > 5:
                # Block trop gros : on rapporte mais on ne reconstitue pas.
                if n_blocks_with_ker2 <= 5:
                    print(f"    [skip] weight {w} | n={n}, ker1={dim_ker1}, "
                          f"ker2={dim_ker2} > 5", flush=True)
                continue

            print(f"    weight {w} | n={n}, ker1={dim_ker1}, "
                  f"ker2={dim_ker2}", flush=True)

            if candidate_idx >= MAX_EXPORT:
                print(f"      cap {MAX_EXPORT} atteint : on stoppe l'export",
                      flush=True)
                break

            # Reconstruire sur QQ
            rank2qq, ker2qq = kernel_qq(Wblock, seqs2, gl_mats,
                                        vDet3, support_det3)
            print(f"      QQ : rank2 = {rank2qq}, dim ker2 = {len(ker2qq)}",
                  flush=True)
            if not ker2qq:
                continue

            for kvec in ker2qq:
                P = reconstruct_polynomial_from_kernel_vec(kvec, Wblock)
                P = normalize_to_integers(P)
                if not P:
                    continue
                fpath = export_candidate_m2(
                    P, candidate_idx, profile, w, kmax_filter=2
                )
                print(f"      EXPORTED -> {fpath} ({len(P)} monomials)",
                      flush=True)
                candidate_idx += 1
                n_exported_this_profile += 1
                if candidate_idx >= MAX_EXPORT:
                    break

        print(f"  --- profile {profile} summary : "
              f"{n_blocks_total} blocks total, "
              f"{n_blocks_with_ker1} have ker1>0, "
              f"{n_blocks_with_ker2} have ker2>0, "
              f"{n_exported_this_profile} exported", flush=True)

    print(f"\n=== TOTAL exported candidates : {candidate_idx} ===", flush=True)


if __name__ == "__main__":
    main()
