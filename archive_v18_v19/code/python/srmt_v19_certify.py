#!/usr/bin/env python3
# =====================================================================
# srmt_v19_certify.py
#
# Certifie en Python (QQ exact, sparse) que P passe testJet(P, kmax).
# Reutilise les meme briques que srmt_v19_sparse_search.py :
#   - derive_quartic (Leibniz sparse via gl_mats)
#   - eval_at_vDet3
#   - pruning par budget (strict, identique au pipeline de filtrage)
# Pour chaque candidate_NNN.m2 :
#   1. Parse Pcandidate = sum coef * cVars#i * cVars#j * cVars#k * cVars#l ;
#   2. Pour k=0,1,...,kmax, parcourt toutes les seq Lie de longueur k et
#      verifie (D^seq P)(v_det3) == 0.
#   3. Sortie : (ok, level_fail) ou ok=True si passe, sinon level_fail.
#
# Optim : on travaille sur P entier (pas un carrier de monomes) et on prune
# AVANT chaque derivation par un budget r_remaining = kmax - step - 1.
#
# Statuts (alignes editorial v18) :
#   QQ-CERTIFIED     : passe testJet(P, kmax) sur QQ exact, kmax >= 4.
#   Computational obs.: kmax < 4.
#   Rejected         : level_fail < kmax atteint.
# =====================================================================

import sys
import os
import re
import time
from fractions import Fraction
from itertools import product

# Reutilise les routines du module sparse search
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srmt_v19_sparse_search import (
    load_m2_data,
    derive_quartic,
    eval_at_vDet3,
    prune_quartic,
)


def parse_candidate_m2(path):
    """
    Parse un fichier candidate_NNN.m2 et retourne dict{tuple_trie: int}.
    Format attendu (parens, signes explicites) :
        Pcandidate = (
               1 * cVars#34 * cVars#72 * cVars#92 * cVars#153
             + 1 * cVars#37 * cVars#65 * cVars#96 * cVars#153
             );
    Tolere : signes - ; coefs entiers.
    """
    with open(path, "r") as f:
        text = f.read()
    # Extrait le bloc entre Pcandidate = ( ... );
    m = re.search(r"Pcandidate\s*=\s*\((.+?)\)\s*;", text, re.DOTALL)
    if not m:
        raise ValueError(f"Pcandidate block not found in {path}")
    body = m.group(1)
    # Une ligne par terme
    P = {}
    for line in body.split("\n"):
        line = line.strip()
        if not line:
            continue
        # extrait signe, coef, indices
        # ex: "1 * cVars#34 * cVars#72 * cVars#92 * cVars#153"
        # ex: "+ 1 * cVars#34 * ..."  ou "- 2 * cVars#..."
        sign = 1
        if line.startswith("+"):
            line = line[1:].strip()
        elif line.startswith("-"):
            sign = -1
            line = line[1:].strip()
        # split sur "*"
        parts = [p.strip() for p in line.split("*")]
        if not parts:
            continue
        # premier = coef entier, suivants = cVars#NNN
        try:
            coef = int(parts[0])
        except ValueError:
            raise ValueError(f"Cannot parse coef in line: {line}")
        idxs = []
        for p in parts[1:]:
            mm = re.match(r"cVars#(\d+)", p)
            if not mm:
                raise ValueError(f"Cannot parse cVars in line: {line}")
            idxs.append(int(mm.group(1)))
        if len(idxs) != 4:
            raise ValueError(f"Expected 4 cVars, got {len(idxs)} in: {line}")
        tup = tuple(sorted(idxs))
        c = sign * coef
        P[tup] = P.get(tup, 0) + c
        if P[tup] == 0:
            del P[tup]
    return P


def derive_polynomial(P, ab, gl_mats, mod=None):
    """Derivation D_ab d'un dict-polynome (4 facteurs, sparse)."""
    return derive_quartic(P, ab, gl_mats, mod=mod)


def test_jet(P, kmax, gl_mats, vDet3, support_det3, verbose=False):
    """
    Verifie (D^seq P)(v_det3) == 0 pour toute seq de longueur 0..kmax.
    Retourne (ok, level_fail).
    On itere par couches : layer_k = liste des polynomes obtenus par
    application de seq de longueur k. Mais on ne stocke pas tous ; on
    parcourt en DFS avec pruning early-exit.

    Optim majeure : pruning AVANT chaque derivation residuelle avec
    budget = kmax - k_courant. NB : pruning supprime les monomes qui ne
    pourront jamais entrer dans le support apres derivees restantes,
    MAIS on cherche a verifier que P(v_det3) ET (D P)(v_det3) ET ...
    sont TOUS zero, donc le pruning est correct a chaque depth k :
    quand on est en depth k et qu'on doit verifier (D^seq P)(v_det3) avec
    seq de longueur k, c'est r_remaining = 0 ; pour les depths suivantes
    on derive encore et on prune.
    """
    pairs = [(a, b) for a in range(9) for b in range(9)]

    # Verif depth 0
    v0 = eval_at_vDet3(P, vDet3, support_det3)
    if v0 != 0:
        if verbose:
            print(f"    fail at depth 0 : {v0}")
        return (False, 0)

    if kmax == 0:
        return (True, 0)

    # DFS ; current = liste de (P_courant, depth_remaining)
    # Plus simple : BFS layer-by-layer.
    layer = [P]
    for depth in range(1, kmax + 1):
        new_layer = []
        n_processed = 0
        for Q in layer:
            for ab in pairs:
                # Optim : prune avant derive avec budget kmax - depth (deriv restantes apres celle-ci)
                # Non : on derive maintenant a depth, donc apres cette derive il reste 0 derives.
                # Mais on veut evaluer EN v_det3 a depth=current. Donc:
                #   apres cette derive, eval direct -> r_remaining = 0
                # Donc on prune Q avant derive avec r_remaining = 1 (1 derive restante = celle-ci).
                # Apres derive, on prune avec 0.
                Qp = prune_quartic(Q, 1, support_det3)
                if not Qp:
                    continue
                Q2 = derive_polynomial(Qp, ab, gl_mats, mod=None)
                if not Q2:
                    continue
                Q2 = prune_quartic(Q2, 0, support_det3)
                if not Q2:
                    continue
                v = eval_at_vDet3(Q2, vDet3, support_det3)
                if v != 0:
                    if verbose:
                        print(f"    fail at depth {depth}, seq=...{ab}, val={v}")
                    return (False, depth)
                # Pour le prochain layer, on a besoin du Q2 derivee NON evaluee (mais avec un budget plus large)
                # On stocke Q2 sans pruning agressif (avec r_remaining = kmax - depth)
                if depth < kmax:
                    Q2_full = derive_polynomial(Qp, ab, gl_mats, mod=None)
                    Q2_full = prune_quartic(Q2_full, kmax - depth, support_det3)
                    if Q2_full:
                        new_layer.append(Q2_full)
                n_processed += 1
        if verbose:
            print(f"  depth {depth} done : {n_processed} (a,b) processed, "
                  f"new_layer size = {len(new_layer)}")
        layer = new_layer
        if not layer:
            # plus rien a deriver, tous les chemins sont eteints par pruning
            if verbose:
                print(f"  layer empty after depth {depth}, all paths zero")
            return (True, kmax)

    return (True, kmax)


def test_jet_recursive(P, kmax, gl_mats, vDet3, support_det3, verbose=False):
    """
    Version DFS recursive plus efficace en memoire : on n'enumère qu'une
    branche a la fois, avec pruning a chaque etape par budget.

    invariant : on appelle dfs(P_local, depth, budget_remaining=kmax-depth)
      - eval P_local en v_det3, doit etre 0
      - si budget==0, retour OK
      - sinon, pour chaque ab in 81 pairs :
          P' = prune(P_local, budget, support); if empty -> branche vide
          P'' = derive(P', ab); if empty -> branche vide
          dfs(P'', depth+1, budget-1)
    """
    pairs = [(a, b) for a in range(9) for b in range(9)]

    # statistiques
    counter = [0, 0]  # [evals, derivs]

    def dfs(Pcur, depth, budget):
        counter[0] += 1
        v = eval_at_vDet3(Pcur, vDet3, support_det3)
        if v != 0:
            return (False, depth)
        if budget == 0:
            return (True, depth)
        # prune Pcur with budget : monomes avec >budget indices hors support
        # ne contribueront jamais.
        Pp = prune_quartic(Pcur, budget, support_det3)
        if not Pp:
            return (True, depth)  # tout est nul a partir d'ici
        for ab in pairs:
            Q = derive_polynomial(Pp, ab, gl_mats, mod=None)
            if not Q:
                continue
            res = dfs(Q, depth + 1, budget - 1)
            if not res[0]:
                return res
        return (True, kmax)

    res = dfs(P, 0, kmax)
    if verbose:
        print(f"    [{counter[0]} evals]")
    return res


def certify_candidate(path, kmax, gl_mats, vDet3, support_det3, verbose=True):
    P = parse_candidate_m2(path)
    name = os.path.basename(path)
    n_terms = len(P)
    if verbose:
        print(f"  {name} : {n_terms} terms", flush=True)
    t1 = time.time()
    ok, level = test_jet_recursive(P, kmax, gl_mats, vDet3, support_det3,
                                   verbose=False)
    dt = time.time() - t1
    status = "QQ-CERTIFIED" if ok and kmax >= 4 else \
             ("Computational-obs" if ok else f"Rejected@{level}")
    if verbose:
        print(f"    -> testJet(P,{kmax}) = ({ok}, {level})   "
              f"[{dt:.2f}s]   {status}", flush=True)
    return (ok, level, dt, status)


def main():
    if len(sys.argv) < 2:
        print("usage : srmt_v19_certify.py <kmax> [n_candidates]", flush=True)
        sys.exit(1)
    kmax = int(sys.argv[1])
    n_max = int(sys.argv[2]) if len(sys.argv) >= 3 else None

    print(f"=== SRMT v19 certification (Python QQ) kmax={kmax} ===", flush=True)
    print("  loading M2 data...", flush=True)
    mons3_exp, gl_mats, vDet3, support_det3 = load_m2_data()
    print(f"  N3=165, |support_det3|={len(support_det3)}", flush=True)

    cdir = "./candidates"
    files = sorted([os.path.join(cdir, f) for f in os.listdir(cdir)
                    if f.startswith("candidate_") and f.endswith(".m2")])
    if n_max is not None:
        files = files[:n_max]

    print(f"  certifying {len(files)} candidates at kmax={kmax}", flush=True)

    results = []
    for f in files:
        ok, level, dt, status = certify_candidate(f, kmax, gl_mats, vDet3,
                                                  support_det3, verbose=True)
        results.append((f, ok, level, dt, status))

    print("\n=== SUMMARY ===", flush=True)
    n_ok = sum(1 for r in results if r[1])
    n_ko = len(results) - n_ok
    print(f"  ok = {n_ok} / {len(results)}, rejected = {n_ko}", flush=True)
    rejection_levels = {}
    for _, ok, lvl, _, _ in results:
        if not ok:
            rejection_levels[lvl] = rejection_levels.get(lvl, 0) + 1
    if rejection_levels:
        print(f"  rejection levels : {sorted(rejection_levels.items())}",
              flush=True)

    # save results
    out = "./results/certify_kmax%d.txt" % kmax
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as g:
        g.write(f"kmax = {kmax}\n")
        for f, ok, lvl, dt, st in results:
            g.write(f"{os.path.basename(f)}  ok={ok}  level={lvl}  "
                    f"dt={dt:.2f}s  status={st}\n")
    print(f"  wrote {out}", flush=True)


if __name__ == "__main__":
    main()
