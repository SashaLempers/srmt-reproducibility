#!/usr/bin/env python3
# =====================================================================
# srmt_v19_full.py
#
# Pipeline complet : pour chaque profil dans {(4,0,0),(3,1,0),(2,2,0)},
# pour chaque bloc poids torique :
#   1. ker1 mod p ; si 0, skip
#   2. ker2 mod p ; si 0, skip
#   3. ker2 sur QQ exact
#   4. Pour chaque kvec : reconstruire P, normaliser, certifier kmax=4
#      en Python QQ exact (testJet sparse + pruning).
#   5. N'exporter en .m2 QUE les SURVIVANTS kmax=4.
#
# Sortie :
#   ./results/srmt_v19_full.txt : rapport
#   ./results/survivors/         : .m2 survivants
# =====================================================================

import sys, os, time
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from srmt_v19_sparse_search import (
    load_m2_data, all_lie_sequences, build_carrier, support_levels,
    build_constraints_matrix, mat_kernel_mod_p,
    kernel_qq, reconstruct_polynomial_from_kernel_vec,
    normalize_to_integers, export_candidate_m2,
    quartic_weight, PRIME,
)
from srmt_v19_certify import test_jet_recursive


def main():
    KMAX = 4
    LOG = "./results/srmt_v19_full.txt"
    SURVIVE_DIR = "./results/survivors"
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    os.makedirs(SURVIVE_DIR, exist_ok=True)

    log_f = open(LOG, "w")
    def log(msg):
        print(msg, flush=True)
        log_f.write(msg + "\n")
        log_f.flush()

    log(f"=== srmt_v19_full.py (kmax={KMAX} certif inline) ===")
    mons3_exp, gl_mats, vDet3, support_det3 = load_m2_data()
    log(f"  N3=165, |support_det3|={len(support_det3)}")

    idx_levels = support_levels(mons3_exp, gl_mats, support_det3, max_level=2)
    for L in sorted(idx_levels):
        log(f"  IDX{L} = {len(idx_levels[L])} indices")

    seqs1 = all_lie_sequences(1)
    seqs2 = all_lie_sequences(2)
    log(f"  #seqs1 = {len(seqs1)}, #seqs2 = {len(seqs2)}")

    profiles = [(4,0,0), (3,1,0), (2,2,0)]

    survivors = []        # liste de (P, profile, weight)
    candidate_idx = 0     # index global pour fichiers .m2

    for profile in profiles:
        W = build_carrier(idx_levels, profile)
        log(f"\n=== profile {profile} : |W| = {len(W)} ===")
        if not W:
            continue

        blocks = defaultdict(list)
        for tup in W:
            w = quartic_weight(tup, mons3_exp)
            blocks[w].append(tup)
        log(f"  {len(blocks)} weight blocks")

        n_total = 0
        n_ker1 = 0
        n_ker2 = 0
        n_ker2_skipped_big = 0
        n_qq_zero = 0
        n_reconstructed = 0
        n_certified_kmax4 = 0
        n_rejected = defaultdict(int)
        t_profile_start = time.time()

        for w, Wblock in sorted(blocks.items()):
            n = len(Wblock)
            if n == 0:
                continue
            n_total += 1

            # jet1
            rows1 = build_constraints_matrix(Wblock, seqs1, gl_mats, vDet3,
                                             support_det3, PRIME)
            _, ker1 = mat_kernel_mod_p(rows1, PRIME)
            if not ker1:
                continue
            n_ker1 += 1

            # jet2
            rows2 = build_constraints_matrix(Wblock, seqs2, gl_mats, vDet3,
                                             support_det3, PRIME)
            _, ker2 = mat_kernel_mod_p(rows2, PRIME)
            if not ker2:
                continue
            n_ker2 += 1

            if len(ker2) > 5:
                n_ker2_skipped_big += 1
                continue

            # QQ kernel
            _, ker2qq = kernel_qq(Wblock, seqs2, gl_mats, vDet3, support_det3)
            if not ker2qq:
                n_qq_zero += 1
                continue

            for kvec in ker2qq:
                P = reconstruct_polynomial_from_kernel_vec(kvec, Wblock)
                P = normalize_to_integers(P)
                if not P:
                    continue
                n_reconstructed += 1
                # Certif inline kmax=KMAX
                ok, level = test_jet_recursive(P, KMAX, gl_mats, vDet3,
                                               support_det3)
                if ok:
                    n_certified_kmax4 += 1
                    survivors.append((P, profile, w))
                    fpath = os.path.join(
                        SURVIVE_DIR, f"survivor_{candidate_idx:03d}.m2")
                    # ecrit a la main car export_candidate_m2 utilise CANDIDATES_DIR
                    write_m2(fpath, P, candidate_idx, profile, w, KMAX)
                    candidate_idx += 1
                    log(f"  *** SURVIVOR weight {w} : {len(P)} terms -> {fpath}")
                else:
                    n_rejected[level] += 1

        dt = time.time() - t_profile_start
        log(f"  --- profile {profile} : {n_total} blocks, "
            f"{n_ker1} ker1>0, {n_ker2} ker2>0, "
            f"{n_ker2_skipped_big} ker2>5 (skipped), "
            f"{n_qq_zero} ker2 mod p but ker2 QQ = 0, "
            f"{n_reconstructed} reconstructed, "
            f"{n_certified_kmax4} CERTIFIED kmax={KMAX}, "
            f"rejected by level: {dict(n_rejected)}  "
            f"[{dt:.1f}s]")

    log(f"\n=== TOTAL : {len(survivors)} survivors at kmax={KMAX} ===")
    for i, (P, prof, w) in enumerate(survivors):
        log(f"  survivor {i:03d} : profile={prof}, weight={w}, terms={len(P)}")

    log_f.close()


def write_m2(fpath, P_dict, idx, profile, weight, kmax_filter):
    """Ecrit un .m2 survivant avec entete clair (statut QQ-CERTIFIED)."""
    with open(fpath, "w") as f:
        f.write(f"-- survivor_{idx:03d}.m2\n")
        f.write(f"-- profile = {profile}, weight = {weight}\n")
        f.write(f"-- STATUS : Computational obs. (passes Python testJet, "
                f"kmax={kmax_filter})\n")
        f.write(f"-- {len(P_dict)} monomials\n")
        f.write("Pcandidate = (\n")
        terms = []
        for tup, coef in sorted(P_dict.items()):
            cs = " * ".join(f"cVars#{x}" for x in tup)
            sign = "+" if coef > 0 else "-"
            terms.append(f"     {sign} {abs(coef)} * {cs}")
        if terms[0].lstrip().startswith("+"):
            terms[0] = terms[0].replace("+ ", "  ", 1)
        f.write("\n".join(terms))
        f.write("\n     );\n")


if __name__ == "__main__":
    main()
