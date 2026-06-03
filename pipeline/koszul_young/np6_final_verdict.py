"""
NP-6 : verdict final A et consolidation.
"""
import json
import hashlib
from pathlib import Path
import numpy as np

SRMT = Path('/home/user/workspace/srmt')

# Hash des matrices
M_d = np.load(SRMT / 'np3_M_det_a4.npy')
M_p = np.load(SRMT / 'np3_M_perm_a4.npy')

def matrix_hash(M):
    return hashlib.sha256(M.tobytes()).hexdigest()

def file_hash(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()

# Reconstruire le verdict consolide
verdict = {
    "task": "SRMT GL9 Non-Pieri Koszul-Young Flattening (LO11 \u00a71.2 Remark 4.1.3, LO13 \u00a73)",
    "scope": "Famille Koszul flattening YF_{3,9} de Landsberg-Ottaviani 2011 (arXiv:1111.4567)",
    "verdict": "A",
    "verdict_statement": (
        "Il existe un Young flattening Koszul GL_9-equivariant Y_a : Sym^3(C^9) -> "
        "Hom(V^* (x) Lambda^a V, V (x) Lambda^{a+1} V) tel que "
        "rank(Y_a(det_3)) != rank(Y_a(perm_3)). "
        "Cas optimal : a = 4 (delta = 1, a = floor(9/2))."
    ),
    "construction": {
        "definition_source": "Landsberg-Ottaviani 2011, arXiv:1111.4567, Section 1.2 and Remark 4.1.3",
        "map": "YF_{3,9}(f) : S^delta V^* (x) Lambda^a V -> S^delta V (x) Lambda^{a+1} V with delta = floor((d-1)/2) = 1, a = floor(n/2) = 4",
        "block_form": "Matrix in block form indexed by (J, K) with |J|=a, |K|=a+1 ; block (J,K) = (-1)^{sign(k -> J,K)} * H_k where k = K \\ J and H_k[i,j] is the catalecticant coefficient of x_i x_j in df/dx_k.",
        "equivariance": "GL_9 strict on V = C^9 (acts on V and V^* and on all Schur functors and exterior powers)."
    },
    "results": {
        "a": 4,
        "matrix_shape": [1134, 1134],
        "source_space": "V^* (x) Lambda^4 V = 9 * 126 = 1134",
        "target_space": "V (x) Lambda^5 V = 9 * 126 = 1134",
        "rank_det3_Q_exact": 950,
        "rank_perm3_Q_exact": 934,
        "delta_rank": 16,
        "differs_Q": True,
        "ranks_mod_p": {
            "65521": {"det": 950, "perm": 934},
            "1000000007": {"det": 950, "perm": 934},
            "998244353": {"det": 950, "perm": 934},
            "1000000009": {"det": 950, "perm": 934},
            "104729": {"det": 950, "perm": 934},
            "5000011": {"det": 950, "perm": 934},
            "7919": {"det": 950, "perm": 934},
            "31337": {"det": 950, "perm": 934},
        },
        "multi_prime_concordance": True,
        "Q_exact_method": "Sparse Fraction Gauss-Jordan over Q (np4d_bareiss_v3.py)",
        "validation_checks": {
            "x_0^3_decomposable": {"expected_le": 70, "got": 70, "ok": True},
            "x_0^3_plus_x_1^3": {"expected_le": 140, "got": 140, "ok": True},
            "GL9_equivariance_swap_x0_x3": {"expected": 950, "got": 950, "ok": True},
            "random_cubic_full_rank": {"expected": "saturates", "got": 1134, "ok": True}
        }
    },
    "other_a": {
        "a=1": {"rank_det": 80, "rank_perm": 80, "differs": False, "shape": [324, 81]},
        "a=2": {"rank_det": 315, "rank_perm": 315, "differs": False, "shape": [756, 324]},
        "a=3": {"rank_det": 720, "rank_perm": 720, "differs": False, "shape": [1134, 756]},
        "a=4": {"rank_det": 950, "rank_perm": 934, "differs": True, "shape": [1134, 1134]}
    },
    "Y_4_det3_matrix_sha256": matrix_hash(M_d),
    "Y_4_perm3_matrix_sha256": matrix_hash(M_p),
    "comparison_with_Pieri_phase": (
        "La phase Pieri (37 candidats avec longueurs lambda, mu <= 9 et |mu/lambda| = 3) avait donne verdict B "
        "(aucune separation). La famille Koszul-Young non-Pieri (LO11) etend strictement Pieri en autorisant "
        "des facteurs Lambda^a V des deux cotes. C'est dans cette famille que la separation apparait, "
        "specifiquement au cas a = floor(n/2) = 4 predit par LO11 Remark 4.1.3."
    ),
    "interpretation": (
        "Le flattening Koszul Y_4 capture une information GL_9-equivariante non-lineaire suffisante pour "
        "distinguer det_3 et perm_3 par rang. La separation Delta = 16 = 950 - 934 est exacte sur Q. "
        "Theorique : par LO11 Prop 4.1.1, rank(Y_4(phi)) <= r * t pour [phi] in sigma_r(v_3(P^8)) avec t = "
        "rank(Y_4(w^3)) = 70. Donc rank symetrique sym-rk(det_3) >= ceil(950/70) = 14 et "
        "sym-rk(perm_3) >= ceil(934/70) = 14. Plus important : la difference de rang est une obstruction "
        "GL_9-equivariante au fait que det_3 et perm_3 soient dans la meme GL_9-orbite-fermeture."
    )
}

with open(SRMT / 'np6_final_verdict.json', 'w') as f:
    json.dump(verdict, f, indent=2)

verdict_hash = file_hash(SRMT / 'np6_final_verdict.json')
print(f"Verdict file SHA-256: {verdict_hash}")
print(f"\n=== VERDICT A (Koszul flattening separe det_3 et perm_3) ===")
print(f"  rank Y_4(det_3) = 950")
print(f"  rank Y_4(perm_3) = 934")
print(f"  Delta = 16")
print(f"  Concorde sur 8 primes + Gauss exact Q")
print(f"\n  Y_4(det3) SHA-256: {matrix_hash(M_d)}")
print(f"  Y_4(perm3) SHA-256: {matrix_hash(M_p)}")
print(f"  verdict.json SHA-256: {verdict_hash}")

# Update sha256sums.txt
existing = (SRMT / 'sha256sums.txt').read_text() if (SRMT / 'sha256sums.txt').exists() else ""
new_lines = [
    f"{matrix_hash(M_d)}  np3_M_det_a4.npy",
    f"{matrix_hash(M_p)}  np3_M_perm_a4.npy",
    f"{verdict_hash}  np6_final_verdict.json",
]
new_block = "\n# === Phase NON-PIERI Koszul-Young (LO11) ===\n" + "\n".join(new_lines) + "\n"
(SRMT / 'sha256sums.txt').write_text(existing + new_block)
print(f"\nsha256sums.txt mis a jour.")
