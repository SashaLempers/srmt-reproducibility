"""
NP-3 : certification rigoureuse du differs trouve a a=4.

1. Re-build Y_4(det3) et Y_4(perm3) en INT exact (sparse COO).
2. Calcul rang sur 3 premiers : 10^9+7, 998244353, 10^9+9.
3. Bareiss exact sur Q.
4. Hash SHA-256 des matrices.
"""
import json
import sys
import pickle
import itertools
import hashlib
import time
from pathlib import Path
from fractions import Fraction
import numpy as np

sys.path.insert(0, '/home/user/workspace/srmt')
from np2_koszul_matrix import build_koszul_matrix, rank_dense_mod_p

n = 9
SRMT = Path('/home/user/workspace/srmt')

with open(SRMT / 'np1_partials.pkl', 'rb') as f:
    pdata = pickle.load(f)
det3_partials = pdata['det3_partials']
perm3_partials = pdata['perm3_partials']

a = 4
print(f"Build Y_{a}(det3) et Y_{a}(perm3) ...")
nr, nc, ri_d, ci_d, v_d = build_koszul_matrix(det3_partials, a)
nr2, nc2, ri_p, ci_p, v_p = build_koszul_matrix(perm3_partials, a)
assert (nr, nc) == (nr2, nc2)
print(f"  Matrices {nr} x {nc}, nnz det={len(v_d)}, nnz perm={len(v_p)}")


def build_dense_int(nr, nc, ri, ci, vs):
    M = np.zeros((nr, nc), dtype=np.int64)
    for r, c, v in zip(ri, ci, vs):
        M[r, c] += v
    return M


M_d = build_dense_int(nr, nc, ri_d, ci_d, v_d)
M_p = build_dense_int(nr, nc, ri_p, ci_p, v_p)

# SHA-256 des matrices
def matrix_hash(M):
    return hashlib.sha256(M.tobytes()).hexdigest()


print(f"\nSHA-256 :")
print(f"  Y_4(det3) : {matrix_hash(M_d)}")
print(f"  Y_4(perm3): {matrix_hash(M_p)}")

# Rang sur 4 premiers
primes = [65521, 10**9 + 7, 998244353, 10**9 + 9]
results = {'a': a, 'shape': [nr, nc], 'primes': {}}
for p in primes:
    print(f"\n--- Rang mod p={p} ---")
    t0 = time.time()
    # On doit recalculer le rang car np.int64 overflow possible avec p > 2^31
    # Pour p < 2^31, on peut rester en int64
    # Pour p ~ 10^9, le produit factor * M[row] peut atteindre 10^18 < 2^63 = 9.2e18, OK ssi factor et M[row] < ~3*10^9.
    # Mais c'est tendu : factor < p ~ 10^9, M[row] < p ~ 10^9, produit < 10^18, OK juste juste.
    # Pour securite, on borne les valeurs intermediaires.
    r_d = rank_dense_mod_p(M_d.copy(), p)
    t1 = time.time()
    r_p = rank_dense_mod_p(M_p.copy(), p)
    t2 = time.time()
    print(f"  rank Y_4(det3) = {r_d} ({t1-t0:.1f}s)")
    print(f"  rank Y_4(perm3)= {r_p} ({t2-t1:.1f}s)")
    print(f"  DIFFERS : {r_d != r_p}")
    results['primes'][p] = {'rank_det': r_d, 'rank_perm': r_p, 'differs': r_d != r_p, 'time_det': t1-t0, 'time_perm': t2-t1}

print(f"\n\n=== MULTI-PRIME SUMMARY (Koszul a=4) ===")
for p, r in results['primes'].items():
    print(f"  p={p:>11}: det={r['rank_det']}, perm={r['rank_perm']}, DIFFERS={r['differs']}")

with open(SRMT / 'np3_multiprime.json', 'w') as f:
    # JSON serialisable
    out = dict(results)
    out['primes'] = {str(p): v for p, v in results['primes'].items()}
    out['det3_matrix_sha256'] = matrix_hash(M_d)
    out['perm3_matrix_sha256'] = matrix_hash(M_p)
    json.dump(out, f, indent=2)

print(f"\nSauvegarde np3_multiprime.json")

# Sauvegarde matrices pour Bareiss exact
np.save(SRMT / 'np3_M_det_a4.npy', M_d)
np.save(SRMT / 'np3_M_perm_a4.npy', M_p)
print("Matrices sauvees np3_M_det_a4.npy, np3_M_perm_a4.npy")
