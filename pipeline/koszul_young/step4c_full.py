"""
Boucle principale FAST: utilise build sparse + dense numpy mod p.
Traite TOUS les candidats faisables avec dim <= 5*10^7 (limite RAM ~400 MB).
"""
import json, sys, time, os, hashlib
import numpy as np
sys.path.insert(0, '/home/user/workspace/srmt')
from step2_schur import partitions, dim_Schur
from step3_pieri import horizontal_strips
from step4b_fast_modular import build_gamma_sparse_arrays, rank_mod_p_dense_numpy, P

with open('/home/user/workspace/srmt/step1_detperm.json') as F:
    artifact = json.load(F)
det3 = artifact['det3']
perm3 = artifact['perm3']

PRIMES_CERTIF = [1000000007, 998244353, 1000000009]
NMAX = 9
A_MAX = 4
DIM_LIMIT = 5 * 10**7
DENSE_RAM_LIMIT = 5 * 10**7  # cellules

CHECKPOINT = '/home/user/workspace/srmt/results_full.json'
LOG = '/home/user/workspace/srmt/main_full_log.txt'

def log(msg):
    with open(LOG, 'a') as f: f.write(msg + '\n')
    print(msg, flush=True)

def rank_mod_int_p(rows, cols, vals, R, C, p):
    """Variante: convertit vals -> int Python pour gros premiers."""
    A = np.zeros((R, C), dtype=np.int64)
    vals_mod = (vals % p).astype(np.int64)
    np.add.at(A, (rows, cols), vals_mod)
    A %= p
    rank = 0; row = 0
    for col in range(C):
        if row >= R: break
        piv = -1
        for r in range(row, R):
            if A[r, col] != 0: piv = r; break
        if piv < 0: continue
        if piv != row:
            A[[row, piv]] = A[[piv, row]]
        inv = pow(int(A[row, col]), p-2, p)
        for r in range(R):
            if r == row: continue
            if A[r, col] != 0:
                factor = int((A[r, col] * inv) % p)
                A[r] = (A[r] - factor * A[row]) % p
        row += 1; rank += 1
    return rank

def enumerate_candidates():
    cands = []
    for a in range(0, A_MAX + 1):
        for lam in partitions(a, max_len=NMAX):
            d_lam = dim_Schur(lam, NMAX) if lam else 1
            if d_lam == 0: continue
            for mu in horizontal_strips(lam, 3, NMAX):
                if len(mu) > NMAX: continue
                d_mu = dim_Schur(mu, NMAX)
                if d_mu == 0: continue
                prod = d_lam * d_mu
                cands.append({
                    'a': a,
                    'lam': list(lam),
                    'mu': list(mu),
                    'd_lam': d_lam,
                    'd_mu': d_mu,
                    'prod': prod,
                })
    return cands

def process(cand):
    lam = tuple(cand['lam']); mu = tuple(cand['mu'])
    prod = cand['prod']
    t0 = time.time()
    try:
        rd = build_gamma_sparse_arrays(lam, mu, det3, n=NMAX)
        rp = build_gamma_sparse_arrays(lam, mu, perm3, n=NMAX)
    except Exception as e:
        cand['error'] = f'build: {e}'
        cand['status'] = 'BUILD_ERROR'
        return cand
    cand['nnz_det'] = int(len(rd[2]))
    cand['nnz_perm'] = int(len(rp[2]))
    cand['shape'] = [int(rd[3]), int(rd[4])]
    cand['t_build_s'] = round(time.time()-t0, 3)
    
    if prod > DENSE_RAM_LIMIT:
        cand['status'] = 'SKIPPED_RAM'
        return cand
    
    t1 = time.time()
    r_det = rank_mod_p_dense_numpy(*rd, p=P)
    r_perm = rank_mod_p_dense_numpy(*rp, p=P)
    cand['rank_det_p1'] = int(r_det)
    cand['rank_perm_p1'] = int(r_perm)
    cand['t_rank_s'] = round(time.time()-t1, 3)
    
    if r_det == r_perm:
        cand['status'] = 'EQUAL_MODULAR'
        return cand
    
    # Difference -> certifier sur 3 primes
    log(f"  !!! Difference detectee. Certif 3-primes.")
    certif = {}
    for p in PRIMES_CERTIF:
        # Pour ces primes (10^9 +), Python int natif via boucle, mais on peut quand meme utiliser numpy int64
        rd2 = build_gamma_sparse_arrays(lam, mu, det3, n=NMAX)
        rp2 = build_gamma_sparse_arrays(lam, mu, perm3, n=NMAX)
        rd_p = rank_mod_int_p(*rd2, p=p)
        rp_p = rank_mod_int_p(*rp2, p=p)
        certif[p] = {'rank_det': int(rd_p), 'rank_perm': int(rp_p), 'diff': rd_p != rp_p}
        log(f"    p={p}: det={rd_p} perm={rp_p}")
    cand['certif_3primes'] = certif
    if all(certif[p]['diff'] for p in PRIMES_CERTIF):
        cand['status'] = 'SUCCESS_CERTIFIED'
    else:
        cand['status'] = 'INCONSISTENT'
    return cand

def main():
    cands = enumerate_candidates()
    cands.sort(key=lambda c: c['prod'])
    log(f"=== {len(cands)} candidats ===")
    data = {'candidates': [], 'started_at': time.time(), 'p_fast': P, 'primes_certif': PRIMES_CERTIF}
    for cand in cands:
        log(f"\n[a={cand['a']} lam={cand['lam']} mu={cand['mu']} prod={cand['prod']:,}]")
        try:
            cand = process(cand)
        except MemoryError:
            cand['status'] = 'OOM'
            log("  OOM")
        except Exception as e:
            cand['status'] = f'ERR: {e}'
            log(f"  ERR: {e}")
        log(f"  status={cand.get('status')} r_det={cand.get('rank_det_p1')} r_perm={cand.get('rank_perm_p1')} t_rank={cand.get('t_rank_s','?')}s")
        data['candidates'].append(cand)
        with open(CHECKPOINT,'w') as F: json.dump(data, F, indent=2)
        if cand.get('status') == 'SUCCESS_CERTIFIED':
            log("=== SUCCESS CERTIFIED ===")
            data['final'] = 'A'
            with open(CHECKPOINT,'w') as F: json.dump(data, F, indent=2)
            return data
    data['ended_at'] = time.time()
    with open(CHECKPOINT,'w') as F: json.dump(data, F, indent=2)
    return data

if __name__ == "__main__":
    if os.path.exists(LOG): os.remove(LOG)
    if os.path.exists(CHECKPOINT): os.remove(CHECKPOINT)
    data = main()
    print("\n=== FIN ===")
    statuses = {}
    for c in data['candidates']:
        s = c.get('status','?'); statuses[s] = statuses.get(s,0)+1
    for s,n in sorted(statuses.items()): print(f"  {s}: {n}")
