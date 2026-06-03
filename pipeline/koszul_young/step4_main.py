"""
Etape 4 - Boucle principale.

Enumere tous les candidats (lambda, mu) avec:
  - |lambda| = a in {0, 1, 2, 3, 4}
  - |mu| = a + 3
  - mu/lambda est une bande horizontale (Pieri)
  - len(lambda) <= 9, len(mu) <= 9 (pour V = C^9)
  - dim(S_lambda(C^9)) * dim(S_mu(C^9)) <= 1e6

Pour chaque candidat:
  1. Construit Gamma(det3), Gamma(perm3).
  2. Calcule rank mod p (p = 65537) pour les deux.
  3. Si differents -> escalade Bareiss exact Q, puis certif 3-primes.
  4. Sinon, archive et continue.

Output: results.json avec tous les candidats et leur rangs.
"""
import json, sys, time, os, hashlib
sys.path.insert(0, '/home/user/workspace/srmt')
from step2_schur import partitions, dim_Schur, SSYT
from step3_pieri import horizontal_strips
from step3b_gamma import build_gamma_matrix, rank_mod_p, rank_int_matrix_bareiss

# Lecture det3, perm3
with open('/home/user/workspace/srmt/step1_detperm.json') as F:
    artifact = json.load(F)
det3 = artifact['det3']
perm3 = artifact['perm3']

# Premiers pour certif multi-prime
PRIMES_CERTIF = [1000000007, 998244353, 1000000009]
P_FAST = 65537

DIM_LIMIT = 1_000_000  # dim_lam * dim_mu <= 10^6
NMAX = 9
A_MAX = 4

CHECKPOINT = '/home/user/workspace/srmt/results.json'
LOG = '/home/user/workspace/srmt/main_log.txt'

def log(msg):
    with open(LOG, 'a') as f: f.write(msg + '\n')
    print(msg, flush=True)

# Liste candidats
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
                    'feasible': prod <= DIM_LIMIT
                })
    return cands

def candidate_key(c):
    return f"a={c['a']} lam={tuple(c['lam'])} mu={tuple(c['mu'])} d_lam={c['d_lam']} d_mu={c['d_mu']} prod={c['prod']:,}"

def process_candidate(cand):
    """Calcule rank modulaire puis (si different) Bareiss + certif."""
    lam = tuple(cand['lam']); mu = tuple(cand['mu'])
    t0 = time.time()
    try:
        M_det, _, _ = build_gamma_matrix(lam, mu, det3, n=NMAX)
        M_perm, _, _ = build_gamma_matrix(lam, mu, perm3, n=NMAX)
    except Exception as e:
        cand['error'] = f'build: {e}'
        return cand
    t_build = time.time() - t0
    
    t1 = time.time()
    r_det_p = rank_mod_p(M_det, P_FAST)
    r_perm_p = rank_mod_p(M_perm, P_FAST)
    t_modular = time.time() - t1
    cand['rank_det_modp65537'] = r_det_p
    cand['rank_perm_modp65537'] = r_perm_p
    cand['t_build_s'] = round(t_build, 2)
    cand['t_modular_s'] = round(t_modular, 2)
    
    if r_det_p == r_perm_p:
        cand['status'] = 'EQUAL_MODULAR'
        return cand
    
    # ESCALADE: differents en modulaire -> Bareiss exact
    log(f"  !! Difference modulaire! lam={lam} mu={mu}: det={r_det_p} perm={r_perm_p}. Escalade Bareiss.")
    t2 = time.time()
    r_det_Q = rank_int_matrix_bareiss(M_det)
    r_perm_Q = rank_int_matrix_bareiss(M_perm)
    t_bareiss = time.time() - t2
    cand['rank_det_Q'] = r_det_Q
    cand['rank_perm_Q'] = r_perm_Q
    cand['t_bareiss_s'] = round(t_bareiss, 2)
    
    if r_det_Q == r_perm_Q:
        cand['status'] = 'EQUAL_Q_MOD_NOISE'
        return cand
    
    # SUCCESS: differents sur Q. Certifier 3-primes
    log(f"  ** Difference Q confirmee. Certif multi-prime.")
    certif = {}
    for p in PRIMES_CERTIF:
        r1 = rank_mod_p(M_det, p)
        r2 = rank_mod_p(M_perm, p)
        certif[p] = {'rank_det': r1, 'rank_perm': r2, 'diff': r1 != r2}
    cand['certif_multi_prime'] = certif
    if all(certif[p]['diff'] for p in PRIMES_CERTIF):
        cand['status'] = 'SUCCESS_CERTIFIED'
        # Sauvegarde matrices
        outdir = f"/home/user/workspace/srmt/CERTIF_lam={'_'.join(map(str,lam))}_mu={'_'.join(map(str,mu))}"
        os.makedirs(outdir, exist_ok=True)
        with open(f'{outdir}/M_det.json','w') as F: json.dump(M_det, F)
        with open(f'{outdir}/M_perm.json','w') as F: json.dump(M_perm, F)
        cand['certif_dir'] = outdir
        cand['M_det_sha256'] = hashlib.sha256(json.dumps(M_det).encode()).hexdigest()
        cand['M_perm_sha256'] = hashlib.sha256(json.dumps(M_perm).encode()).hexdigest()
    else:
        cand['status'] = 'SUCCESS_UNCERTIFIED'
    return cand


def main():
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as F: 
            data = json.load(F)
        done_keys = {f"{tuple(c['lam'])}->{tuple(c['mu'])}" for c in data['candidates']}
        log(f"Reprise: {len(done_keys)} candidats deja faits.")
    else:
        data = {'candidates': [], 'started_at': time.time()}
        done_keys = set()
    
    cands = enumerate_candidates()
    cands.sort(key=lambda c: c['prod'])  # commencer par les petits
    log(f"=== {len(cands)} candidats total ===")
    feas = [c for c in cands if c['feasible']]
    skip = [c for c in cands if not c['feasible']]
    log(f"  Faisables (prod <= {DIM_LIMIT:,}): {len(feas)}")
    log(f"  A skipper (prod > {DIM_LIMIT:,}): {len(skip)}")
    
    for cand in cands:
        key = f"{tuple(cand['lam'])}->{tuple(cand['mu'])}"
        if key in done_keys: continue
        log(f"\n[{candidate_key(cand)}]")
        if not cand['feasible']:
            cand['status'] = 'SKIPPED_SIZE'
            data['candidates'].append(cand)
            done_keys.add(key)
            with open(CHECKPOINT,'w') as F: json.dump(data, F, indent=2)
            continue
        try:
            cand = process_candidate(cand)
        except MemoryError:
            cand['status'] = 'SKIPPED_OOM'
            log("  OOM -> skipped")
        except Exception as e:
            cand['status'] = f'ERROR: {e}'
            log(f"  ERROR: {e}")
        log(f"  -> status={cand.get('status')}, r_det_p={cand.get('rank_det_modp65537')}, r_perm_p={cand.get('rank_perm_modp65537')}")
        data['candidates'].append(cand)
        done_keys.add(key)
        with open(CHECKPOINT,'w') as F: json.dump(data, F, indent=2)
        if cand.get('status') == 'SUCCESS_CERTIFIED':
            log("=== SUCCESS CERTIFIED. Stopping main loop. ===")
            data['final_verdict'] = 'A'
            with open(CHECKPOINT,'w') as F: json.dump(data, F, indent=2)
            return data
    
    data['ended_at'] = time.time()
    with open(CHECKPOINT,'w') as F: json.dump(data, F, indent=2)
    return data


if __name__ == "__main__":
    data = main()
    print("\n=== FIN ===")
    statuses = {}
    for c in data['candidates']:
        s = c.get('status','?')
        statuses[s] = statuses.get(s,0)+1
    for s,n in sorted(statuses.items()): print(f"  {s}: {n}")
