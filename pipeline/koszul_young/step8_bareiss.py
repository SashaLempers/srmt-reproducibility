"""
Bareiss fraction-free rank pour matrices entieres.
Utilise des entiers Python natifs (precision arbitraire).
Algorithme: Bareiss two-step pivot, calcul du rang par detection des lignes nulles.

Implementation: matrice 2D list[list[int]]. Operations entieres exactes.
La cle est: A[r][c] := (A[piv][piv]*A[r][c] - A[r][piv]*A[piv][c]) / prev_pivot
ou prev_pivot est le pivot precedent (1 au debut).
"""
import json, sys, time, os, hashlib
sys.path.insert(0, '/home/user/workspace/srmt')
import numpy as np
from step3b_gamma import build_gamma_matrix
from step4b_fast_modular import build_gamma_sparse_arrays
from step2_schur import SSYT


def bareiss_rank(M):
    """Rang exact sur Q (= Z[1/p]) via Bareiss fraction-free.
    M: list[list[int]] (sera modifiee).
    Retourne le rang."""
    if not M or not M[0]: return 0
    R = len(M); C = len(M[0])
    # Conversion en liste de listes d'entiers
    A = [list(r) for r in M]
    prev = 1
    rank = 0
    pivot_row = 0
    for col in range(C):
        if pivot_row >= R: break
        # Trouve pivot non nul
        piv = -1
        for r in range(pivot_row, R):
            if A[r][col] != 0:
                piv = r; break
        if piv < 0: continue
        if piv != pivot_row:
            A[pivot_row], A[piv] = A[piv], A[pivot_row]
        pivot_val = A[pivot_row][col]
        # Mise a jour Bareiss pour les lignes != pivot_row
        for r in range(R):
            if r == pivot_row: continue
            arc = A[r][col]
            if arc == 0:
                # Encore Bareiss: A[r][c'] = (pivot_val * A[r][c'] - 0) / prev = pivot_val * A[r][c'] / prev
                # Division exacte garantie par Bareiss
                row_r = A[r]
                for c in range(col+1, C):
                    val = pivot_val * row_r[c]
                    row_r[c] = val // prev
            else:
                row_r = A[r]
                row_piv = A[pivot_row]
                for c in range(col+1, C):
                    val = pivot_val * row_r[c] - arc * row_piv[c]
                    row_r[c] = val // prev
                row_r[col] = 0
        prev = pivot_val
        pivot_row += 1
        rank += 1
    return rank


def bareiss_rank_numpy(M):
    """Version numpy avec int Python. Plus rapide grâce a vectorisation."""
    if not M or not M[0]: return 0
    # On utilise dtype=object pour entiers Python (precision arbitraire)
    A = np.array(M, dtype=object)
    R, C = A.shape
    prev = 1
    rank = 0
    pivot_row = 0
    for col in range(C):
        if pivot_row >= R: break
        # Pivot
        nonzero = np.where(A[pivot_row:, col] != 0)[0]
        if len(nonzero) == 0: continue
        piv = pivot_row + int(nonzero[0])
        if piv != pivot_row:
            A[[pivot_row, piv]] = A[[piv, pivot_row]]
        pivot_val = A[pivot_row, col]
        # Mise a jour des autres lignes
        for r in range(R):
            if r == pivot_row: continue
            arc = A[r, col]
            # A[r, col+1:] = (pivot_val * A[r, col+1:] - arc * A[pivot_row, col+1:]) // prev
            A[r, col+1:] = (pivot_val * A[r, col+1:] - arc * A[pivot_row, col+1:]) // prev
            A[r, col] = 0
        prev = pivot_val
        pivot_row += 1
        rank += 1
    return rank


def build_dense_matrix(rows, cols, vals, R, C):
    """Construit matrice dense int Python depuis COO."""
    M = [[0]*C for _ in range(R)]
    for r, c, v in zip(rows.tolist(), cols.tolist(), vals.tolist()):
        M[r][c] += int(v)
    return M


def transpose(M):
    if not M: return []
    R = len(M); C = len(M[0])
    return [[M[r][c] for r in range(R)] for c in range(C)]


def main():
    with open('/home/user/workspace/srmt/step1_detperm.json') as F:
        a = json.load(F)
    det3 = a['det3']; perm3 = a['perm3']
    
    # Charge resultats modulaires
    with open('/home/user/workspace/srmt/results_full.json') as F:
        prev = json.load(F)
    
    # Filtre candidats <= 5*10^6
    CAP = 5_000_000
    targets = [c for c in prev['candidates'] if c['prod'] <= CAP]
    targets.sort(key=lambda c: c['prod'])
    
    results = []
    OUT = '/home/user/workspace/srmt/bareiss_results.json'
    LOG = '/home/user/workspace/srmt/bareiss_log.txt'
    if os.path.exists(LOG): os.remove(LOG)
    
    def log(msg):
        with open(LOG, 'a') as f: f.write(msg + '\n')
        print(msg, flush=True)
    
    # Reprise si checkpoint
    if os.path.exists(OUT):
        with open(OUT) as F:
            existing = json.load(F)
        results = existing['candidates']
        done_keys = {(tuple(r['lam']), tuple(r['mu'])) for r in results}
        log(f"Reprise: {len(done_keys)} deja faits.")
    else:
        done_keys = set()
    
    log(f"=== Bareiss exact sur {len(targets)} candidats (prod <= {CAP:,}) ===")
    
    for ci, c in enumerate(targets):
        if (tuple(c['lam']), tuple(c['mu'])) in done_keys:
            continue
        lam = tuple(c['lam']); mu = tuple(c['mu'])
        prod = c['prod']
        log(f"\n[{ci+1}/{len(targets)}] a={c['a']} lam={lam} mu={mu} prod={prod:,}")
        
        # Construire les matrices via build sparse + dense (plus rapide que build_gamma_matrix)
        try:
            rd = build_gamma_sparse_arrays(lam, mu, det3, n=9)
            rp = build_gamma_sparse_arrays(lam, mu, perm3, n=9)
        except Exception as e:
            log(f"  BUILD ERROR: {e}")
            results.append({**c, 'bareiss_status': f'BUILD_ERR: {e}'})
            continue
        
        R = rd[3]; C = rd[4]
        # OPTIMISATION: on transpose si plus de lignes que de colonnes (Bareiss complexity O(min(R,C)*R*C))
        # En fait Bareiss fait min(R,C) iterations de pivot, mais reduction sur toute la matrice.
        # Si R >> C, on prefere travailler avec la transposee.
        log(f"  shape {R} x {C}; nnz_det={len(rd[2])} nnz_perm={len(rp[2])}")
        
        M_det = build_dense_matrix(*rd)
        M_perm = build_dense_matrix(*rp)
        
        # Pour le rang: rang(M) = rang(M^T). Choisir la dimension la plus petite pour col-iter.
        # Notre Bareiss itere sur col, donc col-iter = C iterations. Si C < R: tel quel; sinon transposer.
        if C > R:
            M_det = transpose(M_det); M_perm = transpose(M_perm)
            log(f"  Transpose: now {len(M_det)} x {len(M_det[0])}")
        
        t0 = time.time()
        try:
            r_det = bareiss_rank(M_det)
            log(f"  rank_Q(det) = {r_det}  (t = {time.time()-t0:.1f}s)")
        except MemoryError:
            log("  OOM on det Bareiss")
            results.append({**c, 'bareiss_status': 'OOM_DET'})
            continue
        
        t1 = time.time()
        try:
            r_perm = bareiss_rank(M_perm)
            log(f"  rank_Q(perm) = {r_perm}  (t = {time.time()-t1:.1f}s)")
        except MemoryError:
            log("  OOM on perm Bareiss")
            results.append({**c, 'bareiss_status': 'OOM_PERM', 'bareiss_rank_det': r_det})
            continue
        
        out = {**c, 'bareiss_rank_det': r_det, 'bareiss_rank_perm': r_perm,
               'bareiss_t_det_s': round(time.time()-t0, 2),
               'matches_modular': r_det == c['rank_det_p1'] and r_perm == c['rank_perm_p1']}
        if r_det == r_perm:
            out['bareiss_status'] = 'EQUAL_Q'
        else:
            out['bareiss_status'] = 'DIFFERENT_Q'
            log(f"  !!!!!!!! DIFFERENCE EXACTE Q DETECTEE !!!!!!!!")
        results.append(out)
        log(f"  matches_modular={out['matches_modular']}")
        
        with open(OUT,'w') as F:
            json.dump({'candidates': results, 'cap': CAP}, F, indent=2, default=str)
    
    log(f"\n=== Termine. {len(results)} traites ===")
    statuses = {}
    for r in results:
        s = r.get('bareiss_status','?'); statuses[s] = statuses.get(s,0)+1
    for s,n in sorted(statuses.items()): log(f"  {s}: {n}")
    
    # Verifier matches_modular sur tous
    all_match = all(r.get('matches_modular', False) for r in results if 'bareiss_rank_det' in r)
    log(f"\n  Tous concordent avec modulaire? {all_match}")


if __name__ == "__main__":
    main()
