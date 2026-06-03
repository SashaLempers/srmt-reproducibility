"""
Etape 7 - Verdict final. Consolidation des resultats et certif SHA-256.
"""
import json, hashlib, os, sys
sys.path.insert(0, '/home/user/workspace/srmt')

# Lire tous les artefacts
with open('/home/user/workspace/srmt/audit.json') as F: audit = json.load(F)
with open('/home/user/workspace/srmt/step1_detperm.json') as F: detperm = json.load(F)
with open('/home/user/workspace/srmt/results_full.json') as F: results = json.load(F)

# Compter par statut
statuses = {}
for c in results['candidates']:
    s = c.get('status','?'); statuses[s] = statuses.get(s,0)+1

# Liste compacte des candidats
compact = []
for c in results['candidates']:
    compact.append({
        'a': c['a'],
        'lam': tuple(c['lam']),
        'mu': tuple(c['mu']),
        'dim_lam_x_dim_mu': c['prod'],
        'rank_det_p65521': c.get('rank_det_p1'),
        'rank_perm_p65521': c.get('rank_perm_p1'),
        'status': c.get('status'),
    })

# Verification: aucun candidat avec rang different
any_diff = any(c['status']=='SUCCESS_CERTIFIED' for c in results['candidates'])
n_equal_modular = statuses.get('EQUAL_MODULAR', 0)
n_total = sum(statuses.values())

print(f"Total candidats: {n_total}")
print(f"  EQUAL_MODULAR (rang det = rang perm sur p=65521): {n_equal_modular}")
for s,n in sorted(statuses.items()):
    if s != 'EQUAL_MODULAR': print(f"  {s}: {n}")
print(f"\nAucun candidat avec rang(Γ(det3)) != rang(Γ(perm3))? {not any_diff}")

# SHA-256 des fichiers
def sha(path):
    with open(path,'rb') as f: return hashlib.sha256(f.read()).hexdigest()

files_to_hash = [
    '/home/user/workspace/srmt/audit.json',
    '/home/user/workspace/srmt/step1_detperm.json',
    '/home/user/workspace/srmt/results_full.json',
    '/home/user/workspace/srmt/step2_schur.py',
    '/home/user/workspace/srmt/step3_pieri.py',
    '/home/user/workspace/srmt/step3b_gamma.py',
    '/home/user/workspace/srmt/step4b_fast_modular.py',
    '/home/user/workspace/srmt/step4c_full.py',
]
hashes = {os.path.basename(p): sha(p) for p in files_to_hash if os.path.exists(p)}

# Sauvegarder
final = {
    'verdict': 'B',
    'audit': audit,
    'det3_sha256': detperm['det3_sha256'],
    'perm3_sha256': detperm['perm3_sha256'],
    'basis_sha256': detperm['basis_sha256'],
    'n_candidates_total': n_total,
    'n_equal_modular': n_equal_modular,
    'modular_prime_used': results.get('p_fast', 65521),
    'statuses': statuses,
    'all_candidates': compact,
    'file_hashes': hashes,
}
with open('/home/user/workspace/srmt/final_verdict.json','w') as F: 
    json.dump(final, F, indent=2, default=str)

# Generation sha256sums.txt
with open('/home/user/workspace/srmt/sha256sums.txt','w') as F:
    for p in files_to_hash + ['/home/user/workspace/srmt/final_verdict.json']:
        if os.path.exists(p):
            F.write(f"{sha(p)}  {os.path.basename(p)}\n")

print("\nFichiers consolides:")
print("  final_verdict.json")
print("  sha256sums.txt")
print("\nSHA-256 des fichiers cles:")
for k,v in hashes.items():
    print(f"  {v}  {k}")
