"""
Etape 3 - Construction Gamma_{lambda, mu} via la regle de Pieri.

Pour V = C^n, f in Sym^3(V) avec base monomiale {x_a x_b x_c : a <= b <= c, a,b,c in {0..n-1}}.
Pour lambda de longueur <= n, base de S_lambda(V) = SSYT(lambda, n).

Multiplication Pieri:
  Pour T SSYT de forme lambda et monomial m = x_a x_b x_c (multiset),
  T * m = somme sur tous les SSYT T' de forme mu (mu/lambda = horizontal strip de taille 3)
  obtenus en inserant les valeurs {a,b,c} (multiset, indexes 1..n) en horizontal strip dans T,
  i.e. ajout d'au plus une case par colonne, contenu mass {a+1,b+1,c+1} (decalage 0->1-indexed),
  resultat etant un SSYT (rows weakly increasing, columns strictly increasing).

L'application Gamma_{lambda, mu}(f) : S_lambda V -> S_mu V s'ecrit dans les bases SSYT(lambda,n) -> SSYT(mu,n)
comme: pour chaque T in SSYT(lambda,n), 
  Gamma(f)(T) = somme sur monomes x_a x_b x_c (avec coef f_{abc}) de la projection de T*x_axbxc sur S_mu.

On stocke Gamma(f) comme matrice M de taille |SSYT(mu,n)| x |SSYT(lambda,n)|, entiere.

Convention 1-indexée des SSYT (alphabet {1..n}). Le polynome f est donné par ses 165 coefs sur base
monomiale x_a x_b x_c indexée 0..8. Conversion: 0-index a -> 1-index a+1.

Validation: pour lambda = (), S_lambda V = C (un seul SSYT vide).
  Gamma(f)(empty SSYT) = somme sur monomes de coef * SSYT obtenu en mettant ce monome comme SSYT de forme (3).
  Donc Gamma_{(), (3)}(f) : C -> S_{(3)}V = Sym^3V est juste f lui meme (col unique = vecteur 165 = coefs de f).
  Rank = 1 si f != 0, 0 sinon.
"""
import json, sys, time, hashlib
from itertools import combinations
sys.path.insert(0, '/home/user/workspace/srmt')
from step2_schur import partitions, conjugate, dim_Schur, SSYT, count_SSYT

# Lecture de det3 et perm3
with open('/home/user/workspace/srmt/step1_detperm.json') as F:
    artifact = json.load(F)
det3 = artifact['det3']
perm3 = artifact['perm3']

# Base monomiale de Sym^3(C^9): triples (a,b,c) 0<=a<=b<=c<=8
mono_basis = []
for a in range(9):
    for b in range(a,9):
        for c in range(b,9):
            mono_basis.append((a,b,c))
assert len(mono_basis) == 165
mono_index = {m:k for k,m in enumerate(mono_basis)}

# --- Operation de Pieri ---
# horizontal_strip(lam, k): genere tous les mu obtenus en ajoutant k cases en bande horizontale a lam.
# Une bande horizontale: 0 ou 1 case par colonne du diagramme. 
# mu_i >= lam_i for all i, mu_i <= lam_{i-1} (with lam_0 = infinity for i=0).
# Equivalent: at most one box added per column.

def horizontal_strips(lam, k, max_len):
    """Yields all partitions mu obtained from lam by adding a horizontal strip of size k.
       max_len: limite la longueur de mu (mu's de longueur > max_len rejetes)."""
    lam_list = list(lam) + [0]
    # Number of rows in mu = max(len(lam), ?), and we may add a new row at index len(lam) if lam_{len-1} > 0
    # We extend lam with trailing zeros up to max_len+1
    L = max(len(lam_list), max_len)
    lam_ext = list(lam) + [0]*(L - len(lam))  # length L
    # mu_i: at row i, choose extra = mu_i - lam_i >= 0
    # Constraints: mu_0 >= mu_1 >= ... >= mu_{L-1} >= 0
    #              mu_i <= lam_{i-1} for i >= 1  (horizontal strip: cannot add box in column where row above hasn't a box)
    #              Total extras sum to k
    
    def rec(i, remaining, mu_so_far):
        if remaining == 0:
            # Pad mu_so_far with copies of lam_ext[j] for j>=i (no boxes added in those rows)
            mu = mu_so_far + lam_ext[i:]
            # Trim trailing zeros
            while mu and mu[-1] == 0: mu.pop()
            if len(mu) > max_len: return
            yield tuple(mu)
            return
        if i >= L: return
        lam_i = lam_ext[i]
        # max for mu_i:
        if i == 0:
            up = lam_i + remaining  # no constraint above
        else:
            up = min(mu_so_far[i-1], lam_ext[i-1])  # mu_i <= mu_{i-1} AND mu_i <= lam_{i-1} (horiz strip)
        lo = lam_i
        for mu_i in range(lo, up + 1):
            add = mu_i - lam_i
            if add > remaining: break
            yield from rec(i+1, remaining - add, mu_so_far + [mu_i])
    yield from rec(0, k, [])

# Tests Pieri
def test_horizontal_strip():
    # lam = (), k = 3, max_len = 9: partitions of 3 with parts at most ... no constraint other than horizontal.
    # mu has length 1 only (single row): (3,), since horizontal strip from empty = single row
    res = sorted(set(horizontal_strips((), 3, 9)))
    assert res == [(3,)], res
    # lam = (1,), k = 3: mu = lam + horizontal strip. mu_0 >= mu_1, mu_1 <= lam_0 = 1.
    # Extras: e0+e1 = 3 with mu_0 = 1+e0, mu_1 = e1 <= 1.
    # (e0,e1) in {(3,0),(2,1)}. So mu in {(4,), (3,1)}.
    res = sorted(set(horizontal_strips((1,), 3, 9)))
    assert res == [(3,1),(4,)], res
    # lam = (1,1), k = 3: mu_0 = 1+e0, mu_1 = 1+e1 <= mu_0 and <= lam_0 = 1, so e1 = 0; mu_2 = e2 <= lam_1 = 1.
    # e0 + e2 = 3, e2 <= 1, so (3,0,0) -> (4,1), (2,0,1) -> (3,1,1).
    res = sorted(set(horizontal_strips((1,1), 3, 9)))
    assert res == [(3,1,1),(4,1)], res
    # lam = (2,1), k = 3: classic. mu obtained adding 3 boxes in horizontal strip to (2,1).
    # Reference: (2,1) * h_3 = s_{5,1} + s_{4,2} + s_{4,1,1} + s_{3,2,1} + s_{3,3} + s_{2,2,2}? Let's see.
    # Actually s_{(2,1)} * h_3 = sum of s_mu over mu/lam horizontal strip of size 3.
    # mu must contain (2,1). Possibilities found by enumeration.
    res = sorted(set(horizontal_strips((2,1), 3, 9)))
    expected = sorted([(5,1),(4,2),(4,1,1),(3,3),(3,2,1),(3,2,1),(2,2,2)])
    # Need to recompute manually: mu_0 = 2+e0 >= mu_1 = 1+e1 <= 2 (so e1 <= 1); mu_2 = e2 <= 1.
    # e0+e1+e2=3, e1 in {0,1}, e2 in {0,1}.
    # (3,0,0)->(5,1); (2,1,0)->(4,2); (2,0,1)->(4,1,1); (1,1,1)->(3,2,1); ... wait need e1 <= 1. 
    # Also (1,0, ?) with e2<=1: (1,0,1)? sum=2, not 3. e0=1 means sum left = 2 distributed in e1,e2 each <=1: e1=e2=1, total (1,1,1)->(3,2,1).
    # (0,1,?): sum left=2 with e1=1, e2 <=1, e2=2? no.
    # We also have mu_0 >= mu_1: with e0=0, mu_0=2; mu_1=1+e1; if e1=1, mu_1=2 OK.
    # (0,1,?): need 2 more in e2... but e2<=1. Impossible.
    # So muset: (5,1),(4,2),(4,1,1),(3,2,1) only? Hmm let me re-examine.
    # Wait I also need to consider e0 + e1 + e2 partitioning, and we can ALSO add box in column 3 etc. The "horizontal strip" allows at most one new box per column.
    # Let me draw (2,1):
    #   X X
    #   X
    # Columns: col 0 has 2 boxes, col 1 has 1 box, col 2 has 0 boxes, etc.
    # Adding horizontal strip of size 3: at most one box per column. Boxes added must be in valid positions (each in a row right-adjacent to existing or extending).
    # Possible additions: cells (i,j) with j >= lam_i (extending row i to the right). New cells: at row 0, j>=2; row 1, j>=1; row 2, j>=0.
    # Horizontal strip: each column has at most one new cell. So if we add (0,j0), (1,j1), (2,j2), they must be in distinct columns.
    # (0,j0): j0>=2 column j0. (1,j1): j1>=1 column j1. (2,j2): j2>=0 column j2.
    # Choose how many cells per row: r0 cells in row 0 (consecutive starting from col 2), r1 in row 1 (starting col 1), r2 in row 2 (starting col 0).
    # r0+r1+r2=3. Columns occupied by row 0: {2,3,...,2+r0-1}; by row 1: {1,2,...,r1}; by row 2: {0,...,r2-1}.
    # Disjoint columns: row 0 starts at col 2; row 1 starts at col 1. So row 0 col 2..2+r0-1 and row 1 col 1..r1. Disjoint iff r1 <= 1 (else row 1 reaches col 2, conflict with row 0 if r0>=1).
    # And row 1 vs row 2: row 1 cols 1..r1 vs row 2 cols 0..r2-1. Disjoint iff r2 <= 1.
    # Cases:
    # r0=3,r1=0,r2=0: mu=(5,1). cols used: 2,3,4.
    # r0=2,r1=1,r2=0: mu=(4,2). cols: 2,3 and 1. OK.
    # r0=2,r1=0,r2=1: mu=(4,1,1). cols: 2,3 and 0.
    # r0=1,r1=1,r2=1: mu=(3,2,1). cols: 2 and 1 and 0. OK (r1<=1, r2<=1).
    # r0=1,r1=2,r2=0: r1=2 means row 1 occupies cols 1,2, conflict with row 0 col 2. Invalid.
    # r0=0,r1=3,r2=0: row 1 occupies cols 1,2,3. r1=3 > 1 if r0 >=1... but r0=0. Conflict with row 0? No row 0 is empty. Row 2 r2=0. So OK actually! But wait my constraint was r1<=1 only if r0>=1. Re-check: row 1 cols 1..r1, row 0 cols 2..2+r0-1. If r0=0, row 0 is empty, no conflict. r1 can be 3. mu=(2,4)? But then mu_0 < mu_1. That violates partition. So invalid as partition.
    # r0=0,r1=2,r2=1: mu=(2,3,1). Not a partition. Invalid.
    # r0=0,r1=1,r2=2: r2=2 means row 2 cols 0,1, conflict with row 1 col 1. Invalid (unless r1=0).
    # r0=0,r1=0,r2=3: row 2 cols 0,1,2. mu=(2,1,3). Not partition.
    # r0=2,r1=1,r2=0: done.
    # So valid + partition: (5,1), (4,2), (4,1,1), (3,2,1). 4 candidates.
    # My code outputs:
    print(res)
    assert res == [(3,2,1),(4,1,1),(4,2),(5,1)], res
    print("horizontal_strips OK")

test_horizontal_strip()
print()
print("=== Demo Pieri pour lambda = (2,1) ⊢ 3, k=3 ===")
for mu in horizontal_strips((2,1), 3, 9):
    print(f"  mu = {mu}, dim S_mu(C^9) = {dim_Schur(mu,9)}")
