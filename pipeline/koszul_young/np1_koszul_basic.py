"""
NP-1: Koszul flattening YF_{3,9}(f) basique (cas Pieri-hook de LO11 §1.2 Remarque 4.1.3).

Pour f in Sym^3(C^9), delta = 1, a in {1,2,3,4}:
   YF_a(f) : V^* (x) Lambda^a V  -->  V (x) Lambda^{a+1} V
   
Matrice block-form : block (I,J) ou |I|=a, |J|=a+1 :
   bloc = +/- (df/dx_k)  si J = I cup {k} (signe = (-1)^{position de k dans J trie})
   bloc = 0              sinon.

Donc la matrice est ((9) x C(9,a))  par ((9) x C(9,a+1)),
   forme: (9 * C(9,a))  x  (9 * C(9,a+1))    ... mais c'est V^* tensor Lambda^a V -> V tensor Lambda^{a+1} V
   donc rows index = (V tensor Lambda^{a+1}V) : 9 * C(9,a+1)
        cols index = (V* tensor Lambda^a V)  : 9 * C(9,a)

Notation index : (i, J) pour cols (i in V*, J in Lambda^a),
                 (j, K) pour rows (j in V,  K in Lambda^{a+1}).

Action (formule LO) : alpha (x) omega -> sum_k (df/dx_k * alpha(e_k)) (x) (e_k wedge omega)
Donc M[(j, K), (i, J)] = somme sur k de  
   (coef de monome x_j dans df/dx_k) * delta(i, k) * coef Koszul (K = e_k wedge J).
   
Puisque alpha = e_i^*, alpha(e_k) = delta_{ik}, donc k = i fixe :
   M[(j, K), (i, J)] = (df/dx_i)_j * sign(K -> i wedge J)
ou sign = +1 si K = {i} cup J avec i a la position canonique et J trie, ajuste pour permutation;
   0 si i in J ou K != {i} cup J.

Pour d=3, df/dx_i est un polynome de degre 2 en 9 vars : donc (df/dx_i)_j (coef de x_j dans df/dx_i)
n'est PAS un scalaire, c'est un polynome lineaire en les 8 autres vars. ATTENTION.

Reformulons : le bloc (J,K) est une matrice 9x9 a coef polynomiaux ? Non, LO dit
"square catalecticants (df/dx_i)_{delta,delta}" :
pour delta=1 ce sont les coefs lineaires de df/dx_i, qui sont les coefs quadratiques de f.
Soit M_i = la matrice 9x9 Hessienne (f_{ij}) i.e. coef de x_i x_j dans f (symetrique).
Alors le bloc K, J est : M_i avec i tel que K = {i} cup J, multiplie par signe.

OK donc plus precisement : delta=1, S^1 V* (x) Lambda^a V -> S^1 V (x) Lambda^{a+1} V :
   (e_i^* (x) e_J)  ->  sum_k  e_k (x) (e_k wedge e_J)  *  H(f)[k,i]
ou H(f)[k,i] est le coef de x_k * x_i (avec multiplicite 2 si k != i, 1 si k=i)... 
attention au coef 1/2 du polynome quadratique df/dx_i.

Pour f = sum_M coef[M] * x^M (M multi-indice de deg 3), 
   df/dx_i = sum_M coef[M] * M[i] * x^{M - e_i}
puis (df/dx_i)[x_k] (coef de x_k dans le polynome quadratique df/dx_i) = 
   = somme sur M tq M >= e_i+e_k de coef[M] * M[i] * (M-e_i)[k]_{coef monome canonique}
Trop alambique. Faisons-le numeriquement.
"""
import json
import sys
import os
import itertools
from pathlib import Path
from fractions import Fraction

sys.path.insert(0, '/home/user/workspace/srmt')
from step2_schur import partitions, dim_Schur

# Charger det3 et perm3 en base multiset monomiale
SRMT = Path('/home/user/workspace/srmt')
with open(SRMT / 'step1_detperm.json') as f:
    s1 = json.load(f)

# La base utilise les variables x_{ij} indexees de 0 a 8 (3x3 matrice flatten ligne-majeur).
# step1 utilise la convention : variable index = 3*i + j pour ligne i, col j (0-indexe).
# Base : multisets de taille 3 dans [0,9), donc C(9+3-1, 3) = 165 monomes.

n = 9
d = 3

# Reconstruire la base monomiale comme dans step1.
mono_basis = []
for i in range(n):
    for j in range(i, n):
        for k in range(j, n):
            mono_basis.append((i, j, k))
assert len(mono_basis) == 165
mono_index = {m: idx for idx, m in enumerate(mono_basis)}

det3_vec = s1['det3']
perm3_vec = s1['perm3']
assert len(det3_vec) == 165 and len(perm3_vec) == 165

# Hessienne : pour un polynome f de degre 3, df/dx_i est de degre 2.
# Coef de x_i*x_j dans f (i<=j) : note H[i][j].
# Si i==j : coef du monome (i,i,?), il y a 3 monomes possibles (i,i,k) k=0..8 avec k != i (et k=i donne (i,i,i)).
# En fait coef de x_i x_j x_k dans f = coef du multiset (i,j,k) trie, par definition de la base.
#
# Pour la Hessienne : H[i][j] = coef de x_i x_j dans (df/dx_? ) ... NON.
# Reformulons : pour la Koszul, on a besoin de la matrice (delta=1, delta=1) catalecticant :
#    Cat(f) : V* -> Sym^2 V* -> ... non, en fait c'est juste la deuxieme derivee :
#    (d^2 f / dx_i dx_j) = forme bilineaire symetrique sur V dont la matrice 9x9 a pour (k,l)-coef
#    le coefficient de x_k x_l dans d^2 f / dx_i dx_j.
# Pour delta=1 : Cat_{1,1}(f) : V* -> V* defini par : alpha -> (d/d alpha)(f) restreint a degre 1 ?
# Non. Reprenons LO11 directement.
#
# YF_{3,9}(f) : delta=1, a=4 typique :
#   S^1 V* (x) Lambda^4 V  ->  S^1 V (x) Lambda^5 V
#
# La construction explicite : 
#   YF(f)(alpha (x) v_1 ^ ... ^ v_a) = sum_k (alpha . partial_k) f * e_k ^ v_1 ^ ... ^ v_a ?
#
# Plus simple : on definit la *catalecticant* Cat_1(f) : V* -> Sym^2 V comme la differentielle
# (apres choix de bases avec forme volume).
# Pour f = sum c_{ijk} x_i x_j x_k (i<=j<=k), 
#   df/dx_p = sum c_{ijk} (mult p in (ijk)) * x_{(ijk) \ p_one}
# c'est un polynome de degre 2 en 9 vars, vu comme element de Sym^2 V dans la base monomiale.
#
# La matrice Koszul ("Koszul-Young flattening") explicite : pour p,q in {0,...,n-1}, 
# bloc B_{pq} = Cat_{1,1}(f)[p, q] = "coef de x_p x_q dans la 'derivee mixte' de f".
# C'est exactement la *Hessienne symetrique* H[p,q] qui est une matrice 9x9 :
#   H[p][q] = coef de x_p x_q dans f, MAIS attention aux multiplicites :
#     si p=q : c'est le coef de x_p^2 dans df/dx_p divise par 2 ? Non, la 2eme derivee de x_p^3 est 6 x_p,
#     coef de x_p dans 6x_p = 6, et le coef "1/6 H" non... rebookons proprement.
#
# Disons : la 2-jet de f en 0 dans la base monomiale (e_p x e_q) selon la convention symetrique :
#   J^2(f)[p,q] = coef du monome (avec ordre) p,q dans f * (mult).
# Ou plus simple : posons H[p,q] := coefficient de x_p x_q dans le polynome 1/2 * d^2 f / dx_p dx_q.
#   d^2(x_a x_b x_c)/dx_p dx_q : non nul ssi {p,q} subset {a,b,c} multiplicativement.
#
# Trop sec. Calculons-le numeriquement par enumeration.

def hess_matrix(f_vec):
    """Renvoie H[p][q] = coefficient (avec convention catalecticant) tel que
       le polynome 'pieces lineaires de df/dx_p en variable q' soit code.
       
       Convention : pour f = sum_{multiset (a,b,c)} f_vec[(a,b,c)] * x_a x_b x_c (forme symetrique),
       df/dx_p = sum_{multiset trie (a,b,c) contenant p} f_vec[(a,b,c)] * mult_p((a,b,c)) * monome (a,b,c) sans un p.
       
       Le 'catalecticant' Cat_{1,1}(f) est la matrice 9x9 H[p][q] = 
         coef de x_q dans df/dx_p, MAIS en regroupant les monomes de degre 2.
       Non, Cat_{1,1} c'est V* -> V tel que :
         f -> map alpha -> "alpha applique sur f" comme element de Sym^2 V.
       Alors Cat_{1,1}(f)[p][q] = coef de x_p x_q dans f * (mult de (p,q) comme monome non-trie).
       
       Convention propre : H[p][q] = coef du monome trie (p,q,?,?) ...
       
       SIMPLIFIONS : H_2(f) = matrice symetrique 9x9 dont (p,q) entry est 
       le coef du monome x_p x_q (vu comme polynome quadratique) dans la derivee partielle df/dx_? 
       pour ? variable. Non c'est tjr ambigu.
       
       Approche pragmatique : la Hessienne Hess(f) au point x est la matrice 9x9 dont (p,q) entry est
       d^2 f / dx_p dx_q (un polynome lineaire en x). Au point x = 0 c'est nul pour f cubique.
       
       En fait pour LO11 §1.2, ce qu'on prend pour delta=1 est juste la matrice 9x9 a coefficients
       polynomiaux d/dx_p evalue comme un polynome lineaire. Ca devient :
         M_a(f)[(j,K),(i,J)] = coef de x_j dans (df/dx_i) * sign(K vs i,J)
                              = mult avec laquelle le monome (i,j,?) apparait
       """
    # Pour chaque triplet (i, j, third) multiset trie, on a un monome m = (a,b,c) avec a<=b<=c.
    # df/dx_i (i fixe) : ce polynome de degre 2 a comme coef de x_j (en convention symetrique)
    # le coefficient du monome (i,j,?) dans f multiplie par la multiplicite avec laquelle (i,j) 
    # apparait dans (i,j,?), divisise par la mult de (?)... 
    #
    # Plus concretement : on veut, pour chaque (i,j) ordered (avec i possiblement = j), une valeur
    # M[i][j] = somme sur k=0..8 de f_vec[sort((i,j,k))] * comb_factor(i,j,k),
    # ou comb_factor est tel que f reconstruit en x_i x_j x_k vaille le bon coef.
    #
    # Soit f(x) = sum_{a<=b<=c} f_vec[(a,b,c)] * x_a x_b x_c.
    # df/dx_i = sum_{a<=b<=c, i in {a,b,c}} f_vec[(a,b,c)] * 
    #             [ (a==i)*x_b*x_c + (b==i)*x_a*x_c + (c==i)*x_a*x_b ]
    # = sum_{(b,c) such that there exists multiset containing i,b,c} ... regroupe par x_b x_c.
    #
    # OK plus simple : df/dx_i en forme totalement symetrique (mais c'est deja un polynome).
    # On l'evalue : pour chaque triplet (i,b,c) avec b<=c (8x9/2 = 45 pairs),
    #   coef de x_b x_c dans df/dx_i = f_vec[sort((i,b,c))] * (mult de i dans (i,b,c)) * (correction sym).
    # Concretement, posons G[i] le polynome de degre 2 en 9 vars : 
    #   G[i] = df/dx_i = sum_{b<=c} c[i,b,c] x_b x_c, c[i,b,c] = ?
    # Pour le monome trie (a,b,c) avec multiplicite, df/dx_i contribue (a==i)*x_b*x_c + (b==i)*x_a*x_c + (c==i)*x_a*x_b.
    # Donc chaque terme f_vec[(a,b,c)] * x_a * x_b * x_c (sans coef multinomial supplementaire ; on a 
    # juste sum coef * x_a x_b x_c) derive en :
    #   (a==i) f_vec[(a,b,c)] x_b x_c + (b==i) f_vec[(a,b,c)] x_a x_c + (c==i) f_vec[(a,b,c)] x_a x_b
    # 
    # Pour collecter le coef de x_b' x_c' (b' <= c') dans df/dx_i : on cherche les triplets (a,b,c) trie tels que
    # ({a,b,c} \ {un i}) (en multiset) = {b', c'}.
    # C'est-a-dire : (a,b,c) trie = sort(i, b', c'). Et chaque match contribue (nb de fois ou i apparait dans
    # (a,b,c) au bon emplacement) * f_vec[(a,b,c)].
    # Donc coef de x_{b'} x_{c'} dans df/dx_i (avec b'<=c') = (nb de fois i apparait dans sort(i,b',c')) * f_vec[sort(i,b',c')].
    # MAIS attention au coef de x_b' x_c' quand b' != c' : le polynome 'x_b' x_c'' est juste un monome,
    # tandis que quand b' = c' on a x_{b'}^2 = monome avec poids 1 dans la base monomiale standard.
    # Donc pour b' != c' : coef = mult(i, sort(i,b',c')) * f_vec[sort(i,b',c')]. Or quand b'!=c',
    # mult(i, sort(i,b',c')) = (i==b')+(i==c')+(... non, i n'est pas dans (b',c')) = 1 (juste le i lui-meme).
    # Hmm. Reprenons : x_b' x_c' (b'!=c') comme monome n'a qu'une seule forme normale ; mult(i,sort(i,b',c'))
    # designe combien de positions de la liste sort(i,b',c') contiennent i. Si i not in {b',c'}, c'est 1.
    # Si i = b' (et b' != c') : c'est 2 (les deux b' = i dans la liste). Etc.
    return None  # placeholder

# OK abandonnons cette derivation a la main et faisons-le programmaitiquement via sympy ou enumeration.

# Construire les polynomes det3, perm3 comme dict {(a,b,c) trie : Fraction(coef)}
def vec_to_poly_dict(vec):
    return {mono_basis[i]: Fraction(v) for i, v in enumerate(vec) if v != 0}

det3_poly = vec_to_poly_dict(det3_vec)
perm3_poly = vec_to_poly_dict(perm3_vec)

print(f"det3 nb monomes non nuls : {len(det3_poly)}")
print(f"perm3 nb monomes non nuls : {len(perm3_poly)}")
print(f"det3 echantillon : {dict(list(det3_poly.items())[:5])}")
print(f"perm3 echantillon : {dict(list(perm3_poly.items())[:5])}")

# Derivee partielle df/dx_i comme dict {(b,c) trie : Fraction(coef)}
# Convention : f = sum coef[(a,b,c)] x_a x_b x_c ou coef est la valeur de f_vec et (a,b,c) trie.
# df/dx_i : on parcourt tous les monomes (a,b,c) contenant i au moins une fois.
# Pour chaque tel monome, df/dx_i de "coef * x_a x_b x_c" = coef * (mult de i dans (a,b,c)) * (monome reste sans un i).

def partial(poly_dict, i):
    """df / dx_i comme dict {(b,c) trie : Fraction}."""
    out = {}
    for (a, b, c), v in poly_dict.items():
        # mult de i dans (a,b,c)
        mult = (a == i) + (b == i) + (c == i)
        if mult == 0:
            continue
        # retire un i et obtient une paire triee
        lst = [a, b, c]
        lst.remove(i)
        key = tuple(sorted(lst))
        out[key] = out.get(key, Fraction(0)) + v * mult
    # nettoie zeros
    return {k: v for k, v in out.items() if v != 0}

# Test : df/dx_0 de det3, devrait etre x_4 x_8 - x_5 x_7 (mineur 2x2 inferieur droit).
# Notre convention : x_{i,j} = variable index 3*i+j. Donc x_{1,1} = idx 4, x_{2,2} = idx 8, x_{1,2}=idx 5, x_{2,1}=idx 7.
# det3 = x_00 x_11 x_22 - x_00 x_12 x_21 - x_01 x_10 x_22 + x_01 x_12 x_20 + x_02 x_10 x_21 - x_02 x_11 x_20
# d/dx_00 det3 = x_11 x_22 - x_12 x_21 = x_4 x_8 - x_5 x_7
df_dx0 = partial(det3_poly, 0)
print(f"d det3/dx_0 = {df_dx0}")
expected = {(4, 8): Fraction(1), (5, 7): Fraction(-1)}
assert df_dx0 == expected, f"MISMATCH: got {df_dx0}, expected {expected}"
print("OK partial(det3, x_0) correcte.")

# d perm3 / d x_0 = x_4 x_8 + x_5 x_7
df_p0 = partial(perm3_poly, 0)
expected_p = {(4, 8): Fraction(1), (5, 7): Fraction(1)}
assert df_p0 == expected_p, f"MISMATCH perm: got {df_p0}"
print("OK partial(perm3, x_0) correcte.")

print("Toutes les derivees partielles d/dx_i pour det3 et perm3 :")
det3_partials = [partial(det3_poly, i) for i in range(n)]
perm3_partials = [partial(perm3_poly, i) for i in range(n)]
print(f"  det3 partials : {[len(p) for p in det3_partials]} monomes (par variable)")
print(f"  perm3 partials : {[len(p) for p in perm3_partials]} monomes")

# Sauvegarder
import pickle
with open(SRMT / 'np1_partials.pkl', 'wb') as f:
    pickle.dump({
        'det3_partials': [{k: int(v) for k, v in p.items()} for p in det3_partials],
        'perm3_partials': [{k: int(v) for k, v in p.items()} for p in perm3_partials],
    }, f)
print("Sauvegarde dans np1_partials.pkl OK")
