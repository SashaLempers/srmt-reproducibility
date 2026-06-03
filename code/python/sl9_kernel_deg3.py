#!/usr/bin/env python3
# =====================================================================
# sl9_kernel_deg3.py
#
# Reproduces Computation 9 of the AMSART paper:
#   dim Sym^3( Sym^3 C^9 )^{SL_9} = 0
#
# Method (exact over Q):
#   V = C^9, W = Sym^3 V  (basis = 165 cubic monomials x_{i<=j<=k}, 1<=i,j,k<=9)
#   The maximal torus T of SL(V) acts diagonally on Sym^3(W).
#   Restrict to the T-fixed subspace (weight (1,1,...,1) over the 9 variables,
#   i.e. degree-3 multisets of cubic monomials whose total exponent vector
#   on x_1..x_9 is (1,1,...,1)).  This T-fixed subspace has dim 280.
#   The 8 simple raising operators e_{a,a+1} (a=1..8) of sl_9 act on Sym^3(W)
#   via the Leibniz rule (derivation extending the action on W).
#   SL_9 invariants = common kernel of all 8 e_{a,a+1} restricted to (Sym^3 W)^T.
#   Stack the 8 images -> an 11200 x 280 rational matrix; its kernel dimension
#   is dim of the invariant space.
#
# Output: prints dims, writes a transcript + ranks.
# =====================================================================

import sys, time, json
from itertools import combinations_with_replacement as cwr
from collections import defaultdict
from fractions import Fraction

import sympy
from sympy import Rational, Matrix, zeros

N = 9  # ambient C^9

# ---------- Basis of W = Sym^3 V : cubic monomials x_{i<=j<=k} -------------
# Represent a cubic monomial by its sorted index triple (i,j,k), 1-based.
cubic_monos = list(cwr(range(1, N + 1), 3))   # 165 of them
assert len(cubic_monos) == 165, len(cubic_monos)
mono_index = {m: t for t, m in enumerate(cubic_monos)}

def exp_vector(triple):
    """Exponent vector on x_1..x_9 for cubic monomial given by index triple."""
    v = [0] * N
    for i in triple:
        v[i - 1] += 1
    return tuple(v)

cubic_exp = [exp_vector(m) for m in cubic_monos]  # length-9 exponent vectors

# ---------- T-fixed subspace of Sym^3(W) ----------------------------------
# An element of Sym^3(W) is a degree-3 multiset of cubic monomials:
#   { m_a, m_b, m_c } with a<=b<=c indices into cubic_monos.
# Its T-weight = sum of the three exponent vectors.
# T-fixed (weight (1,1,...,1)) means the three cubic monomials together use
# each variable x_1..x_9 exactly once  (total degree 9 = 3*3, balanced).
TARGET = (1,) * N

basis_Sym3W = []   # list of sorted triples (a,b,c) of cubic-mono indices
for a, b, c in cwr(range(165), 3):
    e = tuple(cubic_exp[a][t] + cubic_exp[b][t] + cubic_exp[c][t] for t in range(N))
    if e == TARGET:
        basis_Sym3W.append((a, b, c))

dimT = len(basis_Sym3W)
print(f"dim (Sym^3 W)^T = {dimT}  (expected 280)", flush=True)
sym3w_index = {b: t for t, b in enumerate(basis_Sym3W)}

# ---------- Action of e_{a,a+1} on W = Sym^3 V ----------------------------
# e_{a,a+1} = x_a * d/dx_{a+1}  (raising operator: replaces one x_{a+1} by x_a).
# On a cubic monomial (multiset of 3 indices), for each occurrence of (a+1),
# replace it by a; coefficient = (number of (a+1) factors).
def raise_on_cubic(triple, a):
    """Return dict {result_cubic_index: integer coeff} for e_{a,a+1} applied
    to the cubic monomial 'triple' (sorted index triple)."""
    out = defaultdict(int)
    src = a + 1  # the variable index that gets lowered (a is 1-based here)
    cnt = triple.count(src)
    if cnt == 0:
        return out
    # replace ONE src by a:  derivative brings factor = multiplicity of src
    lst = list(triple)
    pos = lst.index(src)
    lst[pos] = a
    res = tuple(sorted(lst))
    out[mono_index[res]] += cnt
    return out

# Precompute raise action on each cubic-mono index, for each a in 1..8
raise_cubic = [None] * 8  # raise_cubic[a-1][cubic_idx] = dict
for a in range(1, N):  # a = 1..8
    table = {}
    for ci, m in enumerate(cubic_monos):
        table[ci] = dict(raise_on_cubic(m, a))
    raise_cubic[a - 1] = table

# ---------- Action of e_{a,a+1} on Sym^3(W) by Leibniz --------------------
# Element = symmetric product of three cubic monomials (multiset {p,q,r}).
# e acts as derivation: e(p*q*r) = e(p)*q*r + p*e(q)*r + p*q*e(r).
def raise_on_sym3w(pqr, a):
    """pqr = sorted triple of cubic-mono indices. Returns dict
    {sorted-triple : int coeff} = e_{a,a+1} applied (in Sym^3 W)."""
    out = defaultdict(int)
    table = raise_cubic[a - 1]
    p, q, r = pqr
    for slot, others in ((p, (q, r)), (q, (p, r)), (r, (p, q))):
        for new_idx, coeff in table[slot].items():
            res = tuple(sorted((new_idx,) + others))
            out[res] += coeff
    return {k: v for k, v in out.items() if v != 0}

# ---------- Build the 11200 x 280 constraint matrix -----------------------
# Columns indexed by the 280 T-fixed basis elements.
# Rows: for each a=1..8, the image lives in a space spanned by target sorted
# triples; collect the union of all such target triples => row index set.
t0 = time.time()
# Convention of the paper: target space = Sym^3(W) (x) n^+ , i.e. for EACH of
# the 8 simple root vectors e_{a,a+1} we record a separate block of rows, one
# per possible target triple. First collect all target triples that can occur
# (across all columns and all a), then index rows by (a, target_triple): this
# gives 8 * 1400 = 11200 formal rows, matching the paper exactly.
all_targets = set()
col_images = []  # per column: list of (a, tgt, coeff)
for col, pqr in enumerate(basis_Sym3W):
    lst = []
    for a in range(1, N):  # a=1..8
        img = raise_on_sym3w(pqr, a)
        for tgt, coeff in img.items():
            all_targets.add(tgt)
            lst.append((a, tgt, coeff))
    col_images.append(lst)

target_list = sorted(all_targets)
ntargets = len(target_list)
target_idx = {t: i for i, t in enumerate(target_list)}
print(f"distinct target triples = {ntargets}  (paper: 1400 ; 8 x 1400 = 11200)", flush=True)

# Row index for (a, tgt): block a occupies rows [a*ntargets : (a+1)*ntargets]
row_keys = {}
for a in range(8):
    for ti in range(ntargets):
        row_keys[(a + 1, target_list[ti])] = a * ntargets + ti

entries = []
for col, lst in enumerate(col_images):
    for (a, tgt, coeff) in lst:
        entries.append((row_keys[(a, tgt)], col, coeff))

nrows = 8 * ntargets
print(f"matrix: {nrows} x {dimT}  (paper states 11200 x 280)", flush=True)
print(f"nonzero entries: {len(entries)}", flush=True)

# ---------- Exact rank over Q (=> kernel dim = 280 - rank) ----------------
M = zeros(nrows, dimT)
for (i, j, c) in entries:
    M[i, j] = c
print("computing exact rank over Q (Bareiss/fraction-free)...", flush=True)
rank = M.rank()
kerdim = dimT - rank
dt = time.time() - t0
print(f"rank_Q = {rank}", flush=True)
print(f"dim ker = {dimT} - {rank} = {kerdim}", flush=True)
print(f"=> dim Sym^3(Sym^3 C^9)^SL_9 = {kerdim}", flush=True)
print(f"[{dt:.1f}s]", flush=True)

result = {
    "computation": "sl9_kernel_deg3_QQ",
    "claim": "dim Sym^3(Sym^3 C^9)^{SL_9} = 0",
    "dim_T_fixed": dimT,
    "matrix_rows": nrows,
    "matrix_cols": dimT,
    "rank_Q": int(rank),
    "dim_ker": int(kerdim),
    "verdict": "PASS" if kerdim == 0 else "FAIL",
    "elapsed_s": round(dt, 2),
}
with open("/home/user/workspace/sl9_kernel_deg3_result.json", "w") as f:
    json.dump(result, f, indent=2)
print("written sl9_kernel_deg3_result.json", flush=True)
