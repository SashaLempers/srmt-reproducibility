#!/usr/bin/env python3
# =====================================================================
# gen_m2_data.py
#
# Generates code/m2_data.txt consumed by the Python certification core
# (srmt_v19_sparse_search.py / srmt_v19_full.py / srmt_v19_certify.py).
#
# This is an exact, dependency-light replacement for the Macaulay2 dump
# step. It reproduces the SAME conventions as srmt_v15_pipeline.m2:
#
#   * Ring  : QQ[y_0..y_8], 9 variables arranged as a 3x3 matrix Y with
#             matIdx(i,j) = 3*i + j  =>  Y = [[y0,y1,y2],[y3,y4,y5],[y6,y7,y8]]
#   * Basis : the 165 degree-3 monomials mons3, ordered by the nested loop
#                 for a in 0..8: for b in a..8: for c in b..8: y_a y_b y_c
#             (exactly mons3OnRing in the .m2 file).
#   * gl_9  : generator E_{ab} acts on S_1 = Sym^3(C^9) as the derivation
#                 y_a * d/dy_b
#             written in the mons3 basis as  c_k -> sum_i M[i,k] c_i.
#   * det3  : det Y, expanded in the mons3 basis -> integer coefficient
#             vector vDet3 of length 165.
#
# Output format (read by load_m2_data):
#     165
#     <165 lines: 9 exponents each, summing to 3>
#     <gl blocks>:  "a b" / triplets "i k v" / "ENDM"  ... then "ENDGL"
#     <165 integers: vDet3>
# =====================================================================

from itertools import product
import sympy as sp

N = 9
DEG = 3
N3 = 165  # dim Sym^3(C^9) = C(11,3)

# 9 variables y_0..y_8
y = sp.symbols('y_0:9')

# --- 1. mons3 basis in EXACT M2 order (nested a<=b<=c loop) -----------
mons3 = []          # list of sympy monomials
mons3_exp = []      # list of exponent tuples (len 9, sum 3)
for a in range(N):
    for b in range(a, N):
        for c in range(b, N):
            mon = y[a] * y[b] * y[c]
            exp = [0] * N
            exp[a] += 1
            exp[b] += 1
            exp[c] += 1
            mons3.append(sp.expand(mon))
            mons3_exp.append(tuple(exp))
assert len(mons3) == N3, len(mons3)

# index lookup: exponent tuple -> basis index
exp_to_idx = {e: i for i, e in enumerate(mons3_exp)}


def poly_to_vec(poly):
    """Expand a degree-3 polynomial into the mons3 coordinate vector."""
    poly = sp.expand(poly)
    pd = sp.Poly(poly, *y)
    vec = [0] * N3
    for monom, coeff in pd.terms():
        # monom is a tuple of 9 exponents
        idx = exp_to_idx[tuple(monom)]
        vec[idx] = int(coeff)
    return vec


# --- 2. gl_9 action: E_{ab} = y_a * d/dy_b on each basis monomial -----
# c_k -> E_{ab} c_k = sum_i M[i,k] c_i, stored as col[k][i] = M[i,k]
gl_mats = {}
for a in range(N):
    for b in range(N):
        col = {}  # k -> {i: coef}
        for k in range(N3):
            mon = mons3[k]
            deriv = y[a] * sp.diff(mon, y[b])
            if deriv == 0:
                continue
            vec = poly_to_vec(deriv)
            entry = {i: v for i, v in enumerate(vec) if v != 0}
            if entry:
                col[k] = entry
        gl_mats[(a, b)] = col

# --- 3. det3 coefficient vector --------------------------------------
Y = sp.Matrix([[y[3 * i + j] for j in range(3)] for i in range(3)])
det3 = sp.expand(Y.det())
vDet3 = poly_to_vec(det3)
assert sum(1 for v in vDet3 if v != 0) == 6, "det3 must have 6 terms"

# --- 4. write m2_data.txt --------------------------------------------
import os
out_path = os.path.join(os.path.dirname(__file__), "python", "m2_data.txt")
out_path = os.path.abspath(out_path)
# also place a copy next to the python scripts (./m2_data.txt convention)
with open(out_path, "w") as f:
    f.write(f"{N3}\n")
    for e in mons3_exp:
        f.write(" ".join(str(x) for x in e) + "\n")
    for (a, b), col in gl_mats.items():
        if not col:
            continue
        f.write(f"{a} {b}\n")
        for k in sorted(col):
            for i in sorted(col[k]):
                f.write(f"{i} {k} {col[k][i]}\n")
        f.write("ENDM\n")
    f.write("ENDGL\n")
    f.write(" ".join(str(v) for v in vDet3) + "\n")

print(f"Wrote {out_path}")
print(f"  N3={N3}, gl generators with nonzero action: "
      f"{sum(1 for c in gl_mats.values() if c)}")
print(f"  det3 support indices: {[i for i,v in enumerate(vDet3) if v]}")
print(f"  det3 coeffs        : {[v for v in vDet3 if v]}")
