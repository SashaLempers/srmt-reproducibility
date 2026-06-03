#!/usr/bin/env python3
# =====================================================================
# hessian_perm3.py
#
# Exact rational (SymPy) certification of the degree-9 Hessian-determinant
# claims for perm_3 and det_3 in 9 variables x_11..x_33, arranged as
#
#        [ x11 x12 x13 ]
#    X = [ x21 x22 x23 ]
#        [ x31 x32 x33 ]
#
# This script is the reproducible computation behind, in the manuscript:
#
#   * Theorem  (Hessian determinant of det_3):
#        det Hess(det_3) = -2 (det_3)^3.
#   * Theorem  (Hessian determinant of perm_3), parts (1)-(4):
#        D := det Hess(perm_3) is NOT proportional to (perm_3)^3;
#        residual R = D + 2(perm_3)^3 has exactly 37 monomials;
#        distinct coefficient ratios D/(perm_3)^3 on common support
#        are {-2, 2/3, 10/3}; D is irreducible over Q (single factor,
#        degree 9, multiplicity 1, 55 monomials, unit -2); D is not a
#        rational cube.
#   * Lemma  (D is squarefree over C) -- THE HEADLINE LEMMA:
#        g := gcd(d D/d x11, ..., d D/d x33) is a nonzero constant,
#        hence D is squarefree over C and in particular NOT a cube
#        over C. Field-independent in characteristic 0.
#   * Remark (independent line test):
#        the restriction of D to a generic rational line has 9 distinct
#        roots (gcd of the univariate restriction with its derivative
#        is constant); a cube would have all roots of multiplicity 3.
#
# All arithmetic is exact over Q. No floating point. random.seed(42).
#
# Usage:
#     python3 hessian_perm3.py
# Exit code 0 iff every claim is verified.
# =====================================================================

import sys
import random
import sympy as sp

random.seed(42)

# --- variables, arranged as a 3x3 matrix --------------------------------
x = sp.symbols('x11 x12 x13 x21 x22 x23 x31 x32 x33')
X = sp.Matrix([[x[0], x[1], x[2]],
               [x[3], x[4], x[5]],
               [x[6], x[7], x[8]]])
VARS = list(x)


def perm3(M):
    """Permanent of a 3x3 matrix."""
    import itertools
    s = 0
    for sigma in itertools.permutations(range(3)):
        prod = 1
        for i in range(3):
            prod *= M[i, sigma[i]]
        s += prod
    return sp.expand(s)


def hessian_det(f):
    """det of the 9x9 Hessian of f w.r.t. VARS, expanded over Q."""
    H = sp.hessian(f, VARS)        # 9x9 symmetric matrix of 2nd partials
    return sp.expand(H.det())


def nterms(poly):
    """Number of monomials of an expanded polynomial."""
    poly = sp.expand(poly)
    if poly == 0:
        return 0
    return len(sp.Poly(poly, *VARS).terms())


results = []  # (name, ok, detail)


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    flag = "PASS" if ok else "FAIL"
    print(f"[{flag}] {name}" + (f"  --  {detail}" if detail else ""))


# =====================================================================
# det_3 baseline:  det Hess(det_3) = -2 (det_3)^3
# =====================================================================
print("=== det_3 : classical Hessian identity ===")
det3 = sp.expand(X.det())
D_det = hessian_det(det3)
target = sp.expand(-2 * det3**3)
check("det Hess(det3) = -2 (det3)^3", sp.expand(D_det - target) == 0)

# =====================================================================
# perm_3 : Theorem (Hessian determinant of perm_3)
# =====================================================================
print("\n=== perm_3 : Hessian determinant D = det Hess(perm3) ===")
perm = perm3(X)
D = hessian_det(perm)                       # degree 9, 9 variables
check("D = det Hess(perm3) is nonzero, homogeneous degree 9",
      D != 0 and sp.Poly(D, *VARS).is_homogeneous and sp.total_degree(D) == 9,
      f"deg = {sp.total_degree(D)}")

# (1) residual R = D + 2 (perm3)^3 has exactly 37 monomials
R = sp.expand(D + 2 * perm**3)
nR = nterms(R)
check("(1) residual R = D + 2 perm3^3 has exactly 37 monomials",
      R != 0 and nR == 37, f"|R| = {nR}")

# (2) distinct coefficient ratios D / perm3^3 on common support = {-2, 2/3, 10/3}
P3 = sp.expand(perm**3)
pd_D = sp.Poly(D, *VARS).as_dict()
pd_P = sp.Poly(P3, *VARS).as_dict()
common = set(pd_D) & set(pd_P)
ratios = {sp.nsimplify(sp.Rational(pd_D[m], pd_P[m])) for m in common}
expected_ratios = {sp.Rational(-2), sp.Rational(2, 3), sp.Rational(10, 3)}
check("(2) distinct ratios D/perm3^3 on common support = {-2, 2/3, 10/3}",
      ratios == expected_ratios, f"ratios = {sorted(ratios, key=float)}")
# 'not proportional': more than one distinct ratio
check("(2) D is not proportional to perm3^3 (>=2 distinct ratios)",
      len(ratios) >= 2, f"#distinct ratios = {len(ratios)}")

# (3) D irreducible over Q: single factor, degree 9, multiplicity 1, 55 monomials, unit -2
unit, factors = sp.factor_list(D, *VARS)
check("(3) D has a single Q-irreducible factor", len(factors) == 1,
      f"#factors = {len(factors)}")
if len(factors) == 1:
    fac, mult = factors[0]
    check("(3) factor multiplicity = 1", mult == 1, f"mult = {mult}")
    check("(3) factor degree = 9", sp.total_degree(fac) == 9,
          f"deg = {sp.total_degree(fac)}")
    check("(3) factor has 55 monomials", nterms(fac) == 55,
          f"|factor| = {nterms(fac)}")
    check("(3) factorization unit = -2", sp.nsimplify(unit) == -2,
          f"unit = {unit}")

# (4) not a rational cube: unique Q-factor has multiplicity 1 (not div by 3),
#     and the unit -2 is not a rational cube.
mult_ok = (len(factors) == 1 and factors[0][1] % 3 != 0)
unit_not_cube = sp.nsimplify(unit) == -2  # -2 has no rational cube root
check("(4) D is not a perfect cube over Q",
      mult_ok and unit_not_cube,
      "unique factor mult=1 (not div by 3); unit -2 not a rational cube")

# =====================================================================
# Lemma (HEADLINE): gcd of the nine first partials of D is constant
#   => D squarefree over C => D not a cube over C.
# =====================================================================
print("\n=== Lemma : gcd of the nine first partials of D is constant ===")
partials = [sp.expand(sp.diff(D, v)) for v in VARS]
check("nine partials all have degree 8",
      all(sp.total_degree(p) == 8 for p in partials),
      "deg(partials) = {" + ",".join(str(sp.total_degree(p)) for p in partials) + "}")

g = partials[0]
for p in partials[1:]:
    g = sp.gcd(g, p)
g = sp.expand(g)
g_is_const = (sp.total_degree(g) == 0) if g != 0 else False
check("g = gcd(dD/dx11, ..., dD/dx33) is a nonzero CONSTANT",
      g != 0 and g_is_const, f"g = {g}")
# Consequence stated explicitly:
check("=> D is squarefree over C (no repeated factor); NOT a cube over C",
      g != 0 and g_is_const,
      "squarefreeness/logarithmic-derivative criterion, char 0 field-independent")

# =====================================================================
# Remark (independent verification): generic line test
#   D restricted to x = a + t b has 9 distinct roots.
# =====================================================================
print("\n=== Remark : independent generic-line test (9 distinct roots) ===")
t = sp.symbols('t')
a = [sp.Rational(random.randint(-5, 5)) for _ in range(9)]
b = [sp.Rational(random.randint(-5, 5)) for _ in range(9)]
subs = {VARS[i]: a[i] + t * b[i] for i in range(9)}
D_line = sp.expand(D.subs(subs))
D_line_poly = sp.Poly(D_line, t)
check("D restricted to a generic rational line has degree 9 in t",
      D_line_poly.degree() == 9, f"deg_t = {D_line_poly.degree()}")
gg = sp.gcd(D_line_poly, D_line_poly.diff(t))
gg_const = (sp.Poly(gg, t).degree() == 0)
check("gcd(D|line, d/dt D|line) is constant => 9 DISTINCT roots",
      gg_const, "a cube would give all roots of multiplicity 3")

# =====================================================================
# summary
# =====================================================================
print("\n=== SUMMARY ===")
n_ok = sum(1 for _, ok, _ in results if ok)
n_tot = len(results)
print(f"{n_ok}/{n_tot} checks passed.")
if n_ok != n_tot:
    print("FAILED checks:")
    for name, ok, detail in results:
        if not ok:
            print(f"   - {name}  ({detail})")
    sys.exit(1)
print("ALL CLAIMS VERIFIED (exact rational arithmetic).")
sys.exit(0)
