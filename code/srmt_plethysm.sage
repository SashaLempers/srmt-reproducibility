# ================================================================
# SRMT v9: SageMath verification of GL_3 plethysms
# Reference: Proposition (Plethysms for GL_3) in SRMT v9,
#            equations (plethysm-sym) and (plethysm-ext).
#
# SageMath version: 9.6 or 9.7
# Expected runtime: < 5 minutes
# Usage: sage srmt_plethysm.sage
# ================================================================

import sys

print("=" * 60)
print("SRMT v9 — SageMath plethysm verification")
print("SageMath version:", sage.version.version)
print("=" * 60)

# Symmetric functions over QQ
Sym = SymmetricFunctions(QQ)
s = Sym.schur()
h = Sym.homogeneous()
e = Sym.elementary()

# ----------------------------------------------------------------
# Helper: keep only partitions of length <= n (GL_n truncation)
# ----------------------------------------------------------------
def truncate_gln(sf_element, n):
    """
    Project to GL_n: discard summands whose partition has > n parts.
    In SageMath's Schur ring, s[lambda] = 0 for GL_n when
    len(lambda) > n.
    """
    result = {}
    for mono in sf_element.monomial_coefficients().items():
        part, coeff = mono
        if len(part) <= n:
            result[part] = coeff
    return s.sum_of_terms(result.items(), distinct=True)

# ----------------------------------------------------------------
# Step 1: Verify Sym^2(S_{(2,1)}(C^3))
# ----------------------------------------------------------------
print("\n--- Sym^2(S_{(2,1)}(C^3)) ---")

s21 = s[Partition([2, 1])]
h2  = h[Partition([2])]

# Plethysm h_2 o s_{2,1}  (= Sym^2 of S_{2,1})
sym2_full = s(h2.plethysm(s21))
sym2_gl3  = truncate_gln(sym2_full, 3)

print("Sym^2(S_{(2,1)}) in GL_infty:", sym2_full)
print("Truncated to GL_3:           ", sym2_gl3)

expected_sym2 = s[4, 2] + s[3, 2, 1] + s[2, 2, 2]
print("Expected:                    ", expected_sym2)

if sym2_gl3 == expected_sym2:
    print("MATCH ✓")
else:
    print("MISMATCH ✗")
    print("  Got -expected:", sym2_gl3 - expected_sym2)
    sys.exit(1)

# Dimension check
def dim_gl3(partition):
    """
    Dimension of S_lambda(C^3) for partition of length <= 3.
    Uses the hook-content formula.
    """
    a = partition[0] if len(partition) > 0 else 0
    b = partition[1] if len(partition) > 1 else 0
    c = partition[2] if len(partition) > 2 else 0
    return ((a - b + 1) * (b - c + 1) * (a - c + 2)) // 2

dims = [dim_gl3([4,2]), dim_gl3([3,2,1]), dim_gl3([2,2,2])]
total = sum(dims)
print(f"  Dimensions: S_{{4,2}} = {dims[0]}, S_{{3,2,1}} = {dims[1]}, S_{{2,2,2}} = {dims[2]}")
print(f"  Total: {dims[0]} + {dims[1]} + {dims[2]} = {total}")
assert total == 36, f"Dimension error: {total} != 36"
print(f"  Check: C(9,2) = {binomial(9,2)}. ✓")

# ----------------------------------------------------------------
# Step 2: Verify bigwedge^2(S_{(2,1)}(C^3))
# ----------------------------------------------------------------
print("\n--- bigwedge^2(S_{(2,1)}(C^3)) ---")

e2 = e[Partition([2])]

# Plethysm e_2 o s_{2,1}  (= bigwedge^2 of S_{2,1})
ext2_full = s(e2.plethysm(s21))
ext2_gl3  = truncate_gln(ext2_full, 3)

print("bigwedge^2(S_{(2,1)}) in GL_infty:", ext2_full)
print("Truncated to GL_3:                ", ext2_gl3)

expected_ext2 = s[4, 1, 1] + s[3, 3] + s[3, 2, 1]
print("Expected:                         ", expected_ext2)

if ext2_gl3 == expected_ext2:
    print("MATCH ✓")
else:
    print("MISMATCH ✗")
    print("  Got - expected:", ext2_gl3 - expected_ext2)
    sys.exit(1)

dims2 = [dim_gl3([4,1,1]), dim_gl3([3,3]), dim_gl3([3,2,1])]
total2 = sum(dims2)
print(f"  Dimensions: S_{{4,1,1}} = {dims2[0]}, S_{{3,3}} = {dims2[1]}, S_{{3,2,1}} = {dims2[2]}")
print(f"  Total: {dims2[0]} + {dims2[1]} + {dims2[2]} = {total2}")
assert total2 == 28, f"Dimension error: {total2} != 28"
print(f"  Check: C(8,2) = {binomial(8,2)}. ✓")

# ----------------------------------------------------------------
# Step 3: Combined consistency
# ----------------------------------------------------------------
print("\n--- Combined tensor square check ---")
grand_total = total + total2
print(f"dim Sym^2 + dim bigwedge^2 = {total} + {total2} = {grand_total}")
assert grand_total == 64, f"Error: {grand_total} != 64"
print(f"Check: 8^2 = 64 = (dim S_{{(2,1)}}(C^3))^2. ✓")

# ----------------------------------------------------------------
# Step 4: Note on stable vs GL_3 plethysms
# ----------------------------------------------------------------
print("\n--- Stable plethysm (GL_infty) for comparison ---")
print("Sym^2(s_{2,1}) in GL_infty :")
print(" ", sym2_full)
print("bigwedge^2(s_{2,1}) in GL_infty:")
print(" ", ext2_full)
print("(Terms with partitions of length > 3 vanish for GL_3.)")

print("\n" + "=" * 60)
print("All plethysm identities for GL_3 VERIFIED.")
print("Equations (plethysm-sym) and (plethysm-ext) of SRMT v9 confirmed.")
print("=" * 60)
