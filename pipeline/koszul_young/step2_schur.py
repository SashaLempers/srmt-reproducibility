"""
Etape 2 - Infrastructure Schur:
  * Partitions, conjugue, longueur de hook, dim(S_lambda(C^n)) via hook content formula.
  * Semi-standard Young tableaux (SSYT) sur l'alphabet {1..n} : base de S_lambda(C^n).
  * Standard Young tableaux (SYT) : dim irrep S_d.
  * Pieri rule : S_lambda * Sym^k = sum over mu obtained by adding k boxes in horizontal strip.
  * Validations: dim hook = nombre de SSYT; hook formula classique = nombre de SYT.

On garde tout en entier (Fraction si nécessaire).
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial
from itertools import product

# --- Partitions ---
def partitions(n, max_len=None):
    """Toutes les partitions de n, de longueur <= max_len."""
    res = []
    def rec(remaining, max_part, current):
        if remaining == 0:
            if max_len is None or len(current) <= max_len:
                res.append(tuple(current))
            return
        for p in range(min(remaining, max_part), 0, -1):
            if max_len is not None and len(current) >= max_len: break
            current.append(p)
            rec(remaining - p, p, current)
            current.pop()
    rec(n, n, [])
    return res

def conjugate(lam):
    if not lam: return ()
    return tuple(sum(1 for p in lam if p >= i+1) for i in range(lam[0]))

# --- Dimension irrep S_d associee a lambda (hook length formula) ---
def hook_length(lam, i, j):
    """Hook length at cell (i,j), 0-indexed."""
    lam_conj = conjugate(lam)
    return (lam[i] - j - 1) + (lam_conj[j] - i - 1) + 1

def dim_Sd(lam):
    """Nombre de SYT = dim de l'irrep de S_n indexee par lambda."""
    n = sum(lam)
    num = factorial(n)
    denom = 1
    for i in range(len(lam)):
        for j in range(lam[i]):
            denom *= hook_length(lam, i, j)
    assert num % denom == 0
    return num // denom

# --- Dimension de S_lambda(C^n) (Weyl dimension formula via hook content) ---
def dim_Schur(lam, n):
    """dim S_lambda(C^n) = prod over cells (i,j) of (n + j - i) / hook(i,j)."""
    if len(lam) > n: return 0
    num = 1
    denom = 1
    for i in range(len(lam)):
        for j in range(lam[i]):
            num *= (n + j - i)
            denom *= hook_length(lam, i, j)
    assert num % denom == 0, (lam, n, num, denom)
    return num // denom

# --- Semi-standard Young tableaux (SSYT) ---
# Une SSYT de forme lam sur {1..n}: entries faiblement croissantes le long des lignes,
# strictement croissantes le long des colonnes.

def SSYT(lam, n):
    """Genere toutes les SSYT de forme lam sur l'alphabet 1..n."""
    if not lam:
        yield ()
        return
    cells = [(i,j) for i in range(len(lam)) for j in range(lam[i])]
    # Ordre lexico (i,j)
    def fill(idx, T):
        if idx == len(cells):
            yield tuple(tuple(row) for row in T)
            return
        i, j = cells[idx]
        # min value: max(T[i][j-1], T[i-1][j]+1, 1)
        lo = 1
        if j > 0: lo = max(lo, T[i][j-1])
        if i > 0: lo = max(lo, T[i-1][j] + 1)
        for v in range(lo, n+1):
            T[i].append(v) if j == len(T[i]) else None
            if j < len(T[i]):
                T[i][j] = v
            else:
                T[i].append(v)
            yield from fill(idx + 1, T)
            T[i].pop()
    # init T comme listes
    T = [[] for _ in range(len(lam))]
    yield from fill(0, T)

def count_SSYT(lam, n):
    return sum(1 for _ in SSYT(lam, n))

def SYT(lam):
    """Generate standard Young tableaux of shape lam, alphabet 1..n with n=|lam|, each value used once."""
    n = sum(lam)
    cells = [(i,j) for i in range(len(lam)) for j in range(lam[i])]
    # Filling: assign 1..n to cells such that increasing in rows and columns
    def fill(values_left, T, positions_filled):
        if not values_left:
            yield tuple(tuple(row) for row in T)
            return
        v = min(values_left)
        # v doit aller dans une case (i,j) telle que (i,j-1) et (i-1,j) deja remplies
        for (i,j) in cells:
            if (i,j) in positions_filled: continue
            if j > 0 and (i,j-1) not in positions_filled: continue
            if i > 0 and (i-1,j) not in positions_filled: continue
            T[i][j] = v
            positions_filled.add((i,j))
            yield from fill(values_left - {v}, T, positions_filled)
            positions_filled.remove((i,j))
            T[i][j] = 0
    T = [[0]*l for l in lam]
    yield from fill(set(range(1,n+1)), T, set())

# --- Tests ---
if __name__ == "__main__":
    print("=== Tests partitions ===")
    for n in range(1,7):
        ps = partitions(n)
        print(f"|p({n})| = {len(ps)}: {ps}")
    print()
    print("=== Tests dim Schur S_lambda(C^9) ===")
    # Reference values via Weyl dimension formula
    refs = {
        ((1,), 9): 9,
        ((2,), 9): 45,
        ((1,1), 9): 36,
        ((3,), 9): 165,
        ((2,1), 9): 240,
        ((1,1,1), 9): 84,
        ((4,), 9): 495,
        ((3,1), 9): 990,
        ((2,2), 9): 540,
        ((2,1,1), 9): 630,
        ((1,1,1,1), 9): 126,
    }
    for (lam, n), expected in refs.items():
        got = dim_Schur(lam, n)
        ok = "OK" if got == expected else "FAIL"
        print(f"  dim S_{lam}(C^{n}) = {got}  (ref {expected})  {ok}")
        assert got == expected, (lam, got, expected)
    print()
    print("=== Tests dim S_lambda(C^n) == #SSYT(lam, n) ===")
    for lam in [(1,),(2,),(1,1),(3,),(2,1),(1,1,1),(2,2),(3,1),(2,1,1)]:
        for n in [3, 5]:
            d_weyl = dim_Schur(lam, n)
            d_count = count_SSYT(lam, n)
            ok = "OK" if d_weyl == d_count else "FAIL"
            print(f"  lam={lam} n={n}: Weyl={d_weyl} SSYT_count={d_count} {ok}")
            assert d_weyl == d_count

    print()
    print("=== Tests dim_Sd via hook length (nombre de SYT) ===")
    syt_refs = {(1,):1,(2,):1,(1,1):1,(3,):1,(2,1):2,(1,1,1):1,(2,2):2,(3,1):3,(4,):1,(1,1,1,1):1,(3,2):5,(2,2,1):5,(4,1):4,(3,1,1):6,(2,1,1,1):4,(5,):1}
    for lam, ref in syt_refs.items():
        got = dim_Sd(lam)
        syt_count = sum(1 for _ in SYT(lam))
        assert got == ref == syt_count, (lam, got, ref, syt_count)
        print(f"  lam={lam}: hook={got} SYT_count={syt_count} ref={ref} OK")
    print("\nStep 2 infrastructure OK")
