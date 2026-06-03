# SRMT Reproducibility — Equivariant obstructions in Sym³(C⁹)

Reproducibility repository for the paper

> **Equivariant obstructions in Sym³(C⁹): an unconditional one-sided Hessian obstruction perm₃ ∉ GL₉·det₃, and an exhaustive low-degree SRMT certification for det₃ vs. perm₃**
> Sasha Lempers — Independent researcher, Annecy-le-Vieux, France
> Contact: lemperssasha@gmail.com

This repository contains the computational artifacts backing the certified claims of the paper.
The file [`NAMING.md`](NAMING.md) maps each routine named in the manuscript to its source file here.

## Headline result

The paper proves, **unconditionally and one-sided**, the GCT-irrelevant non-containment

> perm₃ ∉ closure(GL₉·det₃)   (Hessian obstruction, PROVED-ALG)

via the degree-9 Hessian map H̄(f) = det Hess(f): one has `det Hess(det₃) = −2·(det₃)³`,
whereas `det Hess(perm₃)` is **squarefree and not a cube over C** (the gcd of its nine first
partials is the nonzero constant `2`). No GL₉ separation is claimed in the reverse
(GCT-relevant) direction det₃ ∈ closure(GL₉·perm₃), which is left **open**.

In addition, an exhaustive degree-≤5 SRMT scan finds **no separator**: across all candidate
Pieri/Koszul–Young flattenings in that range the ranks of det₃ and perm₃ **coincide**
(see `audit/bareiss_results_21.json`, `audit/mod_results_16.json`). In particular,
first-order Koszul/Young flattenings do **not** separate this pair at size m = n = 3.

## Repository layout

```
code/python/    Exact rational (SymPy) scripts: Hessian certification, jet certification,
                degree-≤5 sweep, flattening search
code/m2/        Macaulay2 scripts: GL₃×GL₃ catalecticant ranks, SL₉-invariant multiplicities
code/           srmt_plethysm.sage (plethysm s₃[s₃])
logs/           Run logs and machine-readable certificates
macaulay2/      Macaulay2 audit scripts (homology, plethysm/HWV, Hessian, apolarity)
audit/          The six audit-key entries cited in the paper (see mapping below)
archive_v18_v19/  Earlier-version helper scripts and data (not cited by the current paper)
```

## Audit-key → file mapping

The paper cites six audit log entries. Each corresponds to a file in `audit/`:

| Audit key (in paper) | File | Backing artifact(s) | Value(s) |
|---|---|---|---|
| `bareiss_results_21` | [`audit/bareiss_results_21.json`](audit/bareiss_results_21.json) | self-contained (21 entries embedded); `code/python/srmt_v19_certify.py` | 21 QQ-CERTIFIED Bareiss exact-ℚ rank entries (degree-≤5 scan); det₃ and perm₃ ranks **coincide** on all 21 (EQUAL_Q) |
| `mod_results_16` | [`audit/mod_results_16.json`](audit/mod_results_16.json) | self-contained (16 entries embedded); `code/python/srmt_v19_sparse_search.py` | 16 MOD-CERTIFIED entries over GF(65521); det₃ and perm₃ ranks **coincide** |
| `pdim_detm_4` | [`audit/pdim_detm_4.json`](audit/pdim_detm_4.json) | `macaulay2/srmt_gct_combined_lab.m2` | pdim det₃ = 4, pdim perm₃ = 9 (CM dichotomy) |
| `hwv_vlambda_det_12` | [`audit/hwv_vlambda_det_12.json`](audit/hwv_vlambda_det_12.json) | `macaulay2/srmt_gct_combined_lab.m2` | 12-isotype plethysm s₃[s₃]; HWV non-vanishing on det₃ ⇒ I₄ = 0 |
| `hessrank_HL_9_8_9` | [`audit/hessrank_HL_9_8_9.json`](audit/hessrank_HL_9_8_9.json) | `macaulay2/srmt_v19_pipeline.m2`, `srmt_v18_pipeline.m2` | ρ(det₃)=9 (deposited); ρ(P₁)=8, ρ(P₂)=9 EXACT-COMPUTED (see note) |
| `apolar_hf_det_perm` | [`audit/apolar_hf_det_perm.json`](audit/apolar_hf_det_perm.json) | `macaulay2/srmt_v19_verify.m2`, `srmt_v19_pipeline.m2` | H_apolar(det₃) = H_apolar(perm₃) = [1,9,9,1] |

**Note on `hessrank_HL_9_8_9`.** The deposited Macaulay2 Hessian diagnostic certifies ρ(det₃) = 9
(generic Hessian of the determinant). The Hüttenhain–Lairez boundary-point values ρ(P₁) = 8 and
ρ(P₂) = 9 are stated in the paper at status EXACT-COMPUTED and are reproducible from the HL
boundary-point construction; the explicit boundary-point computation is not part of the deposited
diagnostic scripts. The boundary defect (Pair B = (det₃, P₁)) is a corroborating computation in a
pair distinct from the main pair (det₃, perm₃) and is not load-bearing for the Main Theorem.

## Proof-status taxonomy (as in the paper)

- **PROVED-ALG** — symbolic algebraic proof, no computer arithmetic.
- **QQ-CERTIFIED** — exact ℚ computation in Macaulay2/Sage.
- **EXACT-COMPUTED** — exact ℚ computation, independently reproduced by ≥ 2 methods.
- **OBS-N** — numerical observation; not used for any main claim.

## Reproducing

Hessian certification (exact rational, SymPy):

```bash
python3 code/python/hessian_perm3.py     # exit code 0 iff all 16 checks pass
```

Degree-≤5 SRMT scan (exact ℚ jet-certification + sweep):

```bash
python3 code/python/srmt_v19_certify.py
python3 code/python/srmt_v19_sparse_search.py
```

Macaulay2 audit scripts (require Macaulay2 ≥ 1.19 with the `SchurRings` package):

```bash
M2 --script macaulay2/srmt_gct_combined_lab.m2   # pdim + plethysm/HWV
M2 --script macaulay2/srmt_v19_verify.m2         # apolar Hilbert functions
M2 --script macaulay2/srmt_v19_pipeline.m2       # Hessian diagnostic, singular-locus dims
```

## References

- J.M. Landsberg, G. Ottaviani, *Equations for secant varieties of Veronese and other varieties* (2013).
- S. Sam, *PieriMaps* Macaulay2 package (2009).
- J. Hüttenhain, P. Lairez — boundary of the determinant orbit closure.

## License

See [LICENSE](LICENSE).
