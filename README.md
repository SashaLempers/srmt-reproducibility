# SRMT Reproducibility — Bilateral GL₉-Separation of det₃ and perm₃

Reproducibility repository for the paper

> **Bilateral GL₉-Separation of det₃ and perm₃ via Koszul–Young Pieri Flattening, with Degree-Four Ideal Vanishing and Homological Obstructions** (v23, May 2026)
> Sasha Lempers — Independent Researcher, Annecy-le-Vieux, France
> Contact: lemperssasha@gmail.com

This repository contains the computational artifacts backing the certified claims of the paper.
All numerical values cited in the paper trace to a file here; the mapping from audit keys to files is given below.

## Headline result

A GL₉-equivariant Koszul–Young Pieri flattening `Y(π₃)` with `π₃ = (2⁴,1⁴)` separates the two cubics:

| quantity | det₃ | perm₃ | Δ |
|---|---|---|---|
| `rank Y(π₃)(f)` (1134×1134, exact over ℚ) | **950** | **934** | **16** |

This is the direction `det₃ ∉ GL₉·perm₃` (Theorem Y). Combined with the singular-locus obstruction
`σ_G(perm₃)=3 < 5=σ_G(det₃)` (Theorem G, direction `perm₃ ∉ GL₉·det₃`), it gives the **bilateral**
non-containment on Sym³(C⁹) in the non-padded case m = n = 3.

The 950/934 result is in [`pipeline/koszul_young/np6_final_verdict.json`](pipeline/koszul_young/np6_final_verdict.json),
cross-checked modulo four primes (65521, 998244353, 1000000007, 1000000009) — all agree.

## Repository layout

```
pipeline/koszul_young/   Real Koszul–Young flattening pipeline (Python, exact ℚ + modular)
                         np1..np6  : Koszul–Young Y_4 construction → ranks 950/934
                         step1..step8 : Pieri small-|λ| scan (37 candidates)
                         *.json, *_log.txt : results and logs
                         sha256sums.txt : integrity manifest
macaulay2/               Macaulay2 audit scripts (homology, plethysm/HWV, Hessian, apolarity)
audit/                   The six audit-key entries cited in the paper (see mapping below)
```

## Audit-key → file mapping

The paper cites six audit log entries. Each corresponds to a file in `audit/`:

| Audit key (in paper) | File | Backing artifact(s) | Value(s) |
|---|---|---|---|
| `bareiss_results_21` | [`audit/bareiss_results_21.json`](audit/bareiss_results_21.json) | `pipeline/koszul_young/final_verdict.json`, `bareiss_results.json`, `step8_bareiss.py` | 21 QQ-CERTIFIED Bareiss exact-ℚ rank entries (Pieri scan) |
| `mod_results_16` | [`audit/mod_results_16.json`](audit/mod_results_16.json) | `pipeline/koszul_young/final_verdict.json`, `step4b_fast_modular.py` | 16 MOD-CERTIFIED entries over GF(65521) |
| `pdim_detm_4` | [`audit/pdim_detm_4.json`](audit/pdim_detm_4.json) | `macaulay2/srmt_gct_combined_lab.m2` | pdim det₃ = 4, pdim perm₃ = 9 (CM dichotomy) |
| `hwv_vlambda_det_12` | [`audit/hwv_vlambda_det_12.json`](audit/hwv_vlambda_det_12.json) | `macaulay2/srmt_gct_combined_lab.m2` | 12-isotype plethysm s₃[s₃]; HWV non-vanishing on det₃ ⇒ I₄ = 0 |
| `hessrank_HL_9_8_9` | [`audit/hessrank_HL_9_8_9.json`](audit/hessrank_HL_9_8_9.json) | `macaulay2/srmt_v19_pipeline.m2`, `srmt_v18_pipeline.m2` | ρ(det₃)=9 (deposited); ρ(P₁)=8, ρ(P₂)=9 EXACT-COMPUTED (see note) |
| `apolar_hf_det_perm` | [`audit/apolar_hf_det_perm.json`](audit/apolar_hf_det_perm.json) | `macaulay2/srmt_v19_verify.m2`, `srmt_v19_pipeline.m2` | H_apolar(det₃) = H_apolar(perm₃) = [1,9,9,1] |

**Note on `hessrank_HL_9_8_9`.** The deposited Macaulay2 Hessian diagnostic certifies ρ(det₃) = 9
(generic Hessian of the determinant). The Hüttenhain–Lairez boundary-point values ρ(P₁) = 8 and
ρ(P₂) = 9 are stated in the paper at status EXACT-COMPUTED and are reproducible from the HL
boundary-point construction; the explicit boundary-point computation is not part of the deposited
diagnostic scripts. The boundary defect (Pair B = (det₃, P₁)) is a corroborating computation in a
pair distinct from the main pair (det₃, perm₃) and is not load-bearing for the bilateral Main Theorem.

## Proof-status taxonomy (as in the paper)

- **PROVED-ALG** — symbolic algebraic proof, no computer arithmetic.
- **QQ-CERTIFIED** — exact ℚ computation in Macaulay2/Sage.
- **EXACT-COMPUTED** — exact ℚ computation, independently reproduced by ≥ 2 methods.
- **OBS-N** — numerical observation; not used for any main claim.

## Reproducing

Python pipeline (Koszul–Young ranks + Pieri scan):

```bash
cd pipeline/koszul_young
python3 np1_koszul_basic.py && python3 np2_koszul_matrix.py && python3 np3_certif_a4.py
python3 np4_bareiss_a4.py && python3 np5_validate_koszul.py && python3 np6_final_verdict.py
# Pieri scan:
python3 step1_detperm.py ... step8_bareiss.py
```

Macaulay2 audit scripts (require Macaulay2 ≥ 1.19 with the `SchurRings` package):

```bash
M2 --script macaulay2/srmt_gct_combined_lab.m2   # pdim + plethysm/HWV
M2 --script macaulay2/srmt_v19_verify.m2         # apolar Hilbert functions
M2 --script macaulay2/srmt_v19_pipeline.m2       # Hessian diagnostic, singular-locus dims
```

## References

- J.M. Landsberg, G. Ottaviani, *Equations for secant varieties of Veronese and other varieties* (2013).
- C. Farnsworth (2015) — equivariant Pieri-type flattenings; the 950/934 values.
- S. Sam, *PieriMaps* Macaulay2 package (2009).
- J. Hüttenhain, P. Lairez — boundary of the determinant orbit closure.

## License

See [LICENSE](LICENSE).
