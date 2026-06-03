# Correspondence between paper script names and repository files

The manuscript refers to idealized script names. This table maps each
name cited in the paper to the file(s) that actually perform the
computation in this repository.

| Name cited in paper            | Actual file(s) in this repository                         | Status |
|--------------------------------|-----------------------------------------------------------|--------|
| `plethysm.py`                  | `code/srmt_plethysm.sage`, `code/m2/piste1_option_B_schurrings.m2` | present |
| `repth_tools.py`               | helper functions inside `code/python/srmt_v19_*.py`       | present (inlined) |
| `catalecticant_ranks.py`       | `code/m2/srmt_v15_pipeline.m2` (catalecticant ranks)      | present |
| `build_hw_vectors.py`          | `code/python/srmt_v19_sparse_search.py`                   | present |
| `rank_test.py`                 | `code/python/srmt_v19_certify.py` (exact Q jet test)      | present |
| `koszul_young_test.py`         | `code/python/route_flattening_search.py` (C1-C2 complete; C3-C5 skeletons) | partial |
| `hessian_step_B1_B3.py`        | `code/python/hessian_perm3.py` (det Hess(det3) baseline + perm3 setup) | present |
| `hessian_step_B4b.py`          | `code/python/hessian_perm3.py` (residual R, 37 monomials, ratios, factor_list) | present |
| `hessian_step_B5.py`           | `code/python/hessian_perm3.py` (gcd-of-partials squarefree lemma + line test) | present |

## Hessian certification

`code/python/hessian_perm3.py` is an exact rational (SymPy) script that
verifies, with no floating point:

- `det Hess(det3) = -2 (det3)^3`                                 (Thm. Hess-det)
- residual `R = det Hess(perm3) + 2 perm3^3` has exactly 37 monomials   (Thm. Hess-perm (1))
- distinct coefficient ratios `D/perm3^3` on common support = `{-2, 2/3, 10/3}` (Thm. Hess-perm (2))
- `D = det Hess(perm3)` is Q-irreducible: one factor, degree 9, multiplicity 1, 55 monomials, unit -2 (Thm. Hess-perm (3))
- `D` is not a rational cube                                     (Thm. Hess-perm (4))
- **gcd of the nine first partials of `D` is the nonzero constant `2`**,
  hence `D` is squarefree over C and NOT a cube over C           (Lemma squarefree, HEADLINE)
- generic-line restriction of `D` has 9 distinct roots (independent check) (Remark)

Run it with `python3 code/python/hessian_perm3.py` (exit code 0 iff all
16 checks pass). The machine-readable result is in
`logs/hessian/hessian_perm3_certificate.json`.

## Earlier-version helper scripts

The directory `archive_v18_v19/` contains helper and data-generation scripts
from earlier versions of the project (e.g. `gen_m2_data.py`, `gen_hessian_cert.py`,
`srmt_v17_pipeline.m2`, `route_alpha_perm3.py`, `m2_data.txt`). They are not
cited by the current paper and are retained for provenance only.

## Remaining partial item

- **Flattening defects C3-C5.** `code/python/route_flattening_search.py` documents
  C1 and C2 as complete; C3, C4 and C5 are skeleton stubs with TODOs.
