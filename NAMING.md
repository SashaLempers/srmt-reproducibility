# Correspondence between paper script names and repository files

The reproducibility section of the manuscript refers to idealized script
names. The actual implementation in this repository is organized differently.
This table maps the names cited in the paper to the files that actually
perform each computation.

| Name cited in paper            | Actual file(s) in this repository                         | Status |
|--------------------------------|-----------------------------------------------------------|--------|
| `plethysm.py`                  | `code/srmt_plethysm.sage`, `code/m2/piste1_option_B_schurrings.m2` | present |
| `repth_tools.py`               | helper functions inside `code/python/srmt_v19_*.py`       | present (inlined) |
| `catalecticant_ranks.py`       | `code/m2/srmt_v15_pipeline.m2` (catalecticant ranks)      | present |
| `build_hw_vectors.py`          | `code/python/srmt_v19_sparse_search.py`                   | present |
| `rank_test.py`                 | `code/python/srmt_v19_certify.py` (exact ℚ jet test)      | present |
| `koszul_young_test.py`         | `code/python/route_flattening_search.py` (C1–C2 complete; C3–C5 are skeletons) | partial |
| `hessian_step_B1_B3.py`        | — no dedicated script —                                    | **MISSING** |
| `hessian_step_B4b.py`          | — no dedicated script —                                    | **MISSING** |
| `hessian_step_B5.py`           | — no dedicated script —                                    | **MISSING** |

## Known gaps

1. **Hessian squarefree lemma.** The lemma asserting that `det Hess(perm₃)`
   is squarefree over C, that the gcd of its nine first partials is constant,
   and that it is "not a cube", has no dedicated exact computation here. In the
   current code the Hessian appears only as mod-p numeric rank sampling and as
   the C1 catalecticant in the flattening script. An exact discriminant /
   factorisation script over ℚ should be added (or the claim re-tagged as
   purely algebraic).

2. **Flattening defects C3–C5.** `route_flattening_search.py` documents C1 and
   C2 as complete; C3, C4 and C5 are skeleton stubs with TODOs.

3. **`m2_data.txt`.** Required by the Python core; regenerate from the M2
   pipeline (see README).
