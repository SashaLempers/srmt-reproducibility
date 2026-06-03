# SRMT — Reproducibility code for "Equivariant obstructions in Sym³(C⁹)"

Author: **Sasha Lempers** — independent researcher.
Contact: lemperssasha@gmail.com

This repository contains the computational scripts and certificates that
support the computational claims of the manuscript on equivariant
obstructions and the symmetric rank / apolarity analysis of `det₃` and
`perm₃` inside Sym³(C⁹).

## Layout

```
code/
  python/
    srmt_v19_full.py            Core pipeline: SRMT sparse search + driver
    srmt_v19_certify.py         Exact ℚ jet-certification testJet(P,4)
    srmt_v19_sparse_search.py   Sparse gl₉ search (reads code/m2_data.txt)
    route_alpha_perm3.py        Exploratory GPU / mod-p degree-5 search
    route_flattening_search.py  Catalecticant / Young-flattening defect search
  m2/
    srmt_v15_pipeline.m2        LOCK & AUDIT release; generates the
                                lock_audit certificate (Macaulay2)
    srmt_v17_pipeline.m2        Combined GCT lab (standalone)
    piste1_option_B_schurrings.m2  SL₉-invariant multiplicity via plethysm
  srmt_plethysm.sage            GL₃ plethysm verification (SageMath)
logs/
  lock_audit/
    srmt_v19_lock_audit_certificate.json   QQ-certified results + provenance
LICENSE
```

## Required external data

`code/python/srmt_v19_sparse_search.py` (and the two scripts that import it)
expect a file **`code/m2_data.txt`** that is dumped by the Macaulay2 pipeline
(sparse gl₉ generators, monomial supports, vector of det₃). This file is **not
yet included** — regenerate it from the M2 pipeline before running the Python
core. See `NAMING.md` for the correspondence between the script names used in
the paper and the actual files here.

## Environment

- Macaulay2 ≥ 1.19.1 (certificates produced with seed 42, prime 32003)
- Python ≥ 3.10
- SageMath (for the `.sage` plethysm check)

## License

See `LICENSE`.
