# SRMT — Pipeline complet GL₉ Young flattening det₃ vs perm₃

Auteur : Sasha Lempers — paper *The Schur-Rank Multiplication Tensor* (v19/v20)
Calculs effectués 2026-05-12, sandbox Linux Python 3.12.

## Synthèse des deux phases

| Phase | Famille | Verdict | Δ rang |
|---|---|---|---|
| **Pieri** | Γ_{λ,μ}(f) : S_λV → S_μV, μ/λ bande horizontale de taille 3 | **B** | 0 sur 37 candidats |
| **Non-Pieri / Koszul-Young** | Y_a(f) : V*⊗Λ^a V → V⊗Λ^{a+1} V (LO11 §1.2 Rem 4.1.3) | **A** | **16 à a=4** |

**Résultat principal** : rank Y_4(det₃) = 950, rank Y_4(perm₃) = 934 sur Q exact.

## Navigation des fichiers

### Phase Pieri (Verdict B)

- `step1_detperm.py` / `step1_detperm.json` — det₃, perm₃ comme vecteurs Q^165 (base monomiale multiset)
- `step2_schur.py` — partitions, dim_Schur, SSYT, hook length, validations Weyl
- `step3_pieri.py` — règle de Pieri (horizontal strips)
- `step3b_gamma.py` — construction Γ_{λ,μ}(f) en base SSYT
- `step4_main.py` — boucle Bareiss/Fraction initiale (trop lente)
- `step4b_fast_modular.py` — optimisation numpy int64 dense mod p
- `step4c_full.py` — boucle complète sur 37 candidats mod p=65521
- `step7_verdict.py` — agrégation verdict Pieri
- `step8_bareiss.py` — Bareiss exact Q (21 candidats prod ≤ 5·10⁶)
- `results_full.json` — ranks modulaires des 37 candidats Pieri
- `bareiss_results.json` — ranks Q exacts des 21 candidats faisables
- `final_verdict.json` — verdict B consolidé phase Pieri

### Phase Non-Pieri / Koszul-Young (Verdict A) ← RÉSULTAT PRINCIPAL

- `np1_koszul_basic.py` — construction dérivées partielles ∂f/∂x_i, validation det/perm
- `np1_partials.pkl` — dérivées sauvegardées
- `np2_koszul_matrix.py` — **construction Y_a(f) LO11 §1.2 Rem 4.1.3**, calcul rang mod p pour a∈{1,2,3,4}
- `np2_koszul_principal.json` — résultats ranks modulaires (révèle differs à a=4)
- `np3_certif_a4.py` — re-construction matrices a=4 + multi-prime (4 primes)
- `np3_M_det_a4.npy`, `np3_M_perm_a4.npy` — **matrices 1134×1134 binaires** (preuves brutes)
- `np3_multiprime.json` — ranks sur 4 primes (10⁹+7, 998244353, 10⁹+9, 65521)
- `np4_bareiss_a4.py` — première tentative Bareiss (BUGGUÉE — invariant cassé par skips)
- `np4b_debug_bareiss.py` — diagnostic du bug (4 primes supplémentaires + numpy float)
- `np4c_bareiss_v2.py` — Bareiss avec détection violation invariant (confirme bug)
- `np4d_bareiss_v3.py` — **Gauss-Jordan sparse Fraction exact Q** (correct, 0.3s)
- `np4d_gauss_exact.json` — rangs Q exacts : 950 et 934
- `np5_validate_koszul.py` — **4 tests de validation construction Koszul**
  - Y_4(x₀³) = 70 = C(8,4) ✓
  - Y_4(x₀³+x₁³) = 140 ✓
  - GL₉-équivariance sous swap variables ✓
  - Cubique aléatoire → 1134 (full) ✓
- `np6_final_verdict.py` / `np6_final_verdict.json` — **VERDICT A** consolidé

### Métadonnées

- `sha256sums.txt` — empreintes de tous les fichiers de preuve
- `audit.json` — environnement sandbox
- `biblio_search.json` — recherche bibliographique LO11 / LO13 / Han-Ju-Kim

## Reproduction

```bash
cd srmt/
python step1_detperm.py          # base de calcul
python np1_koszul_basic.py        # dérivées
python np2_koszul_matrix.py       # scan a=1..4 modulaire
python np3_certif_a4.py           # multi-prime certif a=4
python np4d_bareiss_v3.py         # Gauss exact Q
python np5_validate_koszul.py     # 4 tests
python np6_final_verdict.py       # verdict A consolidé
```

## Chaîne de certification a=4

1. **numpy.linalg.matrix_rank (float64)** : 950 / 934
2. **8 primes différents** (mod p) : 950 / 934 — concordance totale
3. **Gauss-Jordan sparse Fraction (Q exact)** : 950 / 934
4. **Tests de sanity** : décomposable rang sym 1 → Y_4 rang 70 = C(8,4), confirmé par LO11 Prop 4.1.1
5. **GL₉-équivariance** : rang invariant sous permutation des variables

## Référence canonique

**Landsberg, J.M. ; Ottaviani, G.** (2011) *Equations for secant varieties of Veronese and other varieties*, arXiv:1111.4567, Section 1.2 Remarque 4.1.3.

Pour V = ℂⁿ, d = 2δ+1 impair, a = ⌊n/2⌋ :
$$YF_{d,n}(\\phi) : S^\\delta V^* \\otimes \\Lambda^a V \\;\\longrightarrow\\; S^\\delta V \\otimes \\Lambda^{a+1}V$$

Pour d=3, n=9 : δ=1, a=4 → matrice 1134×1134 qui sépare det₃ et perm₃ par Δ=16.
