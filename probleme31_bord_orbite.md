# Problème 31 — Bord de l'orbite GL₉·perm₃ et coïncidence des multiplicités

**Auteur :** Sasha Lempers · Annecy-le-Vieux · juin 2026  
**Contact :** lemperssasha@gmail.com  
**Statut :** Document de travail — *à valider avant tout pivot de stratégie*

---

## Légende des statuts épistémiques

> **[PROUVÉ]** = théorème établi en littérature ou preuve fournie ici  
> **[CALCULÉ]** = vérifié numériquement/finiment, pas de preuve générale  
> **[CONJECTURÉ]** = hypothèse de travail, non prouvée  
> **[LACUNE]** = question identifiée, ouverte

---

## 1. Contexte et mission

L'étape précédente (voir `probleme31_synthese_dimensions_superieures.md`) a établi :

- **[CALCULÉ]** d=3 (GL₉) : m_perm(λ) ≥ m_det_full(λ) pour les 30 partitions de 9 (10 avec m_det_full > 0, max = 2, cas serré : égalité sur [2,2,1,1,1,1,1]).
- **[CALCULÉ]** d=4 (GL₁₆) : m_perm(λ) ≥ m_det_full(λ) pour les 185 partitions de 16 à m_det_full > 0 (max = 14, cas serré : égalité sur [1¹⁶]).
- **[PROUVÉ, uniquement pour dim V_λ = 1]** : m_adhérence = m_orbite (collapse lemma).
- **[LACUNE]** : Pour les irreps de dim > 1, les résultats ci-dessus sont des *multiplicités d'orbite*. Le bord ∂(GL₉·perm₃) = adhérence(GL₉·perm₃) ∖ GL₉·perm₃ pourrait en principe contribuer des multiplicités supplémentaires.

**Mission (Option B) :** Démontrer ou réfuter que les composantes du bord ne créent aucune nouvelle obstruction pour une famille identifiée d'irreps λ, ou établir que m_adhérence = m_orbite pour une telle famille.

---

## 2. Premier sous-objectif : structure du bord ∂(GL₉·perm₃)

### 2.1. Résultats connus sur la géométrie de X_perm

#### Dimension de l'orbite

**[PROUVÉ]** Par Kumar [1], Corollaire 5.8 : pour la permanente padded p = x_{1,1}^{m-n} · perm_n dans Sym^m(End(v)^*), avec n = 3 et m = n² = 9 :

$$\dim X_p = m^2(n^2 + 1) - 2n + 1$$

Dans notre cas (perm₃ dans Sym³(ℂ⁹), GL₉ agissant) : la dimension de l'orbite est

$$\dim(GL_9 \cdot \mathrm{perm}_3) = 81 - \dim H_{\mathrm{perm}_3}$$

où H_{perm₃} est le stabilisateur. D'après Proposition 4.1 de Kumar [1], pour n ≥ 3, l'isotropie de perm_n sous GL(End v₁) est le groupe OD des matrices de la forme X ↦ αXβ ou X ↦ αX^t β avec α, β dans le groupe engendré par les matrices de permutation et les matrices diagonales de déterminant 1, ce qui donne pour n=3 un groupe **d'ordre 72** (fini, comme attendu).

Ainsi :
$$\dim(GL_9 \cdot \mathrm{perm}_3) = 81 - 0 = 81$$
(stabilisateur fini ⇒ dimension de l'orbite = dimension du groupe GL₉ = 81)

*Remarque :* dim(Sym³(ℂ⁹)) = C(9+2,3) = 165. L'orbite GL₉·perm₃ est une sous-variété de dimension 81 dans un espace de dimension 165, donc de codimension **84**.

#### Non-normalité de l'adhérence

**[PROUVÉ]** Kumar, Théorème 4.5 [1] : Pour n ≥ 3, la sous-variété GL(End v₁)·perm_n ⊂ Sym^n((End v₁)^*) est **non normale**.

Conséquence immédiate pour notre problème : X_{perm₃} = adhérence(GL₉·perm₃) est **non normale**. En particulier, l'anneau de coordonnées ℂ[X_{perm₃}] est strictement contenu dans sa normalisation ℂ[X̃_{perm₃}], et il n'y a pas égalité ℂ[adhérence] = ℂ[orbite].

**Cette non-normalité est un obstacle fondamental : le théorème de Hartogs algébrique NE s'applique PAS.**

#### Codimension du bord

**[PROUVÉ]** Kumar, Corollaire 5.6 + Lemme 5.7 [1] : Pour 3 ≤ n < m, chaque composante irréductible du bord est de **codimension 1** dans l'adhérence. Par la réduction Theorem 5.2 de Kumar (la décomposition en GL-modules de ℂ[X_p] se ramène à celle de GL(S)·p), ce résultat s'applique à la variété complète.

**Conclusion : le bord ∂(GL₉·perm₃) a ses composantes irréductibles de codimension 1.** La piste « codimension ≥ 2 implique Hartogs » échoue à la base.

### 2.2. Ce que l'on NE sait PAS sur le bord de perm₃

**[LACUNE]** Classification des composantes irréductibles du bord ∂(GL₉·perm₃) : **inconnue dans la littérature**. Comparaison avec le déterminant :

- Pour det₃ : Hüttenhain–Lairez [2] et la thèse de Hüttenhain [3] établissent que ∂(GL₉·det₃) est une union de **exactement deux** clôtures d'orbites, toutes deux de codimension 1. La classification est complète.
- Pour perm₃ : **aucun résultat analogue n'existe dans la littérature.** La classification du bord de l'orbite de la permanente est réputée ouverte.

---

## 3. Deuxième sous-objectif : coïncidence m_adhérence = m_orbite

### 3.1. Le collapse lemma (rappel et limite)

**[PROUVÉ]** *Collapse lemma (dim V_λ = 1)* : Si dim V_λ = 1, toute invariante relative de poids λ est, en restriction à l'orbite dense GL₉·perm₃, un scalaire fois la seule section. L'adhérence ne peut pas rajouter de multiplicité : m_adhérence(λ) = m_orbite(λ).

**[LACUNE]** Ce raisonnement utilise essentiellement que les sections sur l'orbite n'ont qu'une direction possible (le cas dim V_λ = 1), ce qui ne se généralise pas directement à dim > 1.

### 3.2. Piste A : prolongement de Hartogs algébrique — FERMÉE

**Idée initiale :** Si le bord était de codimension ≥ 2, le théorème de Hartogs algébrique (toute section régulière sur l'orbite se prolonge au complémentaire d'une sous-variété de codim ≥ 2 dans une variété normale) donnerait ℂ[adhérence] = ℂ[orbite], d'où m_adhérence = m_orbite pour toutes les irreps.

**Obstruction :** **[PROUVÉ]** Le bord est de codimension 1 (Kumar, Cor. 5.6) et l'adhérence est non normale (Kumar, Thm 4.5). Ces deux faits **invalident simultanément** cette piste :
- La codimension 1 interdit Hartogs (condition nécessaire : codim ≥ 2).
- La non-normalité interdit en outre l'égalité ℂ[adhérence] = ℂ[orbite].

**La Piste A (Hartogs/normalité) est FERMÉE.** Ce n'est pas une conjecture réfutée — c'est un résultat de la littérature (Kumar 2013).

### 3.3. Piste B : généralisation du collapse lemma pour des familles d'irreps

**Piste B.1 — Réduction à la normalisation (Kumar) :**

Kumar établit (Theorem 7.5, Corollary 8.3 de [1]) une description complète de la GL(S₁)-décomposition de ℂ[normalisation(X_perm)] en termes de ℂ[normalisation(X_{perm_n})]. La décomposition de la normalisation est connue *modulo* la connaissance de ℂ[X_{perm_n}].

**Formule :**
$$m_{\mathrm{adhérence}}(\lambda) = \dim V_\lambda^{H_{\mathrm{perm}}} \leq m_{\mathrm{normalisation}}(\lambda) = \dim (V_\lambda \otimes \mathbb{C}[\tilde{X}_{\mathrm{perm}}])^{GL_9}$$

**[CONJECTURÉ — à tester]** La multiplicité dans la normalisation pourrait être calculable pour de petites irreps, ce qui bornerait la contribution du bord.

**[LACUNE]** Calculer m_normalisation(λ) pour les irreps serrées (notamment [2,2,1,1,1,1,1] où m_orbite = 1) nécessiterait la connaissance explicite de la normalisation de X_{perm₃}, qui n'est pas connue en général.

**Piste B.2 — Suite exacte de restriction au bord :**

Pour une composante C ⊂ ∂X_{perm} de codimension 1, les sections sur le bord contribuent potentiellement à m_adhérence − m_orbite. Mais sans classification des composantes C, cette suite ne peut pas être évaluée.

**[LACUNE]** Identifier les composantes du bord ∂(GL₉·perm₃) est un prérequis pour l'approche par suite exacte.

---

## 4. Reformulation précise du verrou restant

Pour conclure **absence d'obstruction pour l'adhérence**, il faut :

$$m_{\mathrm{adhérence}}(\mathrm{perm}, \lambda) \geq m_{\mathrm{adhérence}}(\det, \lambda) \quad \text{pour toutes les irreps } \lambda.$$

Or on sait :
- m_adhérence(perm, λ) ≥ m_orbite(perm, λ) — **trivial**
- m_orbite(perm, λ) ≥ m_orbite(det, λ) — **calculé, validé**
- m_adhérence(det, λ) ≥ m_orbite(det, λ) — **trivial**

Donc le problème se réduit à contrôler :

> **Est-ce que m_adhérence(det, λ) > m_orbite(det, λ) pour certaines λ ?**

Si cette différence est **nulle pour toutes λ** (ce qui serait vrai si X_{det} était normale, mais Kumar Thm 3.8 prouve qu'elle ne l'est pas pour m ≥ 3), les calculs existants suffisent. Sinon, il faut quantifier la différence.

**[LACUNE CRITIQUE]** La quantité m_det_adhérence(λ) − m_det_orbite(λ) pour les irreps de dim > 1 est inconnue. C'est le vrai verrou : sans contrôler cette différence, on ne peut pas conclure même avec tous les calculs d=3, d=4.

---

## 5. Plan d'action proposé (à valider par Sasha)

### Action 1 (prioritaire) : Calculer m_adhérence(det₃, λ) via les composantes du bord connu

**[FAISABLE avec les outils existants]**

Le bord de GL₉·det₃ a exactement deux composantes irréductibles C₁, C₂ (Hüttenhain–Lairez [2]). Pour chaque composante C_i :
- Identifier le polynôme f_i correspondant (de la description dans [2]).
- Calculer le stabilisateur H_{f_i}.
- Calculer dim V_λ^{H_{f_i}} pour les irreps serrées (notamment [2,2,1,1,1,1,1]).
- En déduire si m_adhérence(det, λ) > m_orbite(det, λ) pour ces irreps.

Si la différence est nulle pour toutes les irreps serrées, les calculs §2–§3 de la synthèse précédente deviennent rigoureux pour ces irreps.

### Action 2 : Identifier au moins une composante du bord de perm₃

Utiliser les scripts `degen_groebner.py`, `degen_curves.py`, `degen_fast.py` pour :
- Tester numériquement des courbes d'approche de perm₃ (dégénérescences cubiques, pas encore explorées).
- Identifier un polynôme limite f tel que f ∈ ∂(GL₉·perm₃).
- En calculer le stabilisateur.

*Rappel :* Les scripts existants ont déjà établi que les approches linéaires et quadratiques ne fonctionnent pas (coût numérique minimum ℝ ≈ 1.19, ℂ ≈ 0.90, jamais 0).

### Action 3 : Prouver la coïncidence pour une famille d'irreps via l'argument GIT

Si f ∈ ∂(GL₉·perm₃) est identifié et que son stabilisateur H_f contient H_{perm₃} comme sous-groupe, une invariante relative pour λ sur l'orbite se prolonge automatiquement au bord via l'argument de GIT. Cela donnerait m_adhérence(perm, λ) = m_orbite(perm, λ) pour les irreps telles que V_λ^{H_{perm₃}} → V_λ^{H_f} est surjectif.

---

## 6. Sous-produit potentiellement publiable : mesure de Weyl torsadée

**[CONJECTURÉ — non vérifié dans la littérature]** La formule

$$\mathrm{tr}(\tau \mid V_\lambda^K) = \frac{1}{n!}\,\mathrm{CT}_a\Big[ s_\lambda(\mathrm{eigs}_\tau)\cdot |\Delta(a^2)|^2 \Big]$$

(avec variables longues a_i² et système de racines plié A_{n-1}) est peut-être nouvelle comme calcul d'invariants tordus pour la paire symétrique (SL_n × SL_n, diag(SL_n)) avec automorphisme d'échange. Elle est validée numériquement sur toutes les partitions de 9 et de 16, avec les vérifications [1⁹] → tr=−1 et [9] → tr=+1.

À placer dans la littérature des paires symétriques / restriction de Weyl pliée (Casselman–Osborne, Vogan). **Demande à Sasha si cette piste vaut un effort de vérification bibliographique séparé.**

---

## 7. Résumé du statut

| Question | Statut | Source |
|---|---|---|
| X_{perm₃} non normale | **PROUVÉ** | Kumar [1], Thm 4.5 |
| Bord de codim 1 | **PROUVÉ** | Kumar [1], Cor. 5.6 |
| Hartogs/normalité inapplicable | **PROUVÉ** | Conséquence directe |
| Classification du bord de det₃ | **PROUVÉ** (2 composantes) | Hüttenhain–Lairez [2] |
| Classification du bord de perm₃ | **LACUNE OUVERTE** | — |
| m_orbite(perm) ≥ m_orbite(det) en d=3,4 | **CALCULÉ** | Scripts validés |
| m_adhérence(det, λ) = m_orbite(det, λ) ? | **LACUNE** | À calculer via [2] |
| m_adhérence(perm, λ) = m_orbite(perm, λ) ? | **LACUNE** | Prérequis : bord de perm |
| Absence d'obstruction pour l'adhérence | **OUVERT** | Dépend des 2 lacunes ci-dessus |
| Q31 (det₃ ∈ adhérence(GL₉·perm₃) ?) | **OUVERT** | — |

---

## Références

[1] S. Kumar, *Geometry of orbits of permanents and determinants*, Comment. Math. Helv. **88** (2013) 759–788. — https://kumar.math.unc.edu/papers/kumar62.pdf

[2] J. Hüttenhain, P. Lairez, *Les composantes du bord de la clôture de l'orbite du déterminant*, C. R. Acad. Sci. Paris **354** (2016) 931–935. — https://www.numdam.org/item/10.1016/j.crma.2016.07.002.pdf

[3] J. Hüttenhain, *Geometric Complexity Theory and Orbit Closures of Homogeneous Forms*, Doctoral Thesis, TU Berlin, 2017. — https://d-nb.info/1156010608/34

[4] P. Bürgisser, C. Ikenmeyer, J. Hüttenhain, *Permanent versus determinant: not via saturations*, Proc. AMS **145** (2017). — https://www3.math.tu-berlin.de/algebra/work/gct-sat-det-1.pdf
