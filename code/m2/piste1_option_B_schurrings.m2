-- Option B : multiplicite de l'invariant SL9 via plethysme (SchurRings)
-- Objectif : pour chaque degre d, calculer dim (Sym^d (Sym^3 C^9))^{SL_9}
-- Methode : plethysme s_d[s_3] et extraction du coeff de s_{(k,k,...,k)}
-- Pas de matrice, pas de GC qui sature.

needsPackage "SchurRings"

print "=== OPTION B : invariants SL9 via SchurRings ===";

n = 9;
-- Anneau de Schur sur 9 variables (representations polynomiales de GL_9)
R = schurRing(QQ, s, n)

-- Representation V = Sym^3 C^9 = s_{(3,0,...,0)}
V = s_{3}

-- Fonction : dim des SL_n-invariants dans Sym^d V
-- = somme des multiplicites de s_{(k,k,...,k)} dans Sym^d V, pour tout k >= 0
-- (SL_n-invariants = GL_n-rep. avec tous les poids egaux = multiples du determinant)
dimSLninvariants = (d) -> (
    SymdV := plethysm(s_{d}, V);
    L := listForm SymdV;   -- liste de paires (partition, coef)
    total := 0;
    for pair in L do (
        lam := pair#0;
        coef := pair#1;
        -- on garde les partitions "rectangulaires" (k,k,...,k) a n parts egales
        -- en SchurRings, une partition est stockee tronquee (sans les zeros),
        -- donc (k,k,...,k) avec k>0 apparait comme une liste de n copies de k,
        -- et (0,0,...,0) apparait comme la liste vide {}.
        isRect := (
            (#lam == 0) or
            (#lam == n and all(lam, x -> x == lam#0))
        );
        if isRect then (
            print("  d=" | d | " : mult de " | toString lam | " = " | toString coef);
            total = total + coef
        )
    );
    total
);

print "";
print "--- Calcul des dimensions d'invariants SL9 pour d = 4, 6, 8, 12 ---";
for d in {3, 6, 9, 12, 15} do (
    dimInv := dimSLninvariants(d);
    print("Degre " | d | " : dim (Sym^" | d | " Sym^3 C^9)^{SL_9} = " | toString dimInv);
    print "";
);

print "=== FIN OPTION B ===";
exit 0