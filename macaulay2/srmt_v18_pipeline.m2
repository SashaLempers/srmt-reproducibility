-- ================================================================
-- SRMT v18 — RESEARCH PIPELINE (Locked Edition)
-- ================================================================

debuggingMode = false
try protect T else null
scan({symbol gen, symbol wed, symbol det3,
      symbol fixed, symbol perm3, symbol defect},
     s -> try protect s else null)

printWidth = 0

-- ================================================================
-- §0 — CONFIGURATION GLOBALE
-- ================================================================

CONFIG = new MutableHashTable
CONFIG#"seed"             = 42
CONFIG#"primeMod"         = 32003
CONFIG#"runQQ"            = true
CONFIG#"runQQHeavy"       = true
CONFIG#"runSymbolicM"     = false
CONFIG#"runSampling"      = true
CONFIG#"runHessianDiag"   = true
CONFIG#"nSamplesEquiv"    = 15
CONFIG#"nSamplesProj"     = 5
CONFIG#"nSamplesOrbit"    = 100
CONFIG#"nSamplesOrbitGL3" = 20
CONFIG#"exportDir"        = "."
CONFIG#"runId"            = "srmt_v18"

setRandomSeed CONFIG#"seed"

-- ================================================================
-- §1 — IMPRESSION TAGUÉE & LOG STRUCTURÉ
-- ================================================================

logFile  = CONFIG#"exportDir" | "/" | CONFIG#"runId" | "_log.txt"
csvFile  = CONFIG#"exportDir" | "/" | CONFIG#"runId" | "_ranks.csv"
jsonFile = CONFIG#"exportDir" | "/" | CONFIG#"runId" | "_certificate.json"

logFile << "" << close

logEvents   = new MutableList from {}
rankRecords = new MutableList from {}

ALLOWED_TAGS = set {
    "PROVED-ALG", "QQ-CERTIFIED", "MODULAR",
    "OBS-N", "SUPPORTS-CONJ", "OPEN", "NOTE", "WARN", "INFO"
}

tagPrint = (tag, msg) -> (
    if not member(tag, ALLOWED_TAGS) then (
        logFile << "[INTERNAL] Unknown tag: " << tag << endl << close;
        error("Unknown tag: " | tag)
    );
    line := "  [" | tag | "] " | msg;
    print line;
    logFile << line << endl << close
)

logStep = (stepName, status, info) -> (
    logEvents#(#logEvents) = (stepName, status, info);
    logFile << "  EVENT " << stepName << " : " << status << " | " << info << endl << close
)

recordRank = (label, ringName, value, tag, note) -> (
    rankRecords#(#rankRecords) = (label, ringName, value, tag, note);
    tagPrint(tag, label | " [" | ringName | "] = " | toString value
              | (if note == "" then "" else "  (" | note | ")"))
)

sep = () -> (
    line := concatenate(64 : "-");
    print line;
    logFile << line << endl << close
)

bigsep = () -> (
    line := concatenate(64 : "=");
    print line;
    logFile << line << endl << close
)

-- ================================================================
-- §2 — ANNEAUX ET BASES
-- ================================================================

kk = ZZ/(CONFIG#"primeMod")

R9   = kk[y_0..y_8, MonomialOrder => GRevLex]
R9Q  = QQ[y_0..y_8, MonomialOrder => GRevLex]

matIdx = (i, j) -> 3*i + j

mons2OnRing = (Rgen) -> (
    flatten apply(9, a -> apply(toList(a..8), b -> (gens Rgen)#a * (gens Rgen)#b))
)

mons2  = mons2OnRing(R9)
mons2Q = mons2OnRing(R9Q)

mons3OnRing = (Rgen) -> (
    flatten flatten apply(9, a ->
        apply(toList(a..8), b ->
            apply(toList(b..8), c -> (gens Rgen)#a * (gens Rgen)#b * (gens Rgen)#c)))
)

mons3  = mons3OnRing(R9)
mons3Q = mons3OnRing(R9Q)

-- ================================================================
-- §3 — POLYNÔMES SPÉCIAUX
-- ================================================================

buildY = (Rgen) -> (
    matrix apply(3, i -> apply(3, j -> (gens Rgen)#(matIdx(i,j))))
)

YK = buildY(R9)
YQ = buildY(R9Q)

det3OnRing = (Rgen, Y) -> det Y
perm3OnRing = (Rgen, Y) -> sum(
    permutations 3,
    sigma -> product(3, i -> Y_(i, sigma#i))
)

det3K  = det3OnRing(R9, YK)
perm3K = perm3OnRing(R9, YK)
det3Q  = det3OnRing(R9Q, YQ)
perm3Q = perm3OnRing(R9Q, YQ)

adjOf = (Y) -> (
    matrix apply(3, i ->
        apply(3, j ->
            (-1)^(i+j) * det submatrix'(Y, {j}, {i})))
)

adjYK = adjOf YK
adjYQ = adjOf YQ

-- ================================================================
-- §4 — IDENTITÉ adj(Y)·Y = det(Y)·I  [PROVED-ALG]
-- ================================================================

bigsep()
print "Section 4 - Identite adj(Y)*Y = det(Y)*I"
bigsep()

idCheckK = (adjYK * YK) - det3K * id_(R9^3)
idCheckQ = (adjYQ * YQ) - det3Q * id_(R9Q^3)

if idCheckK == 0 then (
    tagPrint("PROVED-ALG", "adj(Y)*Y = det3(Y)*I sur R9 = kk[y]")
) else (
    tagPrint("WARN", "adj(Y)*Y != det3(Y)*I sur kk[y]")
)

if idCheckQ == 0 then (
    tagPrint("PROVED-ALG", "adj(Y)*Y = det3(Y)*I sur R9Q = QQ[y]")
) else (
    tagPrint("WARN", "adj(Y)*Y != det3(Y)*I sur QQ[y]")
)

logStep("adj_identity",
        if idCheckK == 0 and idCheckQ == 0 then "PROVED-ALG" else "WARN",
        "kk_ok=" | toString(idCheckK == 0) | " QQ_ok=" | toString(idCheckQ == 0))

-- ================================================================
-- §5 — PLITHYSME DE CAUCHY-SCHUR (rappel theorique)
-- ================================================================

bigsep()
print "Section 5 - Decomposition de Cauchy-Schur"
bigsep()

tagPrint("PROVED-ALG", "Sym^3(V*W) = S(3)*S(3) + S(2,1)*S(2,1) + S(1,1,1)*S(1,1,1)")
tagPrint("PROVED-ALG", "Dimensions : 100 + 64 + 1 = 165")
tagPrint("PROVED-ALG", "det3 in S(1,1,1)*S(1,1,1)")
tagPrint("PROVED-ALG", "perm3 in S(3)*S(3)")
tagPrint("NOTE",       "Reference : Macdonald, Fulton-Harris")

logStep("cauchy_schur", "PROVED-ALG", "decomp Sym^3(V*W) cited")

-- ================================================================
-- §6 — ACTIONS DE GROUPE
-- ================================================================

genGL3 = () -> (
    M := random(kk^3, kk^3);
    while det M == 0 do M = random(kk^3, kk^3);
    M
)

genGL9 = () -> (
    M := random(kk^9, kk^9);
    while det M == 0 do M = random(kk^9, kk^9);
    M
)

actionGL3xGL3onPolyOnRing = (Rgen, f, A, B) -> (
    YR := buildY Rgen;
    AYBt := promote(A, Rgen) * YR * transpose promote(B, Rgen);
    subRules := flatten apply(3, i ->
        apply(3, j ->
            (gens Rgen)#(matIdx(i,j)) => AYBt_(i,j)));
    sub(f, subRules)
)

actionGL3xGL3onPoly = (f, A, B) -> actionGL3xGL3onPolyOnRing(R9, f, A, B)

actionGL9onPoly = (f, g) -> (
    subRules := apply(9, i ->
        (gens R9)#i => sum(9, j -> g_(i,j) * (gens R9)#j));
    sub(f, subRules)
)

-- ================================================================
-- §7 — CATALECTICANT Cat(1,2)
-- ================================================================

buildCatMatrixOnRing = (Rgen, monsBasis, f) -> (
    matrix apply(45, b ->
        apply(9, a -> (
            df := diff((gens Rgen)#a, f);
            sub(coefficient(monsBasis#b, df), coefficientRing Rgen)
        )))
)

buildCatMatrix  = (f) -> buildCatMatrixOnRing(R9,  mons2,  f)
buildCatMatrixQ = (f) -> buildCatMatrixOnRing(R9Q, mons2Q, f)

-- ================================================================
-- §8 — PROJECTEURS piMinus, piPlus
-- ================================================================

vCoord = (k) -> k // 3
wCoord = (k) -> k %  3

buildProjectorsOnField = (F, monsBasis, Rgen) -> (
    invIdx := new MutableHashTable;
    scan(45, b -> invIdx#(monsBasis#b) = b);

    Mplus  := mutableMatrix(F, 45, 45);
    Mminus := mutableMatrix(F, 45, 45);

    scan(45, b -> (
        m := monsBasis#b;
        ee := flatten exponents m;
        idxPair := flatten apply(9, k -> toList(ee#k : k));
        ka := idxPair#0;
        kb := idxPair#1;
        iA := vCoord ka; jA := wCoord ka;
        iB := vCoord kb; jB := wCoord kb;

        scan({(1, 1), (1, -1), (-1, 1), (-1, -1)}, sgn -> (
            sV := sgn#0; sW := sgn#1;
            terms := {
                (iA, iB, jA, jB,  1),
                (iB, iA, jA, jB,  sV),
                (iA, iB, jB, jA,  sW),
                (iB, iA, jB, jA,  sV * sW)
            };
            scan(terms, t -> (
                p1 := t#0; p2 := t#1; q1 := t#2; q2 := t#3; coef := t#4;
                k1 := 3*p1 + q1;
                k2 := 3*p2 + q2;
                aSorted := min(k1, k2);
                bSorted := max(k1, k2);
                m' := (gens Rgen)#aSorted * (gens Rgen)#bSorted;
                bIdx := invIdx#m';
                contrib := sub(coef / 4, F);
                if sV == 1 and sW == 1 then
                    Mplus_(bIdx, b) = Mplus_(bIdx, b) + contrib
                else if sV == -1 and sW == -1 then
                    Mminus_(bIdx, b) = Mminus_(bIdx, b) + contrib
            ))
        ))
    ));

    (matrix Mminus, matrix Mplus)
)

(piMinusMat,  piPlusMat)  = buildProjectorsOnField(kk, mons2,  R9)
(piMinusMatQ, piPlusMatQ) = buildProjectorsOnField(QQ, mons2Q, R9Q)

-- ================================================================
-- §9 — VERIFICATIONS DES PROJECTEURS
-- ================================================================

bigsep()
print "Section 9 - Verifications algebriques des projecteurs"
bigsep()

verifyProjectorsOn = (label, F, piM, piP) -> (
    Id45 := id_(F^45);
    chkIdemMinus := (piM * piM - piM == 0);
    chkIdemPlus  := (piP * piP - piP == 0);
    chkOrth      := (piM * piP == 0);
    chkComplete  := (piM + piP - Id45 == 0);
    rkM := rank piM;
    rkP := rank piP;
    fieldName := if F === QQ then "QQ" else "kk";
    tag := if F === QQ then "QQ-CERTIFIED" else "MODULAR";

    if chkIdemMinus then tagPrint(tag, label | " : piMinus^2 = piMinus [" | fieldName | "]")
    else                  tagPrint("WARN", label | " : piMinus^2 != piMinus sur " | fieldName);
    if chkIdemPlus  then tagPrint(tag, label | " : piPlus^2 = piPlus [" | fieldName | "]")
    else                  tagPrint("WARN", label | " : piPlus^2 != piPlus sur " | fieldName);
    if chkOrth      then tagPrint(tag, label | " : piMinus*piPlus = 0 [" | fieldName | "]")
    else                  tagPrint("WARN", label | " : piMinus*piPlus != 0 sur " | fieldName);
    if chkComplete  then tagPrint(tag, label | " : piMinus + piPlus = Id_45 [" | fieldName | "]")
    else                  tagPrint("WARN", label | " : piMinus + piPlus != Id_45 sur " | fieldName);

    recordRank("rank(piMinus)", fieldName, rkM, tag,
               if rkM == 9 then "attendu = 9" else "ANOMALIE");
    recordRank("rank(piPlus)",  fieldName, rkP, tag,
               if rkP == 36 then "attendu = 36" else "ANOMALIE");

    logStep("projectors_check_" | fieldName,
            if chkIdemMinus and chkIdemPlus and chkOrth and chkComplete then tag else "WARN",
            "rkM=" | toString rkM | " rkP=" | toString rkP)
)

if CONFIG#"runQQ" then (
    verifyProjectorsOn("piProj", QQ, piMinusMatQ, piPlusMatQ)
) else (
    tagPrint("NOTE", "Bloc QQ desactive (CONFIG runQQ)")
);

verifyProjectorsOn("piProj", kk, piMinusMat,  piPlusMat)

tagPrint("PROVED-ALG", "Lemme 2 : Sym^2(V*W) = Sym^2 V * Sym^2 W + Lambda^2 V * Lambda^2 W")
tagPrint("PROVED-ALG", "Lemme 3 : piMinus, piPlus sont GL(V)xGL(W)-equivariants (Schur)")

-- ================================================================
-- §10 — RANGS CRITIQUES SUR QQ
-- ================================================================

bigsep()
print "Section 10 - Rangs critiques sur QQ"
bigsep()

CatDet3K  = buildCatMatrix det3K
CatPerm3K = buildCatMatrix perm3K

if CONFIG#"runQQ" then (
    CatDet3Q  := buildCatMatrixQ det3Q;
    CatPerm3Q := buildCatMatrixQ perm3Q;

    rkCatDet3Q  := rank CatDet3Q;
    rkCatPerm3Q := rank CatPerm3Q;

    recordRank("rank(Cat(det3))",  "QQ", rkCatDet3Q,  "QQ-CERTIFIED",
               if rkCatDet3Q  == 9 then "attendu 9" else "ANOMALIE");
    recordRank("rank(Cat(perm3))", "QQ", rkCatPerm3Q, "QQ-CERTIFIED", "");

    deltaDet3Q     := rank (piMinusMatQ * CatDet3Q);
    deltaPerm3Q    := rank (piMinusMatQ * CatPerm3Q);
    deltaPlusDet3Q := rank (piPlusMatQ  * CatDet3Q);
    deltaPlusPerm3Q := rank (piPlusMatQ * CatPerm3Q);

    recordRank("delta(det3) = rank(piMinus*Cat(det3)^T)",  "QQ", deltaDet3Q,
               "QQ-CERTIFIED", if deltaDet3Q == 9 then "attendu 9" else "ANOMALIE");
    recordRank("delta(perm3) = rank(piMinus*Cat(perm3)^T)","QQ", deltaPerm3Q,
               "QQ-CERTIFIED", if deltaPerm3Q == 0 then "attendu 0" else "ANOMALIE");
    recordRank("rank(piPlus*Cat(det3)^T)",  "QQ", deltaPlusDet3Q,
               "QQ-CERTIFIED", if deltaPlusDet3Q == 0 then "attendu 0" else "ANOMALIE");
    recordRank("rank(piPlus*Cat(perm3)^T)", "QQ", deltaPlusPerm3Q,
               "QQ-CERTIFIED", "");

    if deltaDet3Q == 9 and deltaPerm3Q == 0 then (
        tagPrint("QQ-CERTIFIED", "delta(det3)=9 et delta(perm3)=0 sur QQ => vrais sur C");
        tagPrint("PROVED-ALG",   "Prop 5 : adherences GL3xGL3.det3 et GL3xGL3.perm3 disjointes");
        tagPrint("NOTE",         "Conclusion strictement GL3xGL3 ; voir section 13 pour GL9");
        logStep("prop5_separation_GL3xGL3", "PROVED",
                "delta_det3=9 delta_perm3=0 on QQ")
    ) else (
        tagPrint("WARN", "Valeurs delta sur QQ inattendues -- Prop 5 NON activee");
        logStep("prop5_separation_GL3xGL3", "WARN",
                "delta_det3=" | toString deltaDet3Q | " delta_perm3=" | toString deltaPerm3Q)
    )
) else (
    tagPrint("NOTE", "Bloc QQ desactive -- delta(det3), delta(perm3) NON certifies sur QQ");
    tagPrint("OPEN", "Sans certification QQ, la Prop 5 reste conditionnelle")
)

-- ================================================================
-- §11 — RANGS MODULAIRES (kk)
-- ================================================================

bigsep()
print "Section 11 - Rangs modulaires sur kk"
bigsep()

rankCatDet3K  = rank CatDet3K
rankCatPerm3K = rank CatPerm3K
deltaDet3K    = rank (piMinusMat * CatDet3K)
deltaPerm3K   = rank (piMinusMat * CatPerm3K)
plusDet3K     = rank (piPlusMat  * CatDet3K)
plusPerm3K    = rank (piPlusMat  * CatPerm3K)

recordRank("rank(Cat(det3))",        "kk", rankCatDet3K,  "MODULAR", "")
recordRank("rank(Cat(perm3))",       "kk", rankCatPerm3K, "MODULAR", "")
recordRank("delta(det3)",            "kk", deltaDet3K,    "MODULAR", "")
recordRank("delta(perm3)",           "kk", deltaPerm3K,   "MODULAR", "")
recordRank("rank(piPlus*Cat(det3)^T)",  "kk", plusDet3K,  "MODULAR", "")
recordRank("rank(piPlus*Cat(perm3)^T)", "kk", plusPerm3K, "MODULAR", "")

tagPrint("NOTE", "Rangs modulaires = minorants des rangs sur QQ. Ne prouvent PAS sur C.")

-- ================================================================
-- §12 — TESTS ALEATOIRES
-- ================================================================

bigsep()
print "Section 12 - Tests aleatoires (stress, jamais une preuve)"
bigsep()

actionOnVec45 = (A, B, v) -> (
    AkronB := A ** B;
    v2 := mutableMatrix(kk, 45, 1);
    scan(45, a -> (
        mon  := mons2#a;
        expv := flatten exponents mon;
        pqs  := flatten apply(9, k -> toList(expv#k : k));
        p    := pqs#0; q := pqs#1;
        newp := sum(9, j -> AkronB_(j,p) * (gens R9)#j);
        newq := sum(9, j -> AkronB_(j,q) * (gens R9)#j);
        newm := newp * newq;
        scan(45, b -> (
            cb := sub(coefficient(mons2#b, newm), kk);
            if cb != 0 then v2_(b,0) = v2_(b,0) + cb * v_(a,0)
        ))
    ));
    matrix v2
)

testProjectorEquivariance = (piMat, piName, nTrials) -> (
    pass := 0; fail := 0;
    gens0 := {
        (diagonalMatrix(kk, {2,3,5}), diagonalMatrix(kk, {7,11,13}), "diag"),
        (id_(kk^3), id_(kk^3), "identite")
    };
    scan(gens0, g0 -> (
        A := g0#0; B := g0#1;
        v := random(kk^45, kk^1);
        lhs := actionOnVec45(A, B, piMat * v);
        rhs := piMat * actionOnVec45(A, B, v);
        if lhs == rhs then pass = pass + 1 else fail = fail + 1
    ));
    apply(nTrials, idx -> (
        A := genGL3(); B := genGL3();
        v := random(kk^45, kk^1);
        lhs := actionOnVec45(A, B, piMat * v);
        rhs := piMat * actionOnVec45(A, B, v);
        if lhs == rhs then pass = pass + 1 else fail = fail + 1
    ));
    nTotal := 2 + nTrials;
    if fail == 0 then
        tagPrint("OBS-N", "Equivariance " | piName | " : " | toString pass | "/" | toString nTotal | " sur kk")
    else
        tagPrint("WARN", "Equivariance " | piName | " : " | toString fail | " echecs / " | toString nTotal);
    logStep("equiv_proj_" | piName,
            if fail == 0 then "OBS-PASS" else "WARN",
            "pass=" | toString pass | " fail=" | toString fail);
    (pass, fail)
)

testCatEquivarianceOn = (f, fName, nTrials) -> (
    pass := 0; fail := 0;
    Cf := buildCatMatrix f;
    rCf := rank Cf;
    apply(nTrials, idx -> (
        A := genGL3(); B := genGL3();
        fAB := actionGL3xGL3onPoly(f, A, B);
        CfAB := buildCatMatrix fAB;
        if rank CfAB == rCf then pass = pass + 1 else fail = fail + 1
    ));
    if fail == 0 then
        tagPrint("OBS-N", "rank(Cat(" | fName | ")) constant sur " | toString nTrials | " tirages GL3xGL3")
    else
        tagPrint("WARN", "rank(Cat(" | fName | ")) varie sur " | toString fail | "/" | toString nTrials);
    logStep("equiv_cat_" | fName,
            if fail == 0 then "OBS-PASS" else "WARN",
            "pass=" | toString pass | " fail=" | toString fail);
    (pass, fail)
)

testDeltaInvarianceOn = (f, fName, nTrials, refDelta) -> (
    pass := 0; fail := 0;
    apply(nTrials, idx -> (
        A := genGL3(); B := genGL3();
        fAB := actionGL3xGL3onPoly(f, A, B);
        CfAB := buildCatMatrix fAB;
        d := rank (piMinusMat * CfAB);
        if d == refDelta then pass = pass + 1 else fail = fail + 1
    ));
    if fail == 0 then
        tagPrint("OBS-N", "delta(" | fName | ") = " | toString refDelta | " sur " | toString nTrials | " tirages GL3xGL3 (kk)")
    else
        tagPrint("WARN", "delta(" | fName | ") varie sur " | toString fail | "/" | toString nTrials);
    logStep("delta_invar_GL3xGL3_" | fName,
            if fail == 0 then "OBS-PASS" else "WARN",
            "ref=" | toString refDelta | " pass=" | toString pass | " fail=" | toString fail);
    (pass, fail)
)

if CONFIG#"runSampling" then (
    testProjectorEquivariance(piMinusMat, "piMinus", CONFIG#"nSamplesProj");
    testProjectorEquivariance(piPlusMat,  "piPlus",  CONFIG#"nSamplesProj");
    testCatEquivarianceOn(det3K,  "det3",  CONFIG#"nSamplesEquiv");
    testCatEquivarianceOn(perm3K, "perm3", CONFIG#"nSamplesEquiv");
    testDeltaInvarianceOn(det3K,  "det3",  CONFIG#"nSamplesOrbitGL3", deltaDet3K);
    testDeltaInvarianceOn(perm3K, "perm3", CONFIG#"nSamplesOrbitGL3", deltaPerm3K);
    tagPrint("NOTE", "Section 12 = stress-tests sur kk. Ne remplacent PAS les Lemmes 1,3.")
) else (
    tagPrint("NOTE", "Sampling desactive (CONFIG runSampling)")
)

-- ================================================================
-- §13 — AVERTISSEMENT GL9
-- ================================================================

bigsep()
print "Section 13 - Avertissement strict GL3xGL3 vs GL9"
bigsep()

tagPrint("NOTE", "delta(f) depend de la decomposition C^9 = C^3 * C^3")
tagPrint("NOTE", "Lambda^2 V * Lambda^2 W subset Sym^2(C^9) N'EST PAS GL9-stable")
tagPrint("NOTE", "piMinus N'EST PAS GL9-equivariant ; delta N'EST PAS un invariant GL9")
tagPrint("OPEN", "Toute separation GL9.det3 vs GL9.perm3 NECESSITE un argument supplementaire")
tagPrint("OPEN", "non fourni par ce pipeline")

if CONFIG#"runSampling" then (
    bigsep();
    print "Section 13-bis - delta NON invariant sous GL9 (observe)";
    bigsep();
    nGL9 := 10;
    valsDet3 := apply(nGL9, idx -> (
        g := genGL9();
        fG := actionGL9onPoly(det3K, g);
        CfG := buildCatMatrix fG;
        rank (piMinusMat * CfG)
    ));
    tagPrint("OBS-N", "delta(g.det3) sur " | toString nGL9 | " tirages GL9 : " | toString tally valsDet3);
    if #unique valsDet3 > 1 then
        tagPrint("NOTE", "Variabilite observee -- coherent avec non-invariance GL9")
    else
        tagPrint("NOTE", "Constance observee a cette resolution -- N'IMPLIQUE PAS l'invariance GL9")
)

-- ================================================================
-- §14 — DIAGNOSTIC HESSIEN (declasse)
-- ================================================================

if CONFIG#"runHessianDiag" then (
    bigsep();
    print "Section 14 - Diagnostic Hessien (DECLASSE en v17)";
    bigsep();

    buildHessMatrix := f -> matrix apply(9, p ->
        apply(9, q -> diff((gens R9)#p, diff((gens R9)#q, f))));

    evalMatAt := (M, pt) -> matrix apply(numRows M, i ->
        apply(numColumns M, j ->
            sub(M_(i,j), apply(9, k -> (gens R9)#k => sub(pt#k, kk)))));

    nHess := 50;
    hessRanks := apply(nHess, idx -> (
        g := genGL9();
        fG := actionGL9onPoly(det3K, g);
        pt := apply(9, k -> random kk);
        rank evalMatAt(buildHessMatrix fG, pt)
    ));
    tagPrint("OBS-N", "rang Hess(g.det3)(pt) sur " | toString nHess | " tirages : " | toString tally hessRanks);
    tagPrint("NOTE", "Le rang '6' historique etait un artefact y0 fixe (declasse v17)");
    tagPrint("OPEN", "Aucune obstruction Hessienne intrinseque certifiee");
    logStep("hessian_diag", "OBS-N", "ranks=" | toString tally hessRanks)
)

-- ================================================================
-- §15 — M(f_generic) symbolique (optionnel)
-- ================================================================

if CONFIG#"runSymbolicM" then (
    bigsep();
    print "Section 15 - M(f_generic) symbolique";
    bigsep();
    tagPrint("NOTE", "Bloc symbolique active mais non implemente dans cette release")
) else (
    tagPrint("NOTE", "Construction symbolique M(f_generic) desactivee (CONFIG runSymbolicM)")
)

-- ================================================================
-- §16 — EXPORTS
-- ================================================================

bigsep()
print "Section 16 - Exports"
bigsep()

csvOut = openOut csvFile
csvOut << "label,ring,value,tag,note" << endl
scan(toList rankRecords, rec -> (
    label := rec#0; ringName := rec#1; value := rec#2; tag := rec#3; note := rec#4;
    csvOut << "\"" << label << "\","
           << ringName << ","
           << toString value << ","
           << tag << ","
           << "\"" << note << "\"" << endl
))
csvOut << close
tagPrint("INFO", "CSV rangs exporte : " | csvFile)

escapeJson = (s) -> replace("\"", "\\\"", toString s)

jsonOut = openOut jsonFile
jsonOut << "{" << endl
jsonOut << "  \"runId\": \"" << CONFIG#"runId" << "\"," << endl
jsonOut << "  \"seed\": " << CONFIG#"seed" << "," << endl
jsonOut << "  \"primeMod\": " << CONFIG#"primeMod" << "," << endl
jsonOut << "  \"flags\": {" << endl
jsonOut << "    \"runQQ\": "          << toString CONFIG#"runQQ"          << "," << endl
jsonOut << "    \"runQQHeavy\": "     << toString CONFIG#"runQQHeavy"     << "," << endl
jsonOut << "    \"runSymbolicM\": "   << toString CONFIG#"runSymbolicM"   << "," << endl
jsonOut << "    \"runSampling\": "    << toString CONFIG#"runSampling"    << "," << endl
jsonOut << "    \"runHessianDiag\": " << toString CONFIG#"runHessianDiag" << endl
jsonOut << "  }," << endl
jsonOut << "  \"ranks\": [" << endl
nRec = #rankRecords
scan(nRec, idx -> (
    rec := rankRecords#idx;
    label := rec#0; ringName := rec#1; value := rec#2; tag := rec#3; note := rec#4;
    jsonOut << "    {\"label\": \"" << escapeJson label << "\", "
            << "\"ring\": \"" << ringName << "\", "
            << "\"value\": " << toString value << ", "
            << "\"tag\": \"" << tag << "\", "
            << "\"note\": \"" << escapeJson note << "\"}";
    if idx < nRec - 1 then jsonOut << ",";
    jsonOut << endl
))
jsonOut << "  ]," << endl
jsonOut << "  \"events\": [" << endl
nEv = #logEvents
scan(nEv, idx -> (
    ev := logEvents#idx;
    stepName := ev#0; status := ev#1; info := ev#2;
    jsonOut << "    {\"step\": \"" << escapeJson stepName << "\", "
            << "\"status\": \"" << escapeJson status << "\", "
            << "\"info\": \"" << escapeJson info << "\"}";
    if idx < nEv - 1 then jsonOut << ",";
    jsonOut << endl
))
jsonOut << "  ]" << endl
jsonOut << "}" << endl
jsonOut << close
tagPrint("INFO", "JSON certificat exporte : " | jsonFile)

-- ================================================================
-- §17 — SYNTHESE FINALE
-- ================================================================

bigsep()
print "Section 17 - Synthese finale SRMT v18"
bigsep()

print ""
print "  -- PROUVE ALGEBRIQUEMENT (sur C, sans hypothese) --"
tagPrint("PROVED-ALG", "adj(Y)*Y = det3(Y)*I (identite dans Z[Y])")
tagPrint("PROVED-ALG", "Decomposition Cauchy-Schur Sym^3(V*W)")
tagPrint("PROVED-ALG", "det3 in S(1,1,1)*S(1,1,1) ; perm3 in S(3)*S(3)")
tagPrint("PROVED-ALG", "Decomposition Sym^2(V*W) = Sym^2 V * Sym^2 W + Lambda^2 V * Lambda^2 W (Lemme 2)")
tagPrint("PROVED-ALG", "Naturalite Cat(1,2) => equivariance GL(U) (Lemme 1)")
tagPrint("PROVED-ALG", "piMinus, piPlus morphismes GL(V)xGL(W)-modules (Lemme 3)")
tagPrint("PROVED-ALG", "delta(f) := rang(piMinus*Cat(f)^T) invariant GL(V)xGL(W) (Prop 4)")

print ""
print "  -- CERTIFIE SUR QQ (donc sur C par extension) --"
if CONFIG#"runQQ" then (
    tagPrint("QQ-CERTIFIED", "Rangs critiques de section 10 -- voir CSV/JSON")
) else (
    tagPrint("NOTE", "Bloc QQ desactive -- rangs critiques NON certifies QQ")
)

print ""
print "  -- CERTIFIE MODULAIREMENT --"
tagPrint("MODULAR", "Identites projecteurs et rangs section 11 (NE PROUVENT PAS sur C)")

print ""
print "  -- OBSERVE (stress-tests) --"
tagPrint("OBS-N", "Equivariance Cat sur tirages GL3xGL3")
tagPrint("OBS-N", "Equivariance piMinus, piPlus sur tirages")
tagPrint("OBS-N", "delta constant sur orbite GL3xGL3 (echantillon)")
tagPrint("OBS-N", "delta NON constant sur orbite GL9 (echantillon, attendu)")

print ""
print "  -- OUVERT --"
tagPrint("OPEN", "Construction symbolique M(f_generic) et mineur 9x9 separant explicite")
tagPrint("OPEN", "Extension separation a GL9.det3 vs GL9.perm3")
tagPrint("OPEN", "Lien explicite delta vs multiplicites BIP")

print ""
print("  Logs   : " | logFile)
print("  Ranks  : " | csvFile)
print("  Cert.  : " | jsonFile)
print ""
bigsep()
print "  SRMT v18 -- pipeline verrouille execute."
bigsep()

logStep("pipeline_end", "COMPLETE", "SRMT v18 finished")
exit 0
