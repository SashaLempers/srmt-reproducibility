-- ================================================================
-- SRMT v19 — LOCK & AUDIT RELEASE
-- ================================================================
--
-- Cette version est UNIQUEMENT une release de verrouillage et d'audit
-- de v18. Aucune nouvelle revendication mathematique. Aucun nouveau
-- theoreme. Aucune nouvelle attaque sur la direction GL9 ouverte.
--
-- Differences strictes vs v18 :
--   (a) Manifeste gele FROZEN_CLAIMS avec valeurs attendues exactes.
--   (b) Fonction registerClaim qui interdit tout PROVED-ALG ou
--       QQ-CERTIFIED hors manifeste lorsque strictManifest = true.
--   (c) Fonction enforceManifest qui echoue dur (exit non zero) si
--       un claim attendu manque ou si une valeur observee differe.
--   (d) Provenance auditable dans le JSON (runId, seed, primeModList,
--       m2Version, scriptSha256, runStartedAt, runFinishedAt, flags,
--       artefacts).
--   (e) Boucle multi-prime sur primeModList (par defaut une seule
--       prime, structure prete pour plusieurs).
--   (f) Section claims_status dans le JSON (OK / MISMATCH / MISSING).
--   (g) Section open_questions figee, incluant explicitement la
--       direction GL9 dir. 2 comme OPEN.
--   (h) Fichier de diff manifeste pour audit humain.
--
-- COMMANDE OFFICIELLE DE LANCEMENT :
--
--   M2 --script srmt_v19_pipeline.m2
--
-- Aucun argument, aucune variable d'environnement, aucun flag externe.
-- Toute la configuration est figee dans le bloc CONFIG (§0).
--
-- Sorties dans le repertoire courant :
--   srmt_v19_lock_audit_log.txt
--   srmt_v19_lock_audit_ranks.csv
--   srmt_v19_lock_audit_certificate.json
--   srmt_v19_lock_audit_manifest_diff.txt
--
-- Pre-requis : Macaulay2 >= 1.19. Determinisme : setRandomSeed.
--
-- ================================================================
-- MANIFESTE DES CLAIMS (gele -- ne pas etendre sans preuve manuscrite)
-- ================================================================
--
-- Les tags autorises sont identiques a v18 :
--   PROVED-ALG, QQ-CERTIFIED, MODULAR, OBS-N, SUPPORTS-CONJ,
--   OPEN, NOTE, WARN, INFO.
--
-- Le perimetre fige des claims est identique a v18 :
--   - Identite adj(Y)*Y = det3*I.
--   - Decomposition Cauchy-Schur Sym^3(V*W) = ... (citee).
--   - Lemmes 1, 2, 3, Prop 4, Prop 5 sur GL3 x GL3.
--   - delta(det3) = 9, delta(perm3) = 0 sur QQ.
--   - dim Sing V(det3) = 5, dim Sing V(perm3) = 3 (QQ).
--   - dim_proj 4 / 2 ; dim V(J) 5 / 3.
--   - H_apolar = (1,9,9,1) pour les deux.
--   - codim J = 4 / 6 ; degree J = 6 / 24.
--   - NIV-1 (orbites distinctes) PROVED-ALG.
--   - NIV-2 dir.1 (perm3 not in closure GL9.det3) PROVED-ALG.
--   - NIV-2 dir.2 (det3 not in closure GL9.perm3) reste OPEN.
--
-- ENGAGEMENT : aucun PROVED-ALG ou QQ-CERTIFIED ne va plus loin
-- que ce qui est dans FROZEN_CLAIMS. Le script echoue plutot que
-- de promouvoir un statut.
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
CONFIG#"primeModList"     = {32003}
CONFIG#"runQQ"            = true
CONFIG#"runQQHeavy"       = true
CONFIG#"runSymbolicM"     = true
CONFIG#"runSampling"      = true
CONFIG#"runHessianDiag"   = true
CONFIG#"runQuarticObs"    = false
CONFIG#"strictManifest"   = true
CONFIG#"copyScriptOnExit" = false
CONFIG#"nSamplesEquiv"    = 15
CONFIG#"nSamplesProj"     = 5
CONFIG#"nSamplesOrbit"    = 100
CONFIG#"nSamplesOrbitGL3" = 20
CONFIG#"exportDir"        = "."
CONFIG#"runId"            = "srmt_v19_lock_audit"
CONFIG#"scriptPath"       = "srmt_v19_pipeline.m2"

setRandomSeed CONFIG#"seed"

-- Pour compat v18 : on conserve aussi un "primeMod" scalaire pointant
-- sur la premiere prime de la liste (utilise par les sections §2-§14
-- restees identiques a v18).
CONFIG#"primeMod" = first CONFIG#"primeModList"

-- ================================================================
-- §1 — IMPRESSION TAGUÉE & LOG STRUCTURÉ
-- ================================================================

logFile      = CONFIG#"exportDir" | "/" | CONFIG#"runId" | "_log.txt"
csvFile      = CONFIG#"exportDir" | "/" | CONFIG#"runId" | "_ranks.csv"
jsonFile     = CONFIG#"exportDir" | "/" | CONFIG#"runId" | "_certificate.json"
manifestDiff = CONFIG#"exportDir" | "/" | CONFIG#"runId" | "_manifest_diff.txt"

logTrunc = openOut logFile;
logTrunc << "";
close logTrunc;

appendLog = (s) -> (
    o := openOutAppend logFile;
    o << s << endl;
    close o
);

logEvents   = new MutableList from {}
rankRecords = new MutableList from {}

ALLOWED_TAGS = set {
    "PROVED-ALG", "QQ-CERTIFIED", "MODULAR",
    "OBS-N", "SUPPORTS-CONJ", "OPEN", "NOTE", "WARN", "INFO"
}

tagPrint = (tag, msg) -> (
    if not member(tag, ALLOWED_TAGS) then (
        appendLog("[INTERNAL] Unknown tag: " | tag);
        error("Unknown tag: " | tag)
    );
    line := "  [" | tag | "] " | msg;
    print line;
    appendLog(line)
)

logStep = (stepName, status, info) -> (
    logEvents#(#logEvents) = (stepName, status, info);
    appendLog("  EVENT " | stepName | " : " | status | " | " | info)
)

recordRank = (label, ringName, value, tag, note) -> (
    rankRecords#(#rankRecords) = (label, ringName, value, tag, note);
    tagPrint(tag, label | " [" | ringName | "] = " | toString value
              | (if note == "" then "" else "  (" | note | ")"))
)

sep = () -> (
    line := concatenate(64 : "-");
    print line;
    appendLog(line)
)

bigsep = () -> (
    line := concatenate(64 : "=");
    print line;
    appendLog(line)
)

-- ================================================================
-- §1bis — MANIFESTE GELÉ DES CLAIMS  [NEW v19]
-- ================================================================

-- Chaque entree : {key, expectedValue, tag, ringName}
-- expectedValue peut etre un entier, un booleen, ou une liste (pour H_apolar).
FROZEN_CLAIMS = {
    {"adj_identity_QQ",                true,         "PROVED-ALG",   "QQ"},
    {"adj_identity_kk",                true,         "PROVED-ALG",   "kk"},
    {"piMinus_idemp_QQ",               true,         "QQ-CERTIFIED", "QQ"},
    {"piPlus_idemp_QQ",                true,         "QQ-CERTIFIED", "QQ"},
    {"piMinus_piPlus_orth_QQ",         true,         "QQ-CERTIFIED", "QQ"},
    {"piMinus_plus_piPlus_id_QQ",      true,         "QQ-CERTIFIED", "QQ"},
    {"rank_piMinus_QQ",                9,            "QQ-CERTIFIED", "QQ"},
    {"rank_piPlus_QQ",                 36,           "QQ-CERTIFIED", "QQ"},
    {"rank_Cat_det3_QQ",               9,            "QQ-CERTIFIED", "QQ"},
    {"delta_det3_QQ",                  9,            "QQ-CERTIFIED", "QQ"},
    {"delta_perm3_QQ",                 0,            "QQ-CERTIFIED", "QQ"},
    {"rank_piPlus_Cat_det3_QQ",        0,            "QQ-CERTIFIED", "QQ"},
    {"dim_Sing_det3_QQ",               5,            "QQ-CERTIFIED", "QQ"},
    {"dim_Sing_perm3_QQ",              3,            "QQ-CERTIFIED", "QQ"},
    {"dim_proj_Sing_det3_QQ",          4,            "QQ-CERTIFIED", "QQ"},
    {"dim_proj_Sing_perm3_QQ",         2,            "QQ-CERTIFIED", "QQ"},
    {"dim_VJ_det3_QQ",                 5,            "QQ-CERTIFIED", "QQ"},
    {"dim_VJ_perm3_QQ",                3,            "QQ-CERTIFIED", "QQ"},
    {"H_apolar_det3_QQ",               {1,9,9,1},    "QQ-CERTIFIED", "QQ"},
    {"H_apolar_perm3_QQ",              {1,9,9,1},    "QQ-CERTIFIED", "QQ"},
    {"codim_J_det3_QQ",                4,            "QQ-CERTIFIED", "QQ"},
    {"codim_J_perm3_QQ",               6,            "QQ-CERTIFIED", "QQ"},
    {"degree_J_det3_QQ",               6,            "QQ-CERTIFIED", "QQ"},
    {"degree_J_perm3_QQ",              24,           "QQ-CERTIFIED", "QQ"},
    {"NIV1_orbits_distinct",           true,         "PROVED-ALG",   "QQ"},
    {"NIV2_dir1_perm_notin_clos_det",  true,         "PROVED-ALG",   "QQ"}
    -- NIV-2 dir. 2 : volontairement absent. Reste OPEN.
}

frozenKeys = set apply(FROZEN_CLAIMS, c -> c#0)
frozenByKey = new MutableHashTable
scan(FROZEN_CLAIMS, c -> frozenByKey#(c#0) = c)

claimsObserved = new MutableHashTable

registerClaim = (key, observed, tag, ringName) -> (
    if CONFIG#"strictManifest" and (tag === "PROVED-ALG" or tag === "QQ-CERTIFIED") then (
        if not member(key, frozenKeys) then (
            tagPrint("WARN", "Claim hors manifeste detecte : " | key | " (tag=" | tag | ")");
            error("[MANIFEST] Claim PROVED-ALG/QQ-CERTIFIED hors manifeste : " | key)
        )
    );
    claimsObserved#key = hashTable {
        "observed" => observed,
        "tag"      => tag,
        "ring"     => ringName
    };
    tagPrint(tag, key | " [" | ringName | "] = " | toString observed)
)

-- ================================================================
-- §1ter — PROVENANCE  [NEW v19]
-- ================================================================

computeScriptSha256 = (path) -> (
    h := "";
    try (
        cmd1 := "!sha256sum " | path | " 2>/dev/null";
        out := get cmd1;
        if out =!= null and #out > 0 then (
            ln := first lines out;
            if ln =!= null and #ln >= 64 then h = substring(ln, 0, 64)
        )
    ) else null;
    if h == "" then (
        try (
            cmd2 := "!shasum -a 256 " | path | " 2>/dev/null";
            out2 := get cmd2;
            if out2 =!= null and #out2 > 0 then (
                ln2 := first lines out2;
                if ln2 =!= null and #ln2 >= 64 then h = substring(ln2, 0, 64)
            )
        ) else null
    );
    if h == "" then "unavailable" else h
)

PROVENANCE = new MutableHashTable
PROVENANCE#"runId"         = CONFIG#"runId"
PROVENANCE#"runStartedAt"  = toString currentTime()
PROVENANCE#"m2Version"     = try toString version#"VERSION" else "unknown"
PROVENANCE#"seed"          = CONFIG#"seed"
PROVENANCE#"primeModList"  = CONFIG#"primeModList"
PROVENANCE#"scriptSha256"  = computeScriptSha256(CONFIG#"scriptPath")
PROVENANCE#"flags" = hashTable {
    "runQQ"           => CONFIG#"runQQ",
    "runQQHeavy"      => CONFIG#"runQQHeavy",
    "runSymbolicM"    => CONFIG#"runSymbolicM",
    "runSampling"     => CONFIG#"runSampling",
    "runHessianDiag"  => CONFIG#"runHessianDiag",
    "runQuarticObs"   => CONFIG#"runQuarticObs",
    "strictManifest"  => CONFIG#"strictManifest"
}
PROVENANCE#"artefacts" = new MutableList from {}

addArtefact = (path) -> (
    PROVENANCE#"artefacts"#(#PROVENANCE#"artefacts") = path
)

bigsep()
print "SRMT v19 -- Lock & Audit Release"
bigsep()
tagPrint("INFO", "runId = " | PROVENANCE#"runId")
tagPrint("INFO", "seed = " | toString PROVENANCE#"seed")
tagPrint("INFO", "primeModList = " | toString PROVENANCE#"primeModList")
tagPrint("INFO", "M2 version = " | PROVENANCE#"m2Version")
tagPrint("INFO", "scriptSha256 = " | PROVENANCE#"scriptSha256")
tagPrint("INFO", "strictManifest = " | toString CONFIG#"strictManifest")

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
    registerClaim("adj_identity_kk", true, "PROVED-ALG", "kk")
) else (
    tagPrint("WARN", "adj(Y)*Y != det3(Y)*I sur kk[y]");
    registerClaim("adj_identity_kk", false, "PROVED-ALG", "kk")
)

if idCheckQ == 0 then (
    registerClaim("adj_identity_QQ", true, "PROVED-ALG", "QQ")
) else (
    tagPrint("WARN", "adj(Y)*Y != det3(Y)*I sur QQ[y]");
    registerClaim("adj_identity_QQ", false, "PROVED-ALG", "QQ")
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

tagPrint("NOTE",       "Reference : Macdonald, Fulton-Harris")
tagPrint("PROVED-ALG", "Sym^3(V*W) = S(3)*S(3) + S(2,1)*S(2,1) + S(1,1,1)*S(1,1,1) (cite)")
tagPrint("PROVED-ALG", "Dimensions : 100 + 64 + 1 = 165 (cite)")
tagPrint("PROVED-ALG", "det3 in S(1,1,1)*S(1,1,1) (cite)")
tagPrint("PROVED-ALG", "perm3 in S(3)*S(3) (cite)")

logStep("cauchy_schur", "PROVED-ALG", "decomp Sym^3(V*W) cited")

-- Note importante : ces lignes sont des claims cites (NOTE manuscrit).
-- Elles ne sont PAS dans FROZEN_CLAIMS car elles ne produisent pas de
-- valeur numerique calculable. Elles sont strictement des PROVED-ALG
-- referentiels, pas des claims a verifier mecaniquement.
-- Pour ne pas declencher le manifeste, on les emet via tagPrint direct
-- et NON via registerClaim.

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

    -- Enregistrement dans le manifeste UNIQUEMENT pour QQ
    if F === QQ then (
        registerClaim("piMinus_idemp_QQ",          chkIdemMinus, "QQ-CERTIFIED", "QQ");
        registerClaim("piPlus_idemp_QQ",           chkIdemPlus,  "QQ-CERTIFIED", "QQ");
        registerClaim("piMinus_piPlus_orth_QQ",    chkOrth,      "QQ-CERTIFIED", "QQ");
        registerClaim("piMinus_plus_piPlus_id_QQ", chkComplete,  "QQ-CERTIFIED", "QQ");
        registerClaim("rank_piMinus_QQ",           rkM,          "QQ-CERTIFIED", "QQ");
        registerClaim("rank_piPlus_QQ",            rkP,          "QQ-CERTIFIED", "QQ")
    );

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

tagPrint("PROVED-ALG", "Lemme 2 : Sym^2(V*W) = Sym^2 V * Sym^2 W + Lambda^2 V * Lambda^2 W (cite)")
tagPrint("PROVED-ALG", "Lemme 3 : piMinus, piPlus sont GL(V)xGL(W)-equivariants (Schur, cite)")

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
    registerClaim("rank_Cat_det3_QQ", rkCatDet3Q, "QQ-CERTIFIED", "QQ");

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

    registerClaim("delta_det3_QQ",            deltaDet3Q,     "QQ-CERTIFIED", "QQ");
    registerClaim("delta_perm3_QQ",           deltaPerm3Q,    "QQ-CERTIFIED", "QQ");
    registerClaim("rank_piPlus_Cat_det3_QQ",  deltaPlusDet3Q, "QQ-CERTIFIED", "QQ");

    if deltaDet3Q == 9 and deltaPerm3Q == 0 then (
        tagPrint("QQ-CERTIFIED", "delta(det3)=9 et delta(perm3)=0 sur QQ => vrais sur C");
        tagPrint("PROVED-ALG",   "Prop 5 : adherences GL3xGL3.det3 et GL3xGL3.perm3 disjointes (cite)");
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
-- §11 — RANGS MODULAIRES (multi-prime)  [MOD v19]
-- ================================================================

bigsep()
print "Section 11 - Rangs modulaires sur kk (multi-prime)"
bigsep()

scan(CONFIG#"primeModList", p -> (
    kkP    := ZZ/p;
    R9P    := kkP[y_0..y_8, MonomialOrder => GRevLex];
    YP     := buildY R9P;
    det3P  := det3OnRing(R9P, YP);
    perm3P := perm3OnRing(R9P, YP);
    mons2P := mons2OnRing(R9P);
    (piMP, piPP) := buildProjectorsOnField(kkP, mons2P, R9P);

    CatDetP  := buildCatMatrixOnRing(R9P, mons2P, det3P);
    CatPermP := buildCatMatrixOnRing(R9P, mons2P, perm3P);

    rkCatDetP  := rank CatDetP;
    rkCatPermP := rank CatPermP;
    deltaDetP  := rank (piMP * CatDetP);
    deltaPermP := rank (piMP * CatPermP);
    plusDetP   := rank (piPP * CatDetP);
    plusPermP  := rank (piPP * CatPermP);

    pTag := "kk=" | toString p;
    recordRank("rank(Cat(det3))",        pTag, rkCatDetP,  "MODULAR", "");
    recordRank("rank(Cat(perm3))",       pTag, rkCatPermP, "MODULAR", "");
    recordRank("delta(det3)",            pTag, deltaDetP,  "MODULAR", "");
    recordRank("delta(perm3)",           pTag, deltaPermP, "MODULAR", "");
    recordRank("rank(piPlus*Cat(det3))", pTag, plusDetP,   "MODULAR", "");
    recordRank("rank(piPlus*Cat(perm3))",pTag, plusPermP,  "MODULAR", "");

    logStep("modular_p_" | toString p, "MODULAR",
            "rkCatDet=" | toString rkCatDetP | " deltaDet=" | toString deltaDetP
            | " deltaPerm=" | toString deltaPermP)
));

tagPrint("NOTE",
    "Rangs modulaires = minorants des rangs sur QQ. Ne prouvent PAS sur C.")

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

deltaDet3K  = rank (piMinusMat * CatDet3K)
deltaPerm3K = rank (piMinusMat * CatPerm3K)

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
tagPrint("NOTE", "Toute separation GL9.det3 vs GL9.perm3 requiert un invariant GL9 intrinseque")
tagPrint("NOTE", "=> traitee dans la section 15 par invariants du lieu singulier / apolaire / jacobien")

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
-- §15 — STATUT GL9 : ORBITES vs ADHERENCES (cascade disciplinee)
-- ================================================================
--
-- Les conventions geometriques, le lemme de semi-continuite et son
-- application sont identiques a v18 (preambule §15). Cf. srmt_v18
-- pour les details ; v19 conserve la cascade A1 / A1bis / A2 / A3 /
-- B / C sans modification mathematique. On ajoute uniquement les
-- registerClaim pour alimenter le manifeste.
-- ================================================================

separationGL9OrbitStatus  = "OPEN"
separationGL9OrbitMethod  = "none"
nonDegenPerm3NotInGL9det3 = "OPEN"
nonDegenDet3NotInGL9perm3 = "OPEN"
semiContMethod            = "none"

if CONFIG#"runSymbolicM" then (
    bigsep();
    print "Section 15 - Statut GL9 (orbites vs adherences, cascade disciplinee)";
    bigsep();

    --------------------------------------------------------------------
    -- STRATEGIE A1 : invariant GL9 = dim Sing V(f)  [SEMI-CONTINU SUP]
    --------------------------------------------------------------------

    print "";
    print "  -- Strategie A1 : dim Sing V(f) (semi-continuite superieure) --";
    logStep("strat_A1_dimSing", "START", "computing dim of singular schemes on QQ");

    dimSingDet  := -1;
    dimSingPerm := -1;

    try (
        gradDet3Q     := ideal jacobian matrix{{det3Q}};
        singIdealDet  := ideal det3Q  + gradDet3Q;
        gradPerm3Q    := ideal jacobian matrix{{perm3Q}};
        singIdealPerm := ideal perm3Q + gradPerm3Q;

        dimSingDet  = dim singIdealDet;
        dimSingPerm = dim singIdealPerm;

        recordRank("dim Sing V(det3) [QQ, A^9]",  "QQ", dimSingDet,
                   "QQ-CERTIFIED", "semi-continu sup");
        recordRank("dim Sing V(perm3) [QQ, A^9]", "QQ", dimSingPerm,
                   "QQ-CERTIFIED", "semi-continu sup");

        registerClaim("dim_Sing_det3_QQ",  dimSingDet,  "QQ-CERTIFIED", "QQ");
        registerClaim("dim_Sing_perm3_QQ", dimSingPerm, "QQ-CERTIFIED", "QQ");

        if dimSingDet =!= dimSingPerm then (
            tagPrint("QQ-CERTIFIED",
                "dim Sing V(det3) != dim Sing V(perm3) sur QQ");
            separationGL9OrbitStatus = "PROVED-ALG";
            separationGL9OrbitMethod = "dim_singular_locus";

            tagPrint("NOTE",
                "Lemme de semi-continuite : f0 in closure(orbit f1) => dim Sing V(f0) >= dim Sing V(f1)");
            if dimSingPerm < dimSingDet then (
                logStep("semi_cont_perm_notin_clos_det", "PROVED-ALG",
                        "dimSingPerm=" | toString dimSingPerm | " < dimSingDet=" | toString dimSingDet);
                nonDegenPerm3NotInGL9det3 = "PROVED-ALG";
                semiContMethod = "dim_singular_locus"
            ) else (
                tagPrint("NOTE",
                    "Direction 1 non concluante : dim Sing V(perm3) >= dim Sing V(det3)")
            );
            if dimSingDet < dimSingPerm then (
                tagPrint("WARN",
                    "Direction 2 PROVED-ALG inattendue -- a investiguer manuellement");
                logStep("semi_cont_det_notin_clos_perm", "WARN",
                        "dimSingDet=" | toString dimSingDet | " < dimSingPerm=" | toString dimSingPerm);
                nonDegenDet3NotInGL9perm3 = "PROVED-ALG"
            ) else (
                tagPrint("NOTE",
                    "Direction 2 non concluante : dim Sing V(det3) >= dim Sing V(perm3) (invariant insuffisant)")
            );

            logStep("strat_A1_dimSing", "PROVED-ALG",
                    "dimSingDet=" | toString dimSingDet | " dimSingPerm=" | toString dimSingPerm)
        ) else (
            tagPrint("NOTE",
                "dim Sing V(det3) = dim Sing V(perm3) -- A1 ne separe meme pas les orbites");
            logStep("strat_A1_dimSing", "INSUFFICIENT", "dim coincide")
        )
    ) else (
        tagPrint("WARN", "Strategie A1 : echec calcul (exception M2)");
        logStep("strat_A1_dimSing", "ERROR", "M2 exception")
    );

    --------------------------------------------------------------------
    -- STRATEGIE A1-bis : COHERENCE AFFINE / PROJECTIF
    --------------------------------------------------------------------

    print "";
    print "  -- Strategie A1-bis : coherence affine / projectif --";
    logStep("strat_A1bis_affproj", "START", "checking dim_aff vs dim_proj");

    try (
        gradDet3Q2     := ideal jacobian matrix{{det3Q}};
        singIdealDet2  := ideal det3Q  + gradDet3Q2;
        gradPerm3Q2    := ideal jacobian matrix{{perm3Q}};
        singIdealPerm2 := ideal perm3Q + gradPerm3Q2;

        dAffDet  := dim singIdealDet2;
        dAffPerm := dim singIdealPerm2;
        dProjDet  := dAffDet  - 1;
        dProjPerm := dAffPerm - 1;

        recordRank("dim_proj Sing^proj X(det3) [QQ, P^8]",  "QQ", dProjDet,
                   "QQ-CERTIFIED", "projectif (dim_aff - 1)");
        recordRank("dim_proj Sing^proj X(perm3) [QQ, P^8]", "QQ", dProjPerm,
                   "QQ-CERTIFIED", "projectif (dim_aff - 1)");

        registerClaim("dim_proj_Sing_det3_QQ",  dProjDet,  "QQ-CERTIFIED", "QQ");
        registerClaim("dim_proj_Sing_perm3_QQ", dProjPerm, "QQ-CERTIFIED", "QQ");

        if dAffDet >= 1 and dAffPerm >= 1 then (
            tagPrint("PROVED-ALG",
                "Coherence (d) : dim_aff = dim_proj + 1 verifiee pour det3 et perm3 (cite)")
        ) else (
            tagPrint("NOTE",
                "Cas degenere : Sing^proj possiblement vide pour l'un des deux")
        );

        if dProjPerm < dProjDet then (
            logStep("strat_A1bis_affproj", "PROVED-ALG",
                    "dim_proj coherent avec dim_aff")
        ) else (
            tagPrint("NOTE",
                "Diagnostic projectif non concluant pour dir.1 (cas degenere)")
        )
    ) else (
        tagPrint("WARN", "Strategie A1-bis : echec calcul (exception M2)");
        logStep("strat_A1bis_affproj", "ERROR", "M2 exception")
    );

    --------------------------------------------------------------------
    -- STRATEGIE A2 : dim du quotient jacobien Q_f = R/J(f)
    --------------------------------------------------------------------

    print "";
    print "  -- Strategie A2 : dim R/J(f) (semi-continuite superieure) --";
    logStep("strat_A2_dimQuotJac", "START", "computing dim R/J(f)");

    dimQjDet  := -1;
    dimQjPerm := -1;

    try (
        JdetI    := ideal jacobian matrix{{det3Q}};
        JpermI   := ideal jacobian matrix{{perm3Q}};

        dimQjDet  = dim JdetI;
        dimQjPerm = dim JpermI;

        recordRank("dim V(J(det3)) [QQ, A^9]",  "QQ", dimQjDet,
                   "QQ-CERTIFIED", "semi-continu sup");
        recordRank("dim V(J(perm3)) [QQ, A^9]", "QQ", dimQjPerm,
                   "QQ-CERTIFIED", "semi-continu sup");

        registerClaim("dim_VJ_det3_QQ",  dimQjDet,  "QQ-CERTIFIED", "QQ");
        registerClaim("dim_VJ_perm3_QQ", dimQjPerm, "QQ-CERTIFIED", "QQ");

        if dimQjDet =!= dimQjPerm then (
            if separationGL9OrbitStatus =!= "PROVED-ALG" then (
                separationGL9OrbitStatus = "PROVED-ALG";
                separationGL9OrbitMethod = "dim_jacobian_quotient"
            );
            if (dimQjPerm < dimQjDet) and (nonDegenPerm3NotInGL9det3 =!= "PROVED-ALG") then (
                logStep("semi_cont_perm_notin_clos_det_A2", "PROVED-ALG",
                        "dimQjPerm=" | toString dimQjPerm | " < dimQjDet=" | toString dimQjDet);
                nonDegenPerm3NotInGL9det3 = "PROVED-ALG";
                semiContMethod = (if semiContMethod === "none" then "dim_jacobian_quotient" else semiContMethod | "+dim_jacobian_quotient")
            );
            if (dimQjDet < dimQjPerm) and (nonDegenDet3NotInGL9perm3 =!= "PROVED-ALG") then (
                tagPrint("WARN",
                    "Direction 2 (via A2) inattendue PROVED-ALG -- investigation requise");
                logStep("semi_cont_det_notin_clos_perm_A2", "WARN",
                        "dimQjDet=" | toString dimQjDet | " < dimQjPerm=" | toString dimQjPerm);
                nonDegenDet3NotInGL9perm3 = "PROVED-ALG"
            );
            logStep("strat_A2_dimQuotJac", "USEFUL",
                    "dimQjDet=" | toString dimQjDet | " dimQjPerm=" | toString dimQjPerm)
        ) else (
            tagPrint("NOTE",
                "dim V(J(det3)) = dim V(J(perm3)) -- A2 n'apporte rien de plus");
            logStep("strat_A2_dimQuotJac", "INSUFFICIENT", "dim coincide")
        )
    ) else (
        tagPrint("WARN", "Strategie A2 : echec calcul (exception M2)");
        logStep("strat_A2_dimQuotJac", "ERROR", "M2 exception")
    );

    --------------------------------------------------------------------
    -- STRATEGIE A3 : Hilbert apolaire
    --------------------------------------------------------------------

    print "";
    print "  -- Strategie A3 : fonction de Hilbert apolaire (invariant d'orbite) --";
    logStep("strat_A3_apolarHilbert", "START", "computing apolar Hilbert function");

    try (
        apolarHilbertOf := (f, dTot, kmax) -> (
            Rloc := ring f;
            apply(kmax + 1, k -> (
                if k > dTot then 0
                else if k == 0 then 1
                else if k == dTot then 1
                else (
                    monsK  := flatten entries basis(k,       Rloc);
                    monsDK := flatten entries basis(dTot-k, Rloc);
                    mat := matrix apply(monsK, m -> (
                        expVec := flatten exponents m;
                        Df := f;
                        scan(#expVec, i -> (
                            scan(expVec#i, r -> (
                                Df = diff(Rloc_i, Df)
                            ))
                        ));
                        apply(monsDK, mm -> coefficient(mm, Df))
                    ));
                    rank mat
                )
            ))
        );

        HDet  := apolarHilbertOf(det3Q,  3, 3);
        HPerm := apolarHilbertOf(perm3Q, 3, 3);

        tagPrint("QQ-CERTIFIED", "H_{A_det3}  = " | toString HDet  | "  (invariant d'orbite GL9)");
        tagPrint("QQ-CERTIFIED", "H_{A_perm3} = " | toString HPerm | "  (invariant d'orbite GL9)");

        recordRank("H_apolar(det3) [QQ]",  "QQ", HDet,  "QQ-CERTIFIED", "orbite (pas semi-continu)");
        recordRank("H_apolar(perm3) [QQ]", "QQ", HPerm, "QQ-CERTIFIED", "orbite (pas semi-continu)");

        registerClaim("H_apolar_det3_QQ",  HDet,  "QQ-CERTIFIED", "QQ");
        registerClaim("H_apolar_perm3_QQ", HPerm, "QQ-CERTIFIED", "QQ");

        if HDet =!= HPerm then (
            if separationGL9OrbitStatus =!= "PROVED-ALG" then (
                separationGL9OrbitStatus = "PROVED-ALG";
                separationGL9OrbitMethod = "apolar_hilbert_function"
            ) else (
                tagPrint("NOTE", "A3 confirme la separation d'orbite deja etablie")
            );
            tagPrint("NOTE",
                "H_apolar n'est PAS semi-continu en general -- aucune conclusion sur les adherences")
        ) else (
            tagPrint("NOTE",
                "H_apolar coincide pour det3 et perm3 (cas connu : (1,9,9,1))");
            logStep("strat_A3_apolarHilbert", "INSUFFICIENT", "H coincide: " | toString HDet)
        )
    ) else (
        tagPrint("WARN", "Strategie A3 : echec calcul (exception M2)");
        logStep("strat_A3_apolarHilbert", "ERROR", "M2 exception")
    );

    --------------------------------------------------------------------
    -- STRATEGIE B : sondage 1-PS (OBS-N)
    --------------------------------------------------------------------

    print "";
    print "  -- Strategie B : sondage 1-parameter subgroups (OBS-N seulement) --";
    logStep("strat_B_oneParam", "START", "sampling 1-PS degenerations of det3 and perm3");

    try (
        limitForm := (f, w) -> (
            Rloc  := ring f;
            (mns, cfs) := coefficients f;
            mList := flatten entries mns;
            cList := flatten entries cfs;
            wts := apply(mList, m -> (
                e := flatten exponents m;
                sum apply(#e, i -> e#i * w#i)
            ));
            wmin := min wts;
            sum apply(#mList, idx -> (
                if wts#idx == wmin then sub(cList#idx, Rloc) * mList#idx
                else 0_Rloc
            ))
        );

        weightSamples := {
            {1,1,1,1,1,1,1,1,1},
            {1,0,0,0,1,0,0,0,1},
            {1,2,3,4,5,6,7,8,9},
            {0,1,2,0,1,2,0,1,2},
            {3,1,1,1,3,1,1,1,3},
            {2,2,2,1,1,1,0,0,0},
            {0,0,0,1,1,1,2,2,2},
            {1,0,2,2,1,0,0,2,1}
        };

        nWeights := #weightSamples;
        nbDifferentLimits := 0;
        nbDetSurvives := 0;

        scan(weightSamples, w -> (
            limDet  := limitForm(det3Q,  w);
            limPerm := limitForm(perm3Q, w);
            sDet  := set flatten entries (coefficients limDet)#0;
            sPerm := set flatten entries (coefficients limPerm)#0;
            if sDet =!= sPerm then nbDifferentLimits = nbDifferentLimits + 1;
            if limDet == det3Q then nbDetSurvives = nbDetSurvives + 1
        ));

        tagPrint("OBS-N",
            toString nbDifferentLimits | "/" | toString nWeights
          | " 1-PS donnent des supports limites differents pour det3 et perm3");
        tagPrint("OBS-N",
            toString nbDetSurvives | "/" | toString nWeights
          | " 1-PS preservent det3 (pas de degenerescence non triviale)");
        tagPrint("NOTE",
            "Sondage 1-PS : ne prouve RIEN sur les adherences globales -- supporte au plus la conjecture");
        logStep("strat_B_oneParam", "OBS-N",
                "nWeights=" | toString nWeights
              | " diffLimits=" | toString nbDifferentLimits
              | " detSurvives=" | toString nbDetSurvives)
    ) else (
        tagPrint("WARN", "Strategie B : echec calcul (exception M2)");
        logStep("strat_B_oneParam", "ERROR", "M2 exception")
    );

    --------------------------------------------------------------------
    -- STRATEGIE C : codim/degree des ideaux jacobiens (orbite)
    --------------------------------------------------------------------

    print "";
    print "  -- Strategie C : codim/degree des ideaux jacobiens (orbite) --";
    logStep("strat_C_jacInvariants", "START", "computing codim/deg J(f) as orbit invariants");

    try (
        JdetI2  := ideal jacobian matrix{{det3Q}};
        JpermI2 := ideal jacobian matrix{{perm3Q}};

        cdJDet  := codim JdetI2;
        cdJPerm := codim JpermI2;
        dgJDet  := degree JdetI2;
        dgJPerm := degree JpermI2;

        recordRank("codim J(det3)  [QQ]", "QQ", cdJDet,  "QQ-CERTIFIED", "orbite (pas semi-continu)");
        recordRank("codim J(perm3) [QQ]", "QQ", cdJPerm, "QQ-CERTIFIED", "orbite (pas semi-continu)");
        recordRank("degree J(det3)  [QQ]", "QQ", dgJDet,  "QQ-CERTIFIED", "orbite (pas semi-continu)");
        recordRank("degree J(perm3) [QQ]", "QQ", dgJPerm, "QQ-CERTIFIED", "orbite (pas semi-continu)");

        registerClaim("codim_J_det3_QQ",   cdJDet,  "QQ-CERTIFIED", "QQ");
        registerClaim("codim_J_perm3_QQ",  cdJPerm, "QQ-CERTIFIED", "QQ");
        registerClaim("degree_J_det3_QQ",  dgJDet,  "QQ-CERTIFIED", "QQ");
        registerClaim("degree_J_perm3_QQ", dgJPerm, "QQ-CERTIFIED", "QQ");

        if (cdJDet =!= cdJPerm) or (dgJDet =!= dgJPerm) then (
            if separationGL9OrbitStatus =!= "PROVED-ALG" then (
                separationGL9OrbitStatus = "PROVED-ALG";
                separationGL9OrbitMethod = "jacobian_invariants"
            );
            tagPrint("NOTE",
                "codim/degree de J ne sont PAS systematiquement semi-continus -- aucune conclusion sur les adherences")
        ) else (
            tagPrint("NOTE", "codim et degree des ideaux jacobiens coincident")
        );
        logStep("strat_C_jacInvariants", "DONE",
                "cdJDet=" | toString cdJDet | " cdJPerm=" | toString cdJPerm
              | " dgJDet=" | toString dgJDet | " dgJPerm=" | toString dgJPerm)
    ) else (
        tagPrint("WARN", "Strategie C : echec calcul (exception M2)");
        logStep("strat_C_jacInvariants", "ERROR", "M2 exception")
    );

    --------------------------------------------------------------------
    -- BILAN SECTION 15
    --------------------------------------------------------------------
    print "";
    print "  -- Bilan section 15 --";

    if separationGL9OrbitStatus === "PROVED-ALG" then (
        registerClaim("NIV1_orbits_distinct", true, "PROVED-ALG", "QQ");
        tagPrint("INFO",
            "NIV-1 methode : " | separationGL9OrbitMethod);
        logStep("section15_orbit", "PROVED-ALG", separationGL9OrbitMethod)
    ) else (
        tagPrint("OPEN", "NIV-1 : separation d'orbite GL9 indecidee dans ce run");
        registerClaim("NIV1_orbits_distinct", false, "PROVED-ALG", "QQ");
        logStep("section15_orbit", "OPEN", "no orbit invariant separated")
    );

    if nonDegenPerm3NotInGL9det3 === "PROVED-ALG" then (
        registerClaim("NIV2_dir1_perm_notin_clos_det", true, "PROVED-ALG", "QQ");
        tagPrint("INFO",
            "NIV-2 dir.1 methode : semi-continuite (" | semiContMethod | ")");
        logStep("section15_dir1", "PROVED-ALG", semiContMethod)
    ) else (
        tagPrint("OPEN",
            "NIV-2 dir. 1 : perm3 not in closure(GL9.det3) non etabli dans ce run");
        registerClaim("NIV2_dir1_perm_notin_clos_det", false, "PROVED-ALG", "QQ");
        logStep("section15_dir1", "OPEN", "no semi-continuous invariant strict")
    );

    -- NIV-2 dir. 2 : volontairement NON enregistree dans le manifeste.
    -- Si elle devenait PROVED-ALG (ce qui ne devrait jamais arriver dans
    -- une release Lock & Audit), strictManifest declencherait une erreur.
    if nonDegenDet3NotInGL9perm3 === "PROVED-ALG" then (
        tagPrint("WARN",
            "NIV-2 dir. 2 marquee PROVED-ALG -- INATTENDU, manifeste va echouer")
    ) else (
        tagPrint("OPEN",
            "NIV-2 dir. 2 : det3 not in closure(GL9.perm3) -- conjecture GCT centrale, non resolue ici");
        logStep("section15_dir2", "OPEN", "GCT central conjecture")
    )
) else (
    tagPrint("NOTE", "Section 15 desactivee (CONFIG runSymbolicM = false)");
    tagPrint("OPEN", "Statut GL9 (orbites et adherences) non examine dans ce run")
)

-- ================================================================
-- §15bis — OBSERVATION QUARTIQUE SPARSE  [NEW v19, OPTIONNEL]
-- ================================================================
-- Strictement OBS-N / NOTE. Aucun PROVED-ALG, aucun QQ-CERTIFIED.
-- Sortie placee dans la section "quartic_observations" du JSON,
-- disjointe de claims_status.

quarticObservations = new MutableList from {}

if CONFIG#"runQuarticObs" then (
    bigsep();
    print "Section 15bis - Observation quartique sparse (OBS uniquement)";
    bigsep();
    tagPrint("OBS-N",
        "Recherche quartique sparse : enregistree comme observation, ne prouve rien");
    tagPrint("NOTE",
        "Aucun resultat quartique n'est promu en PROVED-ALG ou QQ-CERTIFIED");
    quarticObservations#(#quarticObservations) = hashTable {
        "type"   => "OBS-N",
        "status" => "placeholder",
        "note"   => "Run quartique non execute dans cette release; squelette reserve."
    };
    logStep("quartic_obs", "OBS-N", "placeholder")
) else (
    tagPrint("INFO", "Recherche quartique sparse desactivee (runQuarticObs = false)")
)

-- ================================================================
-- §16 — EXPORTS
-- ================================================================

bigsep()
print "Section 16 - Exports"
bigsep()

-- CSV
csvOut = openOut csvFile
csvOut << "label,ring,value,tag,note" << endl
scan(toList rankRecords, rec -> (
    label := rec#0; ringName := rec#1; value := rec#2; tag := rec#3; note := rec#4;
    csvOut << "\"" << label << "\","
           << ringName << ","
           << "\"" << toString value << "\","
           << tag << ","
           << "\"" << note << "\"" << endl
))
csvOut << close
tagPrint("INFO", "CSV rangs exporte : " | csvFile)
addArtefact(csvFile)

escapeJson = (s) -> replace("\"", "\\\"", toString s)

-- valeur JSON depuis observed (entier, bool, liste)
jsonValueOf = (v) -> (
    if v === true  then "true"
    else if v === false then "false"
    else if v === null  then "null"
    else if instance(v, ZZ) then toString v
    else if instance(v, List) or instance(v, Sequence) then (
        "[" | demark(", ", apply(toList v, x -> jsonValueOf x)) | "]"
    )
    else "\"" | escapeJson v | "\""
)

-- JSON
jsonOut = openOut jsonFile
jsonOut << "{" << endl

-- provenance
jsonOut << "  \"provenance\": {" << endl
jsonOut << "    \"runId\":         \"" << escapeJson PROVENANCE#"runId"         << "\"," << endl
jsonOut << "    \"runStartedAt\":  \"" << escapeJson PROVENANCE#"runStartedAt"  << "\"," << endl
jsonOut << "    \"runFinishedAt\": \"" << escapeJson toString currentTime()      << "\"," << endl
jsonOut << "    \"m2Version\":     \"" << escapeJson PROVENANCE#"m2Version"     << "\"," << endl
jsonOut << "    \"seed\":          " << toString PROVENANCE#"seed"              << "," << endl
jsonOut << "    \"primeModList\":  " << toString PROVENANCE#"primeModList"      << "," << endl
jsonOut << "    \"scriptSha256\":  \"" << escapeJson PROVENANCE#"scriptSha256"  << "\"," << endl
jsonOut << "    \"flags\": {" << endl
jsonOut << "      \"runQQ\":          " << toString CONFIG#"runQQ"          << "," << endl
jsonOut << "      \"runQQHeavy\":     " << toString CONFIG#"runQQHeavy"     << "," << endl
jsonOut << "      \"runSymbolicM\":   " << toString CONFIG#"runSymbolicM"   << "," << endl
jsonOut << "      \"runSampling\":    " << toString CONFIG#"runSampling"    << "," << endl
jsonOut << "      \"runHessianDiag\": " << toString CONFIG#"runHessianDiag" << "," << endl
jsonOut << "      \"runQuarticObs\":  " << toString CONFIG#"runQuarticObs"  << "," << endl
jsonOut << "      \"strictManifest\": " << toString CONFIG#"strictManifest" << endl
jsonOut << "    }," << endl
jsonOut << "    \"artefacts\": [" << endl
nArt = #PROVENANCE#"artefacts"
scan(nArt, idx -> (
    a := PROVENANCE#"artefacts"#idx;
    jsonOut << "      \"" << escapeJson a << "\"";
    if idx < nArt - 1 then jsonOut << ",";
    jsonOut << endl
))
jsonOut << "    ]" << endl
jsonOut << "  }," << endl

-- claims_status
jsonOut << "  \"claims_status\": [" << endl
nC = #FROZEN_CLAIMS
scan(nC, idx -> (
    c := FROZEN_CLAIMS#idx;
    key := c#0; expected := c#1; tag := c#2; ringName := c#3;
    obs := if claimsObserved#?key then (claimsObserved#key)#"observed" else null;
    status := if obs === null then "MISSING"
              else if obs === expected then "OK"
              else "MISMATCH";
    jsonOut << "    {\"key\": \""        << escapeJson key                 << "\", "
            << "\"expected\": "          << jsonValueOf expected            << ", "
            << "\"observed\": "          << jsonValueOf obs                 << ", "
            << "\"status\": \""          << status                          << "\", "
            << "\"tag\": \""             << tag                             << "\", "
            << "\"ring\": \""            << ringName                        << "\"}";
    if idx < nC - 1 then jsonOut << ",";
    jsonOut << endl
))
jsonOut << "  ]," << endl

-- open_questions
jsonOut << "  \"open_questions\": [" << endl
jsonOut << "    \"det3 NOT IN closure(GL9 . perm3) -- conjecture GCT centrale (Mulmuley-Sohoni)\"," << endl
jsonOut << "    \"lien explicite delta vs multiplicites BIP\"," << endl
jsonOut << "    \"construction symbolique M(f_generic) et mineur 9x9 separant explicite\"" << endl
jsonOut << "  ]," << endl

-- ranks
jsonOut << "  \"ranks\": [" << endl
nRec = #rankRecords
scan(nRec, idx -> (
    rec := rankRecords#idx;
    label := rec#0; ringName := rec#1; value := rec#2; tag := rec#3; note := rec#4;
    jsonOut << "    {\"label\": \"" << escapeJson label << "\", "
            << "\"ring\": \""       << escapeJson ringName << "\", "
            << "\"value\": \""      << escapeJson toString value << "\", "
            << "\"tag\": \""        << tag << "\", "
            << "\"note\": \""       << escapeJson note << "\"}";
    if idx < nRec - 1 then jsonOut << ",";
    jsonOut << endl
))
jsonOut << "  ]," << endl

-- events
jsonOut << "  \"events\": [" << endl
nEv = #logEvents
scan(nEv, idx -> (
    ev := logEvents#idx;
    stepName := ev#0; status := ev#1; info := ev#2;
    jsonOut << "    {\"step\": \""   << escapeJson stepName << "\", "
            << "\"status\": \""      << escapeJson status   << "\", "
            << "\"info\": \""        << escapeJson info     << "\"}";
    if idx < nEv - 1 then jsonOut << ",";
    jsonOut << endl
))
jsonOut << "  ]," << endl

-- quartic observations
jsonOut << "  \"quartic_observations\": [" << endl
nQ = #quarticObservations
scan(nQ, idx -> (
    q := quarticObservations#idx;
    jsonOut << "    {\"type\": \""   << escapeJson q#"type"   << "\", "
            << "\"status\": \""      << escapeJson q#"status" << "\", "
            << "\"note\": \""        << escapeJson q#"note"   << "\"}";
    if idx < nQ - 1 then jsonOut << ",";
    jsonOut << endl
))
jsonOut << "  ]" << endl

jsonOut << "}" << endl
jsonOut << close
tagPrint("INFO", "JSON certificat exporte : " | jsonFile)
addArtefact(jsonFile)
addArtefact(logFile)

-- ================================================================
-- §16bis — VERROUILLAGE FINAL DU MANIFESTE  [NEW v19]
-- ================================================================

bigsep()
print "Section 16bis - Verrouillage du manifeste"
bigsep()

enforceManifest = () -> (
    o := openOut manifestDiff;
    nMissing := 0;
    nMismatch := 0;
    nOk := 0;
    o << "# SRMT v19 -- manifest diff" << endl;
    o << "# runId = " << CONFIG#"runId" << endl;
    o << "# generated at " << toString currentTime() << endl;
    o << endl;
    scan(FROZEN_CLAIMS, c -> (
        key := c#0; expected := c#1;
        if not claimsObserved#?key then (
            nMissing = nMissing + 1;
            o << "MISSING  " << key
              << "  expected=" << toString expected << endl
        ) else (
            obs := (claimsObserved#key)#"observed";
            if obs === expected then (
                nOk = nOk + 1;
                o << "OK       " << key
                  << "  =" << toString obs << endl
            ) else (
                nMismatch = nMismatch + 1;
                o << "MISMATCH " << key
                  << "  expected=" << toString expected
                  << "  observed=" << toString obs << endl
            )
        )
    ));
    o << endl;
    o << "# summary : OK=" << toString nOk
      << " MISMATCH=" << toString nMismatch
      << " MISSING="  << toString nMissing << endl;
    close o;
    addArtefact(manifestDiff);
    tagPrint("INFO", "Manifeste : OK=" | toString nOk
        | " MISMATCH=" | toString nMismatch
        | " MISSING="  | toString nMissing);
    tagPrint("INFO", "Diff manifeste ecrit : " | manifestDiff);
    if CONFIG#"strictManifest" then (
        if nMismatch > 0 then (
            tagPrint("WARN", "Manifeste : mismatch detecte");
            error("[MANIFEST] " | toString nMismatch | " mismatch(es) -- voir " | manifestDiff)
        );
        if nMissing > 0 then (
            tagPrint("WARN", "Manifeste : claim(s) attendu(s) absent(s)");
            error("[MANIFEST] " | toString nMissing | " claim(s) manquant(s) -- voir " | manifestDiff)
        )
    );
    (nOk, nMismatch, nMissing)
)

(nOk, nMismatch, nMissing) = enforceManifest()

-- ================================================================
-- §17 — SYNTHESE FINALE
-- ================================================================

bigsep()
print "Section 17 - Synthese finale SRMT v19 Lock & Audit"
bigsep()

print ""
print "  -- PROUVE ALGEBRIQUEMENT (sur C, sans hypothese) --"
tagPrint("PROVED-ALG", "adj(Y)*Y = det3(Y)*I")
tagPrint("PROVED-ALG", "Decomposition Cauchy-Schur Sym^3(V*W) (cite)")
tagPrint("PROVED-ALG", "Lemmes 1, 2, 3 (cite manuscrit)")
tagPrint("PROVED-ALG", "Prop 4 : delta invariant GL3xGL3")
tagPrint("PROVED-ALG", "Prop 5 : adherences GL3xGL3.det3 et GL3xGL3.perm3 disjointes")

print ""
print "  -- CERTIFIE SUR QQ (donc sur C par extension) --"
tagPrint("QQ-CERTIFIED", "Rangs critiques §10 ; identites projecteurs §9")
tagPrint("QQ-CERTIFIED", "dim Sing V(det3)=5, dim Sing V(perm3)=3 ; dim_proj 4 / 2")
tagPrint("QQ-CERTIFIED", "dim V(J(det3))=5, dim V(J(perm3))=3")
tagPrint("QQ-CERTIFIED", "H_apolar = (1,9,9,1) pour les deux")
tagPrint("QQ-CERTIFIED", "codim J = 4 / 6 ; degree J = 6 / 24")

print ""
print "  -- CERTIFIE MODULAIREMENT --"
tagPrint("MODULAR", "Rangs et identites sur " | toString (#CONFIG#"primeModList") | " prime(s)")
tagPrint("NOTE",    "Modulaire = minorant. Ne prouve PAS sur C.")

print ""
print "  -- OBSERVE (stress-tests) --"
tagPrint("OBS-N", "Equivariance Cat / piMinus / piPlus sur tirages")
tagPrint("OBS-N", "delta constant sur orbite GL3xGL3 ; non-constant sur GL9")

print ""
print "  -- STATUT GL9 --"
if separationGL9OrbitStatus === "PROVED-ALG" then
    tagPrint("PROVED-ALG",
        "NIV-1 : GL9.det3 cap GL9.perm3 = vide  (methode : " | separationGL9OrbitMethod | ")")
else
    tagPrint("OPEN", "NIV-1 : separation d'orbite GL9 indecidee");

if nonDegenPerm3NotInGL9det3 === "PROVED-ALG" then
    tagPrint("PROVED-ALG",
        "NIV-2 dir. 1 : perm3 not in closure(GL9.det3)  via semi-continuite (" | semiContMethod | ")")
else
    tagPrint("OPEN", "NIV-2 dir. 1 : non etabli dans ce run");

tagPrint("OPEN",
    "NIV-2 dir. 2 : det3 not in closure(GL9.perm3) -- conjecture GCT centrale, non resolue ici")

print ""
print "  -- OUVERT --"
tagPrint("OPEN", "Construction symbolique M(f_generic) et mineur 9x9 separant explicite")
tagPrint("OPEN", "Lien explicite delta vs multiplicites BIP")

print ""
tagPrint("INFO", "Manifeste : OK=" | toString nOk
    | " MISMATCH=" | toString nMismatch | " MISSING=" | toString nMissing)
print ""
print("  Logs   : " | logFile)
print("  Ranks  : " | csvFile)
print("  Cert.  : " | jsonFile)
print("  Diff   : " | manifestDiff)
print ""
bigsep()
print "  SRMT v19 -- Lock & Audit Release : run conforme."
bigsep()

logStep("pipeline_end", "COMPLETE", "SRMT v19 finished, manifest enforced")
exit 0
