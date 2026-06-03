-- ================================================================
-- SRMT v19 — INDEPENDENT VERIFIER
-- ================================================================
--
-- Mission : recalculer INDEPENDAMMENT les invariants GL9 critiques de
-- la section 15 du pipeline v19, lire le certificat JSON produit par
-- srmt_v19_pipeline.m2, et confirmer la correspondance.
--
-- Ce script ne rejoue PAS l'integralite du pipeline. En particulier
-- il ne recalcule pas :
--   - les projecteurs piMinus / piPlus,
--   - le catalecticant,
--   - delta(det3), delta(perm3),
--   - les stress-tests OBS-N,
--   - la partie modulaire.
-- Ces calculs sont locaux, peu sensibles, et leur duplication
-- n'apporterait rien d'exploitable. Seuls les invariants algebriques
-- structurels de §15 sont recalcules.
--
-- Quantites verifiees :
--   - dim Sing V(det3),  dim Sing V(perm3)
--   - dim V(J(det3)),    dim V(J(perm3))
--   - codim J(det3),     codim J(perm3)
--   - degree J(det3),    degree J(perm3)
--   - H_apolar(det3),    H_apolar(perm3)
--   - identite adj(Y) * Y = det3(Y) * I  (recalcul independant)
--
-- Mode d'emploi :
--   M2 --script srmt_v19_verify.m2
--
-- Codes de sortie :
--   0  : tout OK
--   2  : valeur(s) recomputee(s) divergente(s) du certificat
--   3  : valeur(s) recomputee(s) divergente(s) du manifeste attendu
--   4  : certificat introuvable ou parsing impossible
--   5  : identite adj*Y echoue
-- ================================================================

CERT_PATH   = "./srmt_v19_lock_audit_certificate.json"
REPORT_PATH = "./srmt_v19_verify_report.txt"

-- ----------------------------------------------------------------
-- 0. Valeurs attendues (gelees, identiques au manifeste du pipeline)
-- ----------------------------------------------------------------

-- Les chaines attendues correspondent au format JSON produit par le
-- pipeline (jsonValueOf). Les listes sont rendues comme [1, 9, 9, 1].
-- La comparaison est faite apres normalisation par normStr() qui retire
-- les espaces, donc "{1,9,9,1}" et "{1, 9, 9, 1}" sont equivalents.
EXPECTED = {
    {"dim_Sing_det3_QQ",    "5"},
    {"dim_Sing_perm3_QQ",   "3"},
    {"dim_VJ_det3_QQ",      "5"},
    {"dim_VJ_perm3_QQ",     "3"},
    {"codim_J_det3_QQ",     "4"},
    {"codim_J_perm3_QQ",    "6"},
    {"degree_J_det3_QQ",    "6"},
    {"degree_J_perm3_QQ",  "24"},
    {"H_apolar_det3_QQ",    "[1,9,9,1]"},
    {"H_apolar_perm3_QQ",   "[1,9,9,1]"}
}

-- ----------------------------------------------------------------
-- 1. Lecture brute du certificat JSON
-- ----------------------------------------------------------------
-- Strategie volontairement minimaliste : on ne parse pas le JSON
-- comme un arbre. Pour chaque clef du manifeste, on cherche la
-- valeur de "observed" associee a cette clef en exploitant la
-- mise en forme stable du pipeline (une ligne par claim, ordre
-- key -> expected -> observed -> status).

certText = ""
certPresent = false
try (
    certText = get CERT_PATH;
    certPresent = true
) else (
    certPresent = false
);

if not certPresent or #certText == 0 then (
    print("[VERIFY] ERREUR : impossible de lire " | CERT_PATH);
    exit 4
);

-- Extraction tres simple : pour chaque ligne du certificat, on
-- recupere les triplets ("key": "...", "expected": ..., "observed": ...).
-- On stocke dans certObserved#key = chaine "observed" telle quelle.

certObserved = new MutableHashTable
certExpected = new MutableHashTable
certStatus   = new MutableHashTable

-- Extraction d'un champ "name": <value> ou <value> peut etre :
--   - un nombre : 9, -1
--   - un booleen : true, false
--   - une chaine : "..."
--   - une liste : [1, 9, 9, 1]
-- L'algorithme consomme la valeur jusqu'au prochain delimiteur de
-- niveau zero (',' ou '}') en respectant les paires [] et "".
extractField = (line, fname) -> (
    needle := "\"" | fname | "\":";
    pos := regex(needle, line);
    if pos === null then ""
    else (
        startIdx := pos#0#0 + pos#0#1;  -- juste apres "name":
        n := #line;
        -- skip whitespace
        while startIdx < n and (line#startIdx === " " or line#startIdx === "\t") do
            startIdx = startIdx + 1;
        -- collecter jusqu'au delimiteur de niveau zero ',' ou '}'
        depthBracket := 0;
        inString := false;
        i := startIdx;
        done := false;
        while i < n and not done do (
            c := line#i;
            if inString then (
                if c === "\"" then inString = false;
                i = i + 1
            ) else (
                if c === "\"" then (inString = true; i = i + 1)
                else if c === "[" then (depthBracket = depthBracket + 1; i = i + 1)
                else if c === "]" then (depthBracket = depthBracket - 1; i = i + 1)
                else if depthBracket == 0 and (c === "," or c === "}") then done = true
                else i = i + 1
            )
        );
        -- trim trailing whitespace
        j := i;
        while j > startIdx and (line#(j-1) === " " or line#(j-1) === "\t") do
            j = j - 1;
        raw := substring(line, startIdx, j - startIdx);
        -- si c'est une chaine "...", retirer les guillemets exterieurs
        if #raw >= 2 and raw#0 === "\"" and raw#(#raw - 1) === "\"" then
            substring(raw, 1, #raw - 2)
        else raw
    )
)

certLines = lines certText
scan(certLines, ln -> (
    if match("\"key\":", ln) then (
        keyStr := extractField(ln, "key");
        if keyStr =!= "" then (
            certExpected#keyStr = extractField(ln, "expected");
            certObserved#keyStr = extractField(ln, "observed");
            certStatus#keyStr   = extractField(ln, "status")
        )
    )
))

print("[VERIFY] Certificat charge : " | toString (#keys certObserved) | " claim(s)")

-- ----------------------------------------------------------------
-- 2. Recalcul independant des invariants §15 dans un anneau frais
-- ----------------------------------------------------------------

R9V = QQ[z_0..z_8, MonomialOrder => GRevLex]
matIdx = (i, j) -> 3*i + j
YV = matrix apply(3, i -> apply(3, j -> (gens R9V)#(matIdx(i,j))))
detV  = det YV
permV = sum(permutations 3, sigma -> product(3, i -> YV_(i, sigma#i)))

JdetV  = ideal jacobian matrix{{detV}}
JpermV = ideal jacobian matrix{{permV}}
SdetV  = ideal detV  + JdetV
SpermV = ideal permV + JpermV

dimSingDetV  = dim SdetV
dimSingPermV = dim SpermV
dimVJDetV    = dim JdetV
dimVJPermV   = dim JpermV
codimJDetV   = codim JdetV
codimJPermV  = codim JpermV
degreeJDetV  = degree JdetV
degreeJPermV = degree JpermV

-- Hilbert apolaire (recalcul independant, meme algorithme que v18/v19)
apolarHilbertOf = (f, dTot, kmax) -> (
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
)

HDetV  = apolarHilbertOf(detV,  3, 3)
HPermV = apolarHilbertOf(permV, 3, 3)

recomputed = hashTable {
    "dim_Sing_det3_QQ"   => toString dimSingDetV,
    "dim_Sing_perm3_QQ"  => toString dimSingPermV,
    "dim_VJ_det3_QQ"     => toString dimVJDetV,
    "dim_VJ_perm3_QQ"    => toString dimVJPermV,
    "codim_J_det3_QQ"    => toString codimJDetV,
    "codim_J_perm3_QQ"   => toString codimJPermV,
    "degree_J_det3_QQ"   => toString degreeJDetV,
    "degree_J_perm3_QQ"  => toString degreeJPermV,
    "H_apolar_det3_QQ"   => toString HDetV,
    "H_apolar_perm3_QQ"  => toString HPermV
}

-- ----------------------------------------------------------------
-- 3. Comparaison : recompute vs certificat vs manifeste attendu
-- ----------------------------------------------------------------

-- Normalisation de chaines pour comparaison robuste :
--   - retire espaces et tabulations
--   - homogeneise les delimiteurs de liste : { } -> [ ]
normStr = (s) -> (
    t := replace("[ \t]+", "", toString s);
    t = replace("\\{", "[", t);
    t = replace("\\}", "]", t);
    t
)

o = openOut REPORT_PATH
o << "# SRMT v19 -- independent verifier report" << endl
o << "# generated at " << toString currentTime() << endl
o << "# certificate : " << CERT_PATH << endl
o << endl

nMismatchCert     = 0
nMismatchExpected = 0
nMissing          = 0

scan(EXPECTED, e -> (
    key := e#0;
    expectedManifest := e#1;
    rec := if recomputed#?key then recomputed#key else "";
    cer := if certObserved#?key then certObserved#key else "";

    okRecVsExp  := normStr(rec) == normStr(expectedManifest);
    okCertVsExp := normStr(cer) == normStr(expectedManifest);
    okRecVsCert := normStr(rec) == normStr(cer);

    statusStr :=
        if cer === "" then "MISSING_IN_CERT"
        else if okRecVsExp and okCertVsExp and okRecVsCert then "OK"
        else if not okRecVsExp then "MISMATCH_RECOMPUTE_VS_EXPECTED"
        else if not okCertVsExp then "MISMATCH_CERT_VS_EXPECTED"
        else "MISMATCH_RECOMPUTE_VS_CERT";

    if statusStr === "MISSING_IN_CERT" then nMissing = nMissing + 1
    else if statusStr === "MISMATCH_RECOMPUTE_VS_EXPECTED" then nMismatchExpected = nMismatchExpected + 1
    else if statusStr =!= "OK" then nMismatchCert = nMismatchCert + 1;

    o << statusStr << "  " << key << endl
      << "    expected (manifest) = " << expectedManifest << endl
      << "    observed (cert)     = " << cer << endl
      << "    recomputed          = " << rec << endl
      << endl
))

-- ----------------------------------------------------------------
-- 4. Identite adj(Y)*Y = det(Y)*I (recalcul independant)
-- ----------------------------------------------------------------

adjV = matrix apply(3, i -> apply(3, j ->
    (-1)^(i+j) * det submatrix'(YV, {j}, {i})))
idCheckV = (adjV * YV) - detV * id_(R9V^3)

adjOk = (idCheckV == 0)

if adjOk then (
    o << "OK  adj_identity_QQ_recomputed  : adj(Y)*Y = det3(Y)*I sur QQ[z]" << endl
) else (
    o << "MISMATCH adj_identity_QQ_recomputed" << endl;
    nMismatchExpected = nMismatchExpected + 1
);

-- ----------------------------------------------------------------
-- 5. Verification de la coherence claims_status du certificat
-- ----------------------------------------------------------------

o << endl;
o << "# claims_status (lu depuis le certificat)" << endl
nCertOk       = 0
nCertMismatch = 0
nCertMissing  = 0
scan(keys certStatus, k -> (
    s := certStatus#k;
    if s === "OK" then nCertOk = nCertOk + 1
    else if s === "MISMATCH" then nCertMismatch = nCertMismatch + 1
    else if s === "MISSING" then nCertMissing = nCertMissing + 1
));
o << "  cert OK="        << toString nCertOk
  << " MISMATCH="        << toString nCertMismatch
  << " MISSING="         << toString nCertMissing << endl

-- ----------------------------------------------------------------
-- 6. Resume final
-- ----------------------------------------------------------------

o << endl
o << "# summary" << endl
o << "  recompute_vs_expected_mismatch = " << toString nMismatchExpected << endl
o << "  cert_vs_expected_mismatch      = " << toString nMismatchCert     << endl
o << "  missing_in_cert                = " << toString nMissing          << endl
o << "  adj_identity_recomputed_ok     = " << toString adjOk             << endl
close o

print("[VERIFY] Rapport ecrit : " | REPORT_PATH)
print("[VERIFY] recompute_vs_expected_mismatch = " | toString nMismatchExpected)
print("[VERIFY] cert_vs_expected_mismatch      = " | toString nMismatchCert)
print("[VERIFY] missing_in_cert                = " | toString nMissing)
print("[VERIFY] adj_identity_recomputed_ok     = " | toString adjOk)

if not adjOk then (
    print "[VERIFY] FAIL : identite adj*Y echoue";
    exit 5
);

if nMismatchExpected > 0 then (
    print "[VERIFY] FAIL : recalcul independant differe du manifeste attendu";
    exit 3
);

if nMismatchCert > 0 or nMissing > 0 then (
    print "[VERIFY] FAIL : certificat divergent du recalcul ou claim manquant";
    exit 2
);

print "[VERIFY] OK -- certificat conforme au recalcul independant et au manifeste"
exit 0
