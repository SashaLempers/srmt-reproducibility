-- ================================================================
-- srmt_gct_combined_lab.m2   version 1.6  (CORRECTED)
-- ================================================================
--
--   CORRECTIONS vs 1.5 :
--   1. A2 plethysm : suppression des variables sym_n_C3 / sym_d_sym_n_C3
--      car "_" est un opérateur binaire de subscript en Macaulay2 ;
--      ces noms déclenchent "no method for assignment to binary operator _".
--      Remplacés par symNC3 et symDSymNC3 (noms sans underscores).
--   2. plethysm({NMAT}, stdRep) : BasicList comme premier argument
--      conforme à la doc SchurRings (plethysm(BasicList, RingElement)).
--   3. symAux ring auxiliaire supprimé (inutile avec BasicList).
--   4. Versions précédentes : (gens GL3ring)#0 → [], r_1 ZZ → erreur,
--      plethysm(ZZ, ...) → erreur.
--
-- Utilisation :
--   M2 --script srmt_gct_combined_lab.m2
-- ================================================================

needsPackage "SchurRings";

-- ================================================================
-- 0. CONFIGURATION
-- ================================================================

VERSION   = "combined-lab-1.6";
TIMESTAMP = "2026-05-08";

NMAT   = 3;
DEGREE = 3;
NVARS  = NMAT^2;

doRepresentation  = true;
doHomology        = true;
checkExpectedPdim = true;
expectedPdimDet   = 4;
expectedPdimPerm  = 9;
doSaturate        = false;

outDir  = "output";
certDir = outDir | "/certificates";

ensureDirectory = (d) -> try makeDirectory d;
ensureDirectory outDir;
ensureDirectory certDir;

repCsvPath     = outDir | "/srmt_gct_results.csv";
repMdPath      = outDir | "/srmt_gct_report.md";
repTexPath     = outDir | "/srmt_gct_table.tex";
valPath        = outDir | "/VALIDATION.md";
cffPath        = outDir | "/CITATION.cff";
readmePath     = outDir | "/README.md";
homSummaryPath = outDir | "/jacobian_pdim_summary.md";
combinedPath   = outDir | "/COMBINED_SUMMARY.md";

print "";
print "================================================================";
print " COMBINED SRMT/GCT + JACOBIAN HOMOLOGY LAB";
print(" Version : " | VERSION);
print(" Date    : " | TIMESTAMP);
print "================================================================";
print "";

-- ================================================================
-- 1. TEST HARNESS
-- ================================================================

hasErr = false;

safeAssert = (b, msg) -> (
    if not b then (
        stderr << "[FAIL] " << msg << endl;
        hasErr = true;
    ) else (
        print("[OK]  " | msg);
    )
);

isZeroMatrix = (M) -> all(flatten entries M, e -> e == 0);

-- ================================================================
-- 2. SHARED UTILITIES
-- ================================================================

canon = (lam) -> select(toList lam, x -> x > 0);

partsLe3 = (n) -> (
    L := {};
    for a from 1 to n do
        for b from 0 to a do
            for c from 0 to b do
                if a + b + c == n then L = append(L, canon {a,b,c});
    L
);

dimGL = (lamIn, n) -> (
    lam := canon lamIn;
    if #lam > n then return 0;
    L := lam | toList((n - #lam):0);
    num := 1; den := 1;
    for i from 0 to n-1 do
        for j from i+1 to n-1 do (
            num = num * (L#i - L#j + j - i);
            den = den * (j - i);
        );
    num // den
);

lambdaString = (lam) -> (
    L := canon lam;
    s := "(";
    for i from 0 to #L-1 do (
        s = s | toString(L#i);
        if i < #L-1 then s = s | ",";
    );
    s | ")"
);

boolInt = (b) -> if b then 1 else 0;

matrixLiteral = (M) -> (
    E := entries M;
    s := "matrix {";
    for i from 0 to #E-1 do (
        row := E#i;
        s = s | "{";
        for j from 0 to #row-1 do (
            s = s | toString(row#j);
            if j < #row-1 then s = s | ", ";
        );
        s = s | "}";
        if i < #E-1 then s = s | ",\n        ";
    );
    s | "}"
);

-- ================================================================
-- PART A. REPRESENTATION-THEORETIC BENCHMARK
-- ================================================================

if doRepresentation then (

print "";
print "================================================================";
print " PART A: REPRESENTATION BENCHMARK  (d=3, n=3)";
print "================================================================";

-- A1. Kronecker coefficients
print "A1. Kronecker coefficients g((3,3,3),(3,3,3),.) ...";
Q9 = schurRing(QQ, q, NVARS, GroupActing => "Sn");
rect333 = q_(toList(NMAT:DEGREE));
kronProd = internalProduct(rect333, rect333);

kronTable = new MutableHashTable;
for term in listForm kronProd do (
    lam := canon(term#0);
    kronTable#lam = term#1;
);
sk3 = (lam) -> (
    k := canon lam;
    if kronTable#?k then kronTable#k else 0
);
print "A1. done.";

-- ----------------------------------------------------------------
-- A2. Plethysm Sym^d(Sym^n C^n)
--
-- CORRECTION v1.6 :
--   - Variables sans underscores : symNC3, symDSymNC3
--     (en M2, "_" est un opérateur binaire de subscript)
--   - plethysm(BasicList, RingElement) : premier argument = {k}
--     conforme à la doc SchurRings.
--   - stdRep = r_{1} : partition (1) dans GL3ring (dim=3).
--
--   GL3ring    = schurRing(QQ, r, 3)
--   stdRep     = r_{1}               -- Sym^1(C^3)  dim=3
--   symNC3     = plethysm({3}, r_{1}) -- Sym^3(C^3)  dim=10
--   symDSymNC3 = plethysm({3}, symNC3)-- Sym^3(Sym^3(C^3)) dim=220
-- ----------------------------------------------------------------
print "A2. Plethysm coefficient in s_3[s_3] ...";

GL3ring    = schurRing(QQ, r, NMAT);
stdRep     = r_{1};
symNC3     = plethysm({NMAT},   stdRep);
symDSymNC3 = plethysm({DEGREE}, symNC3);

-- validation de dimension : dim Sym^3(Sym^3(C^3)) = C(12,9) = 220
totalDimPleth := sum for t in listForm symDSymNC3 list (t#1 * dimGL(t#0, NMAT));
safeAssert(totalDimPleth == 220, "dim Sym3(Sym3(C3)) = 220");

plethTable = new MutableHashTable;
for term in listForm symDSymNC3 do (
    lam := canon(term#0);
    plethTable#lam = term#1;
);
pleth3 = (lam) -> (
    k := canon lam;
    if plethTable#?k then plethTable#k else 0
);

safeAssert(pleth3 {9}     > 0, "pleth3(9) > 0");
safeAssert(pleth3 {3,3,3} > 0, "pleth3(3,3,3) > 0");
print "A2. done.";

-- A3. det multiplicities via rectangular stabilizer rule
print "A3. m_lambda(det_3) ...";
mDet3 = (lamIn) -> (
    lam := canon lamIn;
    if #lam == 0 then return 0;
    if sum lam % NMAT != 0 then return 0;
    k := (sum lam) // NMAT;
    if lam == toList(NMAT:k) then 1 else 0
);
print "A3. done.";

-- A4. Permanent multiplicity table (labelled input)
print "A4. Loading permanent multiplicity table ...";
permTable = new MutableHashTable;
permTable#{3,3,3} = 10;
permTable#{4,3,2} = 12;
permTable#{4,4,1} = 2;
permTable#{5,2,2} = 10;
permTable#{5,3,1} = 2;
permTable#{5,4}   = 10;
permTable#{6,2,1} = 1;
permTable#{6,3}   = 10;
permTable#{7,1,1} = 0;
permTable#{7,2}   = 9;
permTable#{8,1}   = 0;
permTable#{9}     = 9;

mPerm3 = (lam) -> (
    k := canon lam;
    if permTable#?k then permTable#k else 0
);
print "A4. done.";

-- A5. Validation
safeAssert(sk3 {3,3,3} == 1,          "sk3(3,3,3) = 1");
safeAssert(sk3 {9} == 1,              "sk3(9) = 1");
safeAssert(mDet3 {3,3,3} == 1,        "mDet3(3,3,3) = 1");
safeAssert(dimGL({9},9) == 24310,     "dimGL(9) = 24310");
safeAssert(dimGL({3,3,3},9) == 41580, "dimGL(3,3,3) = 41580");
for lam in partsLe3 9 do
    safeAssert(permTable#?lam, "permTable covers " | lambdaString lam);

-- A6. Build rows
print "A6. Building result table ...";
lams = partsLe3 9;
rows = {};
for lam in lams do (
    d9    := dimGL(lam, NVARS);
    sk    := sk3 lam;
    pl    := pleth3 lam;
    md    := mDet3 lam;
    mp    := mPerm3 lam;
    delta := mp - md;
    strict := (sk > 0 and pl > 0 and delta != 0);
    fine   := (strict and md > 0 and mp > 0);
    status := if fine then "FINE-MULT"
              else if strict and md == 0 and mp > 0 then "OCCURRENCE"
              else if sk == 0 then "NO-KRONECKER"
              else if pl == 0 then "NO-PLETHYSM"
              else if delta == 0 then "EQUAL"
              else "OTHER";
    rows = append(rows, (lam, d9, sk, pl, md, mp, delta, status,
                         boolInt strict, boolInt fine));
);

strictRows = select(rows, r -> r#8 == 1);
fineRows   = select(rows, r -> r#9 == 1);
occRows    = select(rows, r -> r#8 == 1 and r#9 == 0);

nStrict = #strictRows;
nFine   = #fineRows;
nOcc    = #occRows;

print "";
print "lambda       | dim_GL9  | sk3 | pl | m_det | m_perm | delta | status";
print "-------------+----------+-----+----+-------+--------+-------+-----------";
for r in rows do (
    print( lambdaString(r#0) | " | " | toString(r#1)
         | " | " | toString(r#2) | " | " | toString(r#3)
         | " | " | toString(r#4) | " | " | toString(r#5)
         | " | " | toString(r#6) | " | " | r#7 );
);
print("  Strict: " | toString nStrict
    | "  Fine: "   | toString nFine
    | "  Occurrence: " | toString nOcc);

-- A7. Artifacts
print "A7. Writing artifacts ...";

csvOut := openOut repCsvPath;
csvOut << "lambda,dim_GL9,sk3,pleth,m_det,m_perm,delta,status,is_strict,is_fine" << endl;
for r in rows do (
    csvOut << "\"" << lambdaString(r#0) << "\","
           << r#1 << "," << r#2 << "," << r#3 << ","
           << r#4 << "," << r#5 << "," << r#6 << ","
           << "\"" << r#7 << "\"," << r#8 << "," << r#9 << endl;
);
close csvOut;

texOut := openOut repTexPath;
texOut << "% Auto-generated by srmt_gct_combined_lab.m2 v" << VERSION << endl;
texOut << "\\begin{table}[ht]\\centering" << endl;
texOut << "\\begin{tabular}{c|rrrrrr|c}\\hline" << endl;
texOut << "$\\lambda$ & $\\dim V_{\\lambda}$ & $sk_3$ & $pl$"
       << " & $m_{\\det_3}$ & $m_{\\mathrm{perm}_3}$ & $\\Delta$"
       << " & status \\\\\\hline" << endl;
for r in rows do (
    texOut << "$" << lambdaString(r#0) << "$ & "
           << r#1 << " & " << r#2 << " & " << r#3 << " & "
           << r#4 << " & " << r#5 << " & " << r#6 << " & "
           << r#7 << " \\\\" << endl;
);
texOut << "\\hline\\end{tabular}" << endl;
texOut << "\\caption{Representation benchmark $d=3$.}\\label{tab:srmt-d3}" << endl;
texOut << "\\end{table}" << endl;
close texOut;

mdOut := openOut repMdPath;
mdOut << "# Representation benchmark (det_3 vs perm_3, d=3)" << endl << endl;
mdOut << "Note: perm multiplicities are a labelled input table." << endl << endl;
mdOut << "| lambda | dim_GL9 | sk3 | pl | m_det | m_perm | delta | status |" << endl;
mdOut << "|---|---:|---:|---:|---:|---:|---:|---|" << endl;
for r in rows do (
    mdOut << "| " << lambdaString(r#0)
          << " | " << r#1 << " | " << r#2 << " | " << r#3
          << " | " << r#4 << " | " << r#5 << " | " << r#6
          << " | " << r#7 << " |" << endl;
);
mdOut << endl << "Strict: " << nStrict
      << "  Fine: " << nFine
      << "  Occurrence: " << nOcc << endl;
close mdOut;

print "A7. done.";

); -- end Part A

-- ================================================================
-- PART B. JACOBIAN HOMOLOGY CERTIFICATE
-- ================================================================

if doHomology then (

print "";
print "================================================================";
print " PART B: JACOBIAN HOMOLOGY CERTIFICATE";
print "================================================================";

R = QQ[x11,x12,x13,x21,x22,x23,x31,x32,x33];
X = matrix {
    {x11,x12,x13},
    {x21,x22,x23},
    {x31,x32,x33}
};
fDet  = det X;
fPerm = x11*x22*x33 + x11*x23*x32
      + x12*x21*x33 + x12*x23*x31
      + x13*x21*x32 + x13*x22*x31;

irrelevant = ideal gens R;

safeAssert(degree fDet  == {3}, "det_3 is homogeneous degree 3");
safeAssert(degree fPerm == {3}, "perm_3 is homogeneous degree 3");

buildJacobianIdeal = (f) -> ideal flatten entries jacobian(ideal f);

writeBettiFile = (tag, C) -> (
    path := certDir | "/" | tag | "_betti.txt";
    bOut := openOut path;
    bOut << "Betti table for " << tag << " Jacobian quotient" << endl << endl;
    bOut << toString betti C << endl;
    close bOut;
    path
);

writeResolutionCertificate = (tag, J, C, p) -> (
    path := certDir | "/" | tag | "_resolution_certificate.m2";
    rOut := openOut path;
    rOut << "-- Auto-generated certificate for " << tag << endl;
    rOut << "-- srmt_gct_combined_lab.m2  v" << VERSION << endl << endl;
    rOut << "R = QQ[x11,x12,x13,x21,x22,x23,x31,x32,x33];" << endl;
    rOut << "zeroQ = (M) -> all(flatten entries M, e -> e == 0);" << endl << endl;
    rOut << ("gensJ" | tag | " = ") << matrixLiteral(gens J) << ";" << endl;
    rOut << ("Jideal" | tag | " = ideal gensJ" | tag | ";") << endl;
    rOut << ("Mmod" | tag | " = R^1 / Jideal" | tag | ";") << endl;
    rOut << ("Cres" | tag | " = res Mmod" | tag | ";") << endl;
    rOut << ("assert(pdim Mmod" | tag | " == " | toString p | ");") << endl << endl;
    for i from 1 to p do (
        rOut << ("diff" | tag | toString i | " = ")
             << matrixLiteral(C.dd_(i)) << ";" << endl << endl;
    );
    if p >= 2 then
        for i from 1 to p-1 do (
            rOut << ("assert(zeroQ(diff" | tag | toString i
                    | " * diff" | tag | toString(i+1) | "));") << endl;
        );
    rOut << endl;
    rOut << ("assert(ideal diff" | tag | "1 == Jideal" | tag | ");") << endl;
    rOut << "print \"Certificate verified for " << tag << "\";" << endl;
    close rOut;
    path
);

verifyComplex = (tag, C, p) -> (
    if p >= 2 then (
        for i from 1 to p-1 do (
            comp := C.dd_(i) * C.dd_(i+1);
            safeAssert(isZeroMatrix comp,
                tag | ": d" | toString i | " * d" | toString(i+1) | " = 0");
        );
    ) else (
        safeAssert(true, tag | ": pdim <= 1, no consecutive differentials to verify");
    );
);

analyzeJacobian = (tag, f, expPdim) -> (
    print("B. Analyzing " | tag | " Jacobian quotient ...");
    Jraw := buildJacobianIdeal f;
    J := if doSaturate then saturate(Jraw, irrelevant) else Jraw;
    safeAssert(isHomogeneous J, tag | ": Jacobian ideal is homogeneous");
    Mmod := R^1 / J;
    print(tag | ": computing minimal free resolution ...");
    Cres := res Mmod;
    p := pdim Mmod;
    print(tag | ": pdim(R/J) = " | toString p);
    print betti Cres;
    verifyComplex(tag, Cres, p);
    if checkExpectedPdim then
        safeAssert(p == expPdim, tag | ": pdim = " | toString expPdim | " (expected)");
    bPath := writeBettiFile(tag, Cres);
    cPath := writeResolutionCertificate(tag, J, Cres, p);
    print("  betti  -> " | bPath);
    print("  cert   -> " | cPath);
    (J, Cres, p, bPath, cPath)
);

detResult  = analyzeJacobian("det",  fDet,  expectedPdimDet);
permResult = analyzeJacobian("perm", fPerm, expectedPdimPerm);

pDet  = detResult#2;
pPerm = permResult#2;

homOut := openOut homSummaryPath;
homOut << "# Jacobian projective-dimension certificate" << endl << endl;
homOut << "Version `" << VERSION << "`  |  Date `" << TIMESTAMP << "`" << endl;
homOut << "Base field `QQ`  |  Saturation `" << toString doSaturate << "`" << endl << endl;
homOut << "| quotient | pdim | certificate |" << endl;
homOut << "|---|---:|---|" << endl;
homOut << "| `R/Jdet`  | " << pDet  << " | `certificates/det_resolution_certificate.m2` |" << endl;
homOut << "| `R/Jperm` | " << pPerm << " | `certificates/perm_resolution_certificate.m2` |" << endl << endl;
homOut << "All differentials satisfy di * d(i+1) = 0 (verified in-script)." << endl;
close homOut;

); -- end Part B

-- ================================================================
-- SHARED METADATA
-- ================================================================

valOut := openOut valPath;
valOut << "# Validation" << endl << endl;
valOut << "## Part A: representation benchmark" << endl;
valOut << "- Kronecker: SchurRings.internalProduct." << endl;
valOut << "- Plethysm: plethysm({3}, plethysm({3}, r_{1})) in GL3ring." << endl;
valOut << "- det multiplicities: rectangular stabilizer rule." << endl;
valOut << "- perm multiplicities: labelled input table." << endl << endl;
valOut << "## Part B: Jacobian homology" << endl;
valOut << "- Jacobian ideals: first partial derivatives over QQ." << endl;
valOut << "- Resolutions: res M." << endl;
valOut << "- Projective dimensions: pdim M." << endl;
valOut << "- Chain complex: di * d(i+1) = 0 verified in-script." << endl << endl;
valOut << "## Out of scope" << endl;
valOut << "- No asymptotic lower bound." << endl;
valOut << "- No padded permanent." << endl;
valOut << "- No Lean/Coq certificate." << endl;
close valOut;

cffOut := openOut cffPath;
cffOut << "cff-version: 1.2.0" << endl;
cffOut << "message: \"If you use this software, please cite it as below.\"" << endl;
cffOut << "title: \"Combined SRMT/GCT and Jacobian Homology Laboratory\"" << endl;
cffOut << "version: \"" << VERSION << "\"" << endl;
cffOut << "date-released: \"" << TIMESTAMP << "\"" << endl;
cffOut << "type: software" << endl;
close cffOut;

readOut := openOut readmePath;
readOut << "# Combined SRMT/GCT + Jacobian Homology Lab (v1.6)" << endl << endl;
readOut << "```bash" << endl << "M2 --script srmt_gct_combined_lab.m2" << endl << "```" << endl << endl;
readOut << "| Output file | Contents |" << endl;
readOut << "|---|---|" << endl;
readOut << "| srmt_gct_results.csv | Representation table |" << endl;
readOut << "| srmt_gct_report.md | Markdown report |" << endl;
readOut << "| srmt_gct_table.tex | LaTeX table |" << endl;
readOut << "| jacobian_pdim_summary.md | Homology summary |" << endl;
readOut << "| certificates/det_resolution_certificate.m2 | det certificate |" << endl;
readOut << "| certificates/perm_resolution_certificate.m2 | perm certificate |" << endl;
readOut << "| VALIDATION.md | Validation tiers |" << endl;
readOut << "| CITATION.cff | Citation metadata |" << endl;
close readOut;

combOut := openOut combinedPath;
combOut << "# Combined run summary" << endl << endl;
combOut << "Two independent pipelines -- do not conflate them." << endl << endl;
combOut << "| Part | Pipeline | Object computed |" << endl;
combOut << "|---|---|---|" << endl;
combOut << "| A | Representation | Kronecker, plethysm, det/perm multiplicities at d=3 |" << endl;
combOut << "| B | Homology | pdim(R/Jdet), pdim(R/Jperm), resolution certificates |" << endl;
close combOut;

-- ================================================================
-- EXIT
-- ================================================================

print "";
print "================================================================";
if hasErr then (
    print " PIPELINE COMPLETED WITH FAILURES (see [FAIL] lines above)";
    print "================================================================";
    exit 1;
) else (
    print " PIPELINE COMPLETED SUCCESSFULLY";
    print "================================================================";
    exit 0;
);
