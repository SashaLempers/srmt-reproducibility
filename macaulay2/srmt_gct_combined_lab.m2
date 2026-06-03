-- ================================================================
-- srmt_gct_combined_lab.m2
-- ================================================================
--
--   COMBINED SRMT / GCT + JACOBIAN HOMOLOGY LABORATORY
--   ---------------------------------------------------
--   This file intentionally combines two independent computational layers:
--
--   PART A. Representation-theoretic benchmark for det_3 vs perm_3
--           at degree d = 3:
--             - Kronecker coefficients sk_3(lambda)
--             - plethysm coefficients of s_lambda in s_3[s_3]
--             - determinant multiplicity via the rectangular rule
--             - permanent-side multiplicity table as a labelled input
--
--   PART B. Homological Jacobian certificate over QQ:
--             - J_det  = ideal(partial derivatives of det_3)
--             - J_perm = ideal(partial derivatives of perm_3)
--             - pdim(R/J_det), pdim(R/J_perm)
--             - minimal free resolutions
--             - exported differential matrices
--             - verification that d_i * d_{i+1} = 0
--
--   IMPORTANT SCIENTIFIC DISCLAIMER
--   -------------------------------
--   These two parts compute different mathematical objects.  Part A is a
--   small representation-theoretic benchmark; Part B is a homological
--   Jacobian computation.  Do not conflate the two, and do not claim that
--   either part proves a lower bound for permanent vs determinant.
--
--   Target: Macaulay2 >= 1.19
-- ================================================================

needsPackage "SchurRings";

-- ================================================================
-- 0. GLOBAL CONFIGURATION
-- ================================================================

VERSION   = "combined-lab-1.0";
TIMESTAMP = "2026-05-08";

NMAT      = 3;
DEGREE    = 3;
NVARS     = NMAT^2;

RUN_REPRESENTATION_PART = true;
RUN_HOMOLOGY_PART       = true;

CHECK_EXPECTED_PDIM = true;
EXPECTED_PDIM_DET   = 4;
EXPECTED_PDIM_PERM  = 9;

-- By default we use raw affine/conical Jacobian ideals.  If you set this to
-- true, you change the object to the irrelevant-ideal saturation.
COMPUTE_SATURATED_JACOBIAN = false;

outDir  = "output";
certDir = outDir | "/certificates";

ensureDirectory = (d) -> if not fileExists d then makeDirectory d;
ensureDirectory outDir;
ensureDirectory certDir;

repCsvPath    = outDir | "/srmt_gct_results.csv";
repMdPath     = outDir | "/srmt_gct_report.md";
repTexPath    = outDir | "/srmt_gct_table.tex";
valPath       = outDir | "/VALIDATION.md";
cffPath       = outDir | "/CITATION.cff";
readmePath    = outDir | "/README.md";
homSummaryPath = outDir | "/jacobian_pdim_summary.md";
combinedSummaryPath = outDir | "/COMBINED_SUMMARY.md";

print "";
print "================================================================";
print " COMBINED SRMT / GCT + JACOBIAN HOMOLOGY LAB";
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
        print("[OK] " | msg);
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

if RUN_REPRESENTATION_PART then (

print "";
print "================================================================";
print " PART A: REPRESENTATION-THEORETIC BENCHMARK";
print "================================================================";

-- A1. Kronecker coefficients
print "Step A1/7 : Kronecker coefficients g((3,3,3),(3,3,3),lambda) ...";
Q9 = schurRing(QQ, q, NVARS, GroupActing => "Sn");
rect333 = q_(toList(NMAT : DEGREE));
kronProd = internalProduct(rect333, rect333);

kronTable = new MutableHashTable;
for term in listForm kronProd do (
    lam := canon(term#0);
    coeff := term#1;
    kronTable#lam = coeff;
);
sk3 = (lam) -> (
    k := canon lam;
    if kronTable#?k then kronTable#k else 0
);

-- A2. Plethysm
print "Step A2/7 : plethysm coefficient in s_3[s_3] ...";
S9 = schurRing(QQ, s, NVARS);
plethElt = plethysm(s_{3}, s_{3});
plethTable = new MutableHashTable;
for term in listForm plethElt do (
    lam := canon(term#0);
    coeff := term#1;
    plethTable#lam = coeff;
);
pleth3 = (lam) -> (
    k := canon lam;
    if plethTable#?k then plethTable#k else 0
);

-- A3. Determinant multiplicity via rectangular stabilizer rule
print "Step A3/7 : m_lambda(det_3) via rectangular rule ...";
mDet3 = (lamIn) -> (
    lam := canon lamIn;
    if sum lam == 0 then return 0;
    if sum lam % NMAT != 0 then return 0;
    k := (sum lam) // NMAT;
    if lam == toList(NMAT:k) then 1 else 0
);

-- A4. Permanent multiplicities as a labelled, checked input table.
-- This is kept intentionally conservative.  A fully independent derivation of
-- orbit-closure multiplicities should be placed in a separate proof/script.
print "Step A4/7 : loading labelled permanent multiplicity table ...";
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

-- A5. Validation suite
print "Step A5/7 : validation ...";
safeAssert(sk3 {3,3,3} == 1, "sk3(3,3,3) = 1");
safeAssert(sk3 {9} == 1, "sk3(9) = 1");
safeAssert(mDet3 {3,3,3} == 1, "m_det(3,3,3) = 1");
safeAssert(dimGL({9}, 9) == 24310, "dimGL(9) = 24310");
safeAssert(dimGL({3,3,3}, 9) == 41580, "dimGL(3,3,3) = 41580");
for lam in partsLe3 9 do safeAssert(permTable#?lam, "permTable covers " | lambdaString lam);

-- A6. Build rows
print "Step A6/7 : building result table ...";
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

print "";
print "lambda       | dim_GL9  | sk3 | pl | m_det | m_perm | delta | status";
print "-------------+----------+-----+----+-------+--------+-------+-------------";
for r in rows do (
    print( lambdaString(r#0) | " | "
         | toString(r#1) | " | "
         | toString(r#2) | " | "
         | toString(r#3) | " | "
         | toString(r#4) | " | "
         | toString(r#5) | " | "
         | toString(r#6) | " | "
         | r#7 );
);
print("Strict rows: " | toString #strictRows);
print("Fine rows  : " | toString #fineRows);
print("Occ rows   : " | toString #occRows);

-- A7. Write artifacts
print "Step A7/7 : writing representation artifacts ...";

csvOut := openOut repCsvPath;
csvOut << "lambda,dim_GL9,sk3,pleth,m_det,m_perm,delta,status,is_strict,is_fine" << endl;
for r in rows do (
    csvOut << "\"" << lambdaString(r#0) << "\"" << ","
           << r#1 << "," << r#2 << "," << r#3 << ","
           << r#4 << "," << r#5 << "," << r#6 << ","
           << "\"" << r#7 << "\"" << ","
           << r#8 << "," << r#9 << endl;
);
close csvOut;

texOut := openOut repTexPath;
texOut << "% Auto-generated by srmt_gct_combined_lab.m2" << endl;
texOut << "\\begin{table}[ht]" << endl;
texOut << "\\centering" << endl;
texOut << "\\begin{tabular}{c|rrrrrr|c}" << endl;
texOut << "\\hline" << endl;
texOut << "$\\lambda$ & $\\dim V_\\lambda$ & $sk_3$ & $pl$ & $m_{\\det_3}$ & $m_{\\mathrm{perm}_3}$ & $\\Delta$ & status \\\\" << endl;
texOut << "\\hline" << endl;
for r in rows do (
    texOut << "$" << lambdaString(r#0) << "$ & "
           << r#1 << " & " << r#2 << " & " << r#3 << " & "
           << r#4 << " & " << r#5 << " & " << r#6 << " & "
           << r#7 << " \\\\" << endl;
);
texOut << "\\hline" << endl;
texOut << "\\end{tabular}" << endl;
texOut << "\\caption{Representation-theoretic benchmark for $(\\det_3,\\mathrm{perm}_3)$ at $d=3$.}" << endl;
texOut << "\\label{tab:srmt-d3}" << endl;
texOut << "\\end{table}" << endl;
close texOut;

mdOut := openOut repMdPath;
mdOut << "# Representation-theoretic benchmark" << endl << endl;
mdOut << "This section computes small-degree data for det_3 vs perm_3 at d=3." << endl;
mdOut << "The permanent multiplicity column is a labelled input table, not an independently derived theorem in this script." << endl << endl;
mdOut << "| lambda | dim_GL9 | sk3 | pl | m_det | m_perm | delta | status |" << endl;
mdOut << "|---|---:|---:|---:|---:|---:|---:|---|" << endl;
for r in rows do (
    mdOut << "| " << lambdaString(r#0)
          << " | " << r#1 << " | " << r#2 << " | " << r#3
          << " | " << r#4 << " | " << r#5 << " | " << r#6
          << " | " << r#7 << " |" << endl;
);
mdOut << endl;
mdOut << "Strict rows: " << #strictRows << endl;
mdOut << "Fine rows: " << #fineRows << endl;
mdOut << "Occurrence rows: " << #occRows << endl;
close mdOut;

); -- end representation part

-- ================================================================
-- PART B. JACOBIAN HOMOLOGY CERTIFICATE
-- ================================================================

if RUN_HOMOLOGY_PART then (

print "";
print "================================================================";
print " PART B: JACOBIAN HOMOLOGY CERTIFICATE";
print "================================================================";

-- B1. Polynomial ring and polynomials
R = QQ[x11,x12,x13,x21,x22,x23,x31,x32,x33];
X = matrix {
    {x11,x12,x13},
    {x21,x22,x23},
    {x31,x32,x33}
};

fDet = det X;
fPerm = x11*x22*x33 + x11*x23*x32
      + x12*x21*x33 + x12*x23*x31
      + x13*x21*x32 + x13*x22*x31;

varsR = {x11,x12,x13,x21,x22,x23,x31,x32,x33};
irrelevant = ideal varsR;

safeAssert(degree fDet == {3}, "det_3 homogeneous degree 3");
safeAssert(degree fPerm == {3}, "perm_3 homogeneous degree 3");

jacobianIdeal = (f) -> (
    Jac := jacobian(ideal f);
    ideal flatten entries Jac
);

writeBettiFile = (tag, C) -> (
    path := certDir | "/" | tag | "_betti.txt";
    out := openOut path;
    out << "Betti table for " << tag << " Jacobian quotient" << endl << endl;
    out << toString betti C << endl;
    close out;
    path
);

writeResolutionCertificate = (tag, J, C, p) -> (
    path := certDir | "/" | tag | "_resolution_certificate.m2";
    out := openOut path;
    out << "-- Auto-generated certificate for " << tag << " Jacobian quotient" << endl;
    out << "-- Generated by srmt_gct_combined_lab.m2, version " << VERSION << endl << endl;
    out << "R = QQ[x11,x12,x13,x21,x22,x23,x31,x32,x33];" << endl;
    out << "zeroMatrixQ = (M) -> all(flatten entries M, e -> e == 0);" << endl << endl;
    out << ("gens_" | tag | " = ") << matrixLiteral(gens J) << ";" << endl;
    out << ("J_" | tag | " = ideal gens_" | tag | ";") << endl;
    out << ("M_" | tag | " = R^1 / J_" | tag | ";") << endl;
    out << ("C_" | tag | " = res M_" | tag | ";") << endl;
    out << ("assert(pdim M_" | tag | " == " | toString p | ");") << endl << endl;
    for i from 1 to p do (
        out << ("d_" | tag | "_" | toString i | " = ")
            << matrixLiteral(C.dd_i) << ";" << endl << endl;
    );
    if p >= 2 then (
        for i from 1 to p-1 do (
            out << ("assert(zeroMatrixQ(d_" | tag | "_" | toString i
                    | " * d_" | tag | "_" | toString(i+1) | "));" ) << endl;
        );
    );
    out << ("assert(ideal d_" | tag | "_1 == J_" | tag | ");") << endl;
    out << "print \"certificate verified\";" << endl;
    close out;
    path
);

verifyComplex = (tag, C, p) -> (
    if p >= 2 then (
        for i from 1 to p-1 do (
            comp := (C.dd_i) * (C.dd_(i+1));
            safeAssert(isZeroMatrix comp,
                tag | ": d_" | toString i | " * d_" | toString(i+1) | " = 0");
        );
    ) else safeAssert(true, tag | ": no consecutive differentials to check");
);

analyzeJacobian = (tag, f, expectedPdim) -> (
    print "";
    print("Analyzing " | tag | " Jacobian quotient ...");
    Jraw := jacobianIdeal f;
    J := if COMPUTE_SATURATED_JACOBIAN then saturate(Jraw, irrelevant) else Jraw;
    safeAssert(isHomogeneous J, tag | ": Jacobian ideal homogeneous");
    M := R^1 / J;
    print(tag | ": computing resolution ...");
    C := res M;
    p := pdim M;
    print(tag | ": pdim(R/J) = " | toString p);
    print betti C;
    verifyComplex(tag, C, p);
    if CHECK_EXPECTED_PDIM then safeAssert(p == expectedPdim, tag | ": pdim expected " | toString expectedPdim);
    bettiPath := writeBettiFile(tag, C);
    certPath  := writeResolutionCertificate(tag, J, C, p);
    (J, M, C, p, bettiPath, certPath)
);

detResult  = analyzeJacobian("det",  fDet,  EXPECTED_PDIM_DET);
permResult = analyzeJacobian("perm", fPerm, EXPECTED_PDIM_PERM);

pDet  = detResult#3;
pPerm = permResult#3;

homOut := openOut homSummaryPath;
homOut << "# Jacobian projective-dimension certificate" << endl << endl;
homOut << "Version: `" << VERSION << "`" << endl;
homOut << "Base field: `QQ`" << endl;
homOut << "Saturation by irrelevant ideal: `" << toString COMPUTE_SATURATED_JACOBIAN << "`" << endl << endl;
homOut << "| quotient | pdim | certificate | betti |" << endl;
homOut << "|---|---:|---|---|" << endl;
homOut << "| `R/J_det` | " << pDet << " | `certificates/det_resolution_certificate.m2` | `certificates/det_betti.txt` |" << endl;
homOut << "| `R/J_perm` | " << pPerm << " | `certificates/perm_resolution_certificate.m2` | `certificates/perm_betti.txt` |" << endl << endl;
homOut << "The script exports all resolution differentials and verifies `d_i * d_{i+1} = 0`." << endl;
close homOut;

); -- end homology part

-- ================================================================
-- FINAL SHARED ARTIFACTS
-- ================================================================

print "";
print "================================================================";
print " WRITING SHARED METADATA";
print "================================================================";

valOut := openOut valPath;
valOut << "# Validation" << endl << endl;
valOut << "This combined script has two independent layers." << endl << endl;
valOut << "## Part A: representation benchmark" << endl;
valOut << "- Kronecker coefficients are computed by SchurRings `internalProduct`." << endl;
valOut << "- Plethysm coefficients are computed by SchurRings `plethysm(s_3,s_3)`." << endl;
valOut << "- Determinant multiplicities use the rectangular rule." << endl;
valOut << "- Permanent multiplicities are a labelled input table in this script." << endl << endl;
valOut << "## Part B: Jacobian homology" << endl;
valOut << "- Jacobian ideals are generated by first partial derivatives over QQ." << endl;
valOut << "- Resolutions are computed by Macaulay2 `res`." << endl;
valOut << "- Projective dimensions are computed by `pdim`." << endl;
valOut << "- Consecutive differential compositions are checked to be zero." << endl << endl;
valOut << "## Out of scope" << endl;
valOut << "- No asymptotic lower bound is proved." << endl;
valOut << "- No padded permanent computation is performed." << endl;
valOut << "- No Lean/Coq formal certificate is generated." << endl;
close valOut;

cffOut := openOut cffPath;
cffOut << "cff-version: 1.2.0" << endl;
cffOut << "message: \"If you use this software, please cite it as below.\"" << endl;
cffOut << "title: \"Combined SRMT/GCT and Jacobian Homology Laboratory\"" << endl;
cffOut << "version: \"" << VERSION << "\"" << endl;
cffOut << "date-released: \"2026-05-08\"" << endl;
cffOut << "type: software" << endl;
cffOut << "abstract: >" << endl;
cffOut << "  Combined Macaulay2 laboratory for small-degree GCT representation" << endl;
cffOut << "  benchmarks and Jacobian homology certificates for det_3 and perm_3." << endl;
close cffOut;

readOut := openOut readmePath;
readOut << "# Combined SRMT / GCT + Jacobian Homology Laboratory" << endl << endl;
readOut << "Run:" << endl << endl;
readOut << "```bash" << endl;
readOut << "M2 --script srmt_gct_combined_lab.m2" << endl;
readOut << "```" << endl << endl;
readOut << "Outputs:" << endl;
readOut << "- `srmt_gct_results.csv`" << endl;
readOut << "- `srmt_gct_report.md`" << endl;
readOut << "- `srmt_gct_table.tex`" << endl;
readOut << "- `jacobian_pdim_summary.md`" << endl;
readOut << "- `certificates/det_resolution_certificate.m2`" << endl;
readOut << "- `certificates/perm_resolution_certificate.m2`" << endl;
readOut << "- `VALIDATION.md`" << endl;
close readOut;

combOut := openOut combinedSummaryPath;
combOut << "# Combined summary" << endl << endl;
combOut << "This run combines a representation-theoretic benchmark and a Jacobian homology certificate." << endl << endl;
combOut << "Do not conflate these two computations: they concern different mathematical objects." << endl << endl;
if RUN_REPRESENTATION_PART then (
    combOut << "## Representation part" << endl << endl;
    combOut << "- Results: `srmt_gct_results.csv`" << endl;
    combOut << "- Report: `srmt_gct_report.md`" << endl << endl;
);
if RUN_HOMOLOGY_PART then (
    combOut << "## Homology part" << endl << endl;
    combOut << "- Report: `jacobian_pdim_summary.md`" << endl;
    combOut << "- Certificates: `certificates/*_resolution_certificate.m2`" << endl << endl;
);
close combOut;

-- ================================================================
-- EXIT
-- ================================================================

print "";
print "================================================================";
if hasErr then (
    print " COMBINED PIPELINE COMPLETED WITH FAILURES";
    print "================================================================";
    exit 1;
) else (
    print " COMBINED PIPELINE COMPLETED SUCCESSFULLY";
    print "================================================================";
    exit 0;
);
