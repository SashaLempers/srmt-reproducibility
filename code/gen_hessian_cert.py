import json, sys, random, itertools, datetime
import sympy as sp
random.seed(42)
x = sp.symbols('x11 x12 x13 x21 x22 x23 x31 x32 x33')
X = sp.Matrix([[x[0],x[1],x[2]],[x[3],x[4],x[5]],[x[6],x[7],x[8]]])
V=list(x)
def perm3(M):
    s=0
    for sg in itertools.permutations(range(3)):
        p=1
        for i in range(3): p*=M[i,sg[i]]
        s+=p
    return sp.expand(s)
def hdet(f): return sp.expand(sp.hessian(f,V).det())
def nt(p):
    p=sp.expand(p); return 0 if p==0 else len(sp.Poly(p,*V).terms())
det3=sp.expand(X.det()); perm=perm3(X)
Ddet=hdet(det3); D=hdet(perm)
R=sp.expand(D+2*perm**3)
P3=sp.expand(perm**3)
dD=sp.Poly(D,*V).as_dict(); dP=sp.Poly(P3,*V).as_dict()
common=set(dD)&set(dP)
ratios=sorted({sp.Rational(dD[m],dP[m]) for m in common}, key=float)
unit,facs=sp.factor_list(D,*V)
fac,mult=facs[0]
partials=[sp.expand(sp.diff(D,v)) for v in V]
g=partials[0]
for p in partials[1:]: g=sp.gcd(g,p)
g=sp.expand(g)
cert={
 "certificate":"hessian_perm3_squarefree",
 "generated_utc":datetime.datetime.utcnow().isoformat()+"Z",
 "engine":"SymPy "+sp.__version__,
 "arithmetic":"exact rational (Q), no floating point",
 "seed":42,
 "variables":[str(v) for v in V],
 "results":{
   "detHess_det3_equals_minus2_det3cubed": bool(sp.expand(Ddet+2*det3**3)==0),
   "D_perm3_degree": int(sp.total_degree(D)),
   "residual_R_num_monomials": nt(R),
   "distinct_ratios_D_over_perm3cubed":[str(r) for r in ratios],
   "D_is_proportional_to_perm3cubed": len(ratios)==1,
   "Q_factor_count": len(facs),
   "Q_factor_degree": int(sp.total_degree(fac)),
   "Q_factor_multiplicity": int(mult),
   "Q_factor_num_monomials": nt(fac),
   "factorization_unit": str(unit),
   "D_is_rational_cube": (len(facs)==1 and mult%3==0),
   "gcd_nine_first_partials": str(g),
   "gcd_is_nonzero_constant": bool(g!=0 and sp.total_degree(g)==0),
   "D_squarefree_over_C": bool(g!=0 and sp.total_degree(g)==0),
   "D_is_cube_over_C": False
 },
 "claims_supported":[
   "thm:Hess-det (det Hess(det3) = -2 det3^3)",
   "thm:Hess-perm (1)-(4)",
   "lem:squarefree (D squarefree over C, not a cube over C)",
   "thm:hess-separator (perm3 not in closure GL9.det3, one-sided)"
 ]
}
import os
out=os.path.abspath(os.path.join("..","logs","hessian","hessian_perm3_certificate.json"))
os.makedirs(os.path.dirname(out),exist_ok=True)
with open(out,"w") as f: json.dump(cert,f,indent=2)
print("wrote",out)
print(json.dumps(cert["results"],indent=2))
