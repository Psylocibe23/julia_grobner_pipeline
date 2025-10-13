# verify_quotient_space_dimension.sage
# Usage:  sage scripts/verify_quotient_space_dimension.sage path/to/LEX.txt

from sage.all import *
import sys, re

def parse_header(lines):
    meta = {"variables": None, "field": None, "order": None,
            "shape_position": None, "zero_dim": None}
    for s in lines:
        if not s.startswith("#"):
            break
        if s.startswith("# Variables:"):
            meta["variables"] = [v.strip() for v in s.split(":",1)[1].split(",")]
        elif s.startswith("# Field:"):
            meta["field"] = s.split(":",1)[1].strip()
        elif s.startswith("# Order:"):
            meta["order"] = s.split(":",1)[1].strip()
        elif s.startswith("# Shape position:"):
            meta["shape_position"] = (s.split(":",1)[1].strip().lower()=="true")
        elif s.startswith("# Zero-dimensional:"):
            meta["zero_dim"] = (s.split(":",1)[1].strip().lower()=="true")
    return meta

def build_ring(meta):
    fld = (meta.get("field") or "").upper()
    if fld.startswith("GF("):
        m = re.match(r"GF\((\d+)\)", fld)
        if not m: raise ValueError(f"Unrecognized field spec: {fld}")
        K = GF(int(m.group(1)))
    elif fld in ("QQ","Q","RATIONALS"):
        K = QQ
    else:
        m = re.search(r"(\d+)", fld)
        K = GF(int(m.group(1))) if m else GF(2)

    vars_list = meta.get("variables")
    if not vars_list:
        raise ValueError("Variables not found in header.")
    order = (meta.get("order") or "lex").lower()
    if order not in ("lex","deglex","degrevlex","drl"):
        order = "lex"

    R = PolynomialRing(K, names=vars_list, order=order)
    R.inject_variables()
    return R

def read_basis(filepath, R):
    G = []
    in_basis = False
    with open(filepath, "r") as f:
        for raw in f:
            s = raw.strip()
            if not s: continue
            if s.startswith("# --- Groebner basis"):
                in_basis = True; continue
            if not in_basis or s.startswith("#"):
                continue
            G.append(R(s))
    if not G:
        raise ValueError("No LEX basis found under the '--- Groebner basis' header.")
    return G

def triangular_product_if_any(G):
    exps = {}
    for g in G:
        lm = g.lm()
        degs = lm.degrees()
        nz = [(i,e) for i,e in enumerate(degs) if e!=0]
        if len(nz)!=1:
            return (False, {}, None)
        i,e = nz[0]
        if i in exps and exps[i]!=e:
            return (False, {}, None)
        exps[i] = e
    prod = Integer(1)
    for i in exps:
        prod *= exps[i]
    return (True, exps, prod)

def main():
    if len(sys.argv)!=2:
        print("Usage: sage scripts/verify_quotient_space_dimension.sage path/to/LEX.txt")
        sys.exit(1)
    path = sys.argv[1]
    with open(path,"r") as f:
        lines = f.readlines()

    meta = parse_header(lines)
    R = build_ring(meta)
    G = read_basis(path, R)
    I = R.ideal(G)

    print("========== LEX Basis Dim-Check ==========")
    print(f"File: {path}")
    print(f"Ring: {R}")
    print(f"Field: {meta.get('field')}, Order: {meta.get('order')}")
    if meta.get("shape_position") is not None:
        print(f"Shape position (header): {meta['shape_position']}")
    if meta.get("zero_dim") is not None:
        print(f"Zero-dimensional (header): {meta['zero_dim']}")

    # Sanity: Krull dimension
    try:
        kd = I.dimension()
        print(f"Krull dimension computed: {kd}")
    except Exception as e:
        print(f"Krull dimension: unavailable ({e})")

    # --- Singular path ---
    names = list(R.variable_names())
    names_str = "(" + ",".join(names) + ")"     # <-- IMPORTANT: a single string
    K = R.base_ring()
    char = 0 if not K.is_finite() else K.characteristic()

    singR = singular.ring(char, names_str, 'lp')   # lex == 'lp' in Singular
    singular.setring(singR)

    # Convert polynomials to this ring
    sing_polys = [singular(str(g)) for g in G]
    singI = singular.ideal(sing_polys)
    singStd = singular.std(singI)                  # Gröbner (standard) basis in Singular

    # Vector-space dimension of quotient (for 0-dim ideals)
    dI = int(singular.vdim(singStd))
    print(f"d_I (via Singular std+vdim) = {dI}")

    # Triangular cross-check from leading terms
    tri, exps, prod = triangular_product_if_any(G)
    if tri:
        e_list = [(names[i], exps[i]) for i in sorted(exps)]
        print("Triangular detected from LT(g):", e_list)
        print(f"∏ exponents = {prod}  {'== d_I ✔' if prod==dI else '!= d_I ✖'}")
    else:
        print("LEX basis is not strictly triangular (some LT(g) not a pure power).")

    # Reporting: max total degree among LEX polynomials
    print(f"Max total degree among LEX polynomials: {max(g.total_degree() for g in G)}")

if __name__ == "__main__":
    main()
