#!/usr/bin/env sage
import sys
import os

def parse_input_file(filename):
    """Parse input .in file with:
       line 1: var1, var2, ..., varn
       line 2: p or p^n
       next: equations
    """
    with open(filename) as f:
        lines = [l.strip() for l in f if l.strip()]
    var_names = [v.strip() for v in lines[0].split(",")]
    field_line = lines[1]
    equations = lines[2:]
    # Parse field
    if "^" in field_line:
        base, ext = field_line.split("^")
        p = int(base.strip())
        n = int(ext.strip())
    else:
        p = int(field_line)
        n = 1
    return var_names, p, n, equations

def expand_system(var_names, p, n, polys):
    """Expand each poly over F_{p^n} as n polynomials over F_p in new variables."""
    if n == 1:
        return var_names, polys

    Fp = GF(p)
    Fpn = GF(p**n, name="a")
    a = Fpn.gen()
    # Original ring over Fpn
    PR_K = PolynomialRing(Fpn, var_names)
    vars_K = PR_K.gens_dict()
    # New variables for base field system
    new_var_names = []
    for v in var_names:
        new_var_names.extend([f"{v}_{j}" for j in range(n)])
    PR_p = PolynomialRing(Fp, new_var_names)

    # Substitution: each x maps to x_0 + a*x_1 + ... + a^{n-1}*x_{n-1}
    subst = {}
    for idx, v in enumerate(var_names):
        expr = sum(PR_p(f"{v}_{j}") * a**j for j in range(n))
        subst[vars_K[v]] = expr

    # Vector space extraction for coordinate expansion
    VS = Fpn.vector_space()
    if isinstance(VS, tuple):
        V = VS[1]
    else:
        V = VS

    # Expand each equation and extract all basis coordinates
    expanded_polys = []
    for poly_s in polys:
        poly_K = PR_K(poly_s)
        poly_sub = poly_K.subs(subst)
        # We'll use .dict() for full monomial expansion.
        coord_polys = [PR_p(0) for _ in range(n)]
        for mono, coeff in poly_sub.dict().items():
            # coeff is in Fpn, expand as vector in Fp^n
            cvec = V(coeff)
            for i, c in enumerate(cvec):
                if c != 0:
                    # Build monomial in PR_p
                    mono_terms = []
                    for j, exp in enumerate(mono):
                        if exp != 0:
                            mono_terms.extend([new_var_names[j]]*exp)
                    if mono_terms:
                        monostr = "*".join(mono_terms)
                        exprstr = f"{c}*{monostr}" if c != 1 else monostr
                    else:
                        exprstr = f"{c}"
                    # Add term to i-th coordinate poly
                    coord_polys[i] += PR_p(exprstr)
        # Add all n coordinate polys as equations
        for cp in coord_polys:
            if cp != 0:
                expanded_polys.append(str(cp))

    return new_var_names, expanded_polys



def write_output_file(var_names, p, n, polys, output_path):
    """Write expanded system in .in format for Julia pipeline."""
    with open(output_path, "w") as f:
        f.write(", ".join(var_names) + "\n")
        f.write(str(p) + "\n")
        for poly in polys:
            f.write(poly + "\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: sage scripts/expand_field_extension_to_base.sage inputfile.in")
        sys.exit(1)
    infile = sys.argv[1]
    var_names, p, n, polys = parse_input_file(infile)
    print(f"Parsed: Variables: {var_names}, Field: GF({p}^{n})" if n > 1 else f"Field: GF({p})")
    print(f"Original system has {len(var_names)} variables, {len(polys)} polynomials.")

    new_var_names, new_polys = expand_system(var_names, p, n, polys)
    print(f"Expanded system has {len(new_var_names)} variables, {len(new_polys)} equations.")

    # Output file: data/originalfilename_expanded.in
    out_name = os.path.splitext(os.path.basename(infile))[0] + "_expanded.in"
    out_path = os.path.join("data", out_name)
    write_output_file(new_var_names, p, 1, new_polys, out_path)
    print(f"Expanded system written to {out_path}")

if __name__ == "__main__" or "sage" in __name__:
    main()
