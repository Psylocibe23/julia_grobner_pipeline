#!/usr/bin/env sage
###############################################################################
# expand_field_extension_to_base.sage
#
# Given a polynomial system over GF(p^n), rewrite it as an equivalent system
# over GF(p), introducing n variables for each original variable.
# This is required for Gröbner basis computations over the base field and is
# standard in algebraic cryptanalysis when only base field arithmetic is available.
###############################################################################

import sys
import os

###############################################################################
# 1. Parse input file: variable names, field, equations
###############################################################################
def parse_input_file(filename):
    """
    Parse .in file with:
      - Line 1: var1, var2, ..., varn
      - Line 2: p or p^n
      - Next lines: polynomial equations
    Returns: (var_names, p, n, equations)
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

###############################################################################
# 2. Expansion logic: replace each GF(p^n) variable with n base field variables
###############################################################################
def expand_system(var_names, p, n, polys):
    """
    Given polynomials in variables x1,...,xm over GF(p^n), expand the system
    to a set of equations over GF(p), introducing n variables for each xi,
    representing the coordinates of xi in a fixed basis of GF(p^n) over GF(p).
    
    Returns:
      - new_var_names: ["x1_0", ..., "x1_{n-1}", ..., "xm_{n-1}"]
      - expanded_polys: list of equations over GF(p) in the new variables
    """
    if n == 1:
        # Trivial case: system is already over GF(p)
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

    # Substitution: x_i → x_i_0 + a*x_i_1 + ... + a^{n-1}*x_i_{n-1}
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

    # For each input polynomial, expand over the base field
    expanded_polys = []
    for poly_s in polys:
        poly_K = PR_K(poly_s)
        # Substitute each x → linear combination of base field vars
        poly_sub = poly_K.subs(subst)
        # Decompose into n coordinate polynomials over Fp
        coord_polys = [PR_p(0) for _ in range(n)]
        for mono, coeff in poly_sub.dict().items():
            # coeff is in GF(p^n): expand as n-tuple over GF(p)
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

###############################################################################
# 3. Output: Write expanded system in .in format for Julia pipeline, etc.
###############################################################################
def write_output_file(var_names, p, n, polys, output_path):
    """
    Write expanded system in .in format:
      - First line: comma-separated variable names
      - Second line: base field (p)
      - Following lines: equations
    """
    with open(output_path, "w") as f:
        f.write(", ".join(var_names) + "\n")
        f.write(str(p) + "\n")
        for poly in polys:
            f.write(poly + "\n")

###############################################################################
# 4. Main routine
###############################################################################
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
