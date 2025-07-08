import sys

def parse_field_vars_basis(filename):
    """
    Parse variables and field characteristic from the header of the Groebner basis result file.
    Accepts both "# Field characteristic: p" and "# Field: GF(p)".
    Returns (vars, p, lines) where vars is a list of variable names, p is the characteristic, and lines are the basis polys.
    """
    with open(filename) as f:
        lines = [l.strip() for l in f if l.strip() and not l.strip().startswith('#')]
    # Try to parse variables/field from comments
    vars_line, char_line, gf_line = None, None, None
    with open(filename) as f:
        for line in f:
            if line.startswith("# Variables:"):
                vars_line = line.strip().split(":",1)[1].strip()
            if line.startswith("# Field characteristic:"):
                char_line = line.strip().split(":",1)[1].strip()
            if line.startswith("# Field: GF("):
                gf_line = line.strip().split("GF(",1)[1].split(")")[0].strip()
    if vars_line:
        vars = [v.strip() for v in vars_line.split(",")]
    else:
        # Fallback: try to parse from polynomials
        poly_vars = set()
        for poly in lines:
            for v in poly.replace("^"," ").replace("*"," ").replace("+"," ").replace("-"," ").replace("("," ").replace(")"," ").split():
                if v.isalpha(): poly_vars.add(v)
        vars = sorted(poly_vars)
    # Prefer explicit field header, fallback to 2 if missing
    if char_line:
        p = char_line
    elif gf_line:
        p = gf_line
    else:
        p = "2" # fallback, but should warn!
        print("Warning: Field characteristic not found, defaulting to 2.")
    return vars, p, lines

def main():
    if len(sys.argv) < 2:
        print("Usage: sage extract_solutions_from_lex.sage <basis_file>")
        sys.exit(1)
    basisfile = sys.argv[1]
    vars, p, polys = parse_field_vars_basis(basisfile)
    print(f"Inferred variables: {vars}")
    print(f"Field characteristic: {p}")
    # Extension field detection (handles GF(p^k) if needed)
    if "^" in str(p):
        q, k = p.split("^")
        q = int(q.strip())
        k = int(k.strip())
        F = GF(q**k, 'a')  # primitive element 'a'
    else:
        F = GF(int(p))
    R = PolynomialRing(F, vars)
    polynomials = [R(poly) for poly in polys]
    I = R.ideal(polynomials)
    print(f"Solving system in {R}...")

    try:
        sols = I.variety()
    except Exception as e:
        print("variety() failed, trying brute force...")
        import itertools
        sols = []
        for vals in itertools.product(F, repeat=len(vars)):
            subst = dict(zip(vars, vals))
            if all(poly(**subst) == 0 for poly in polynomials):
                sols.append(subst)

    # Double-check brute force for all fields (robust for small/moderate fields)
    if not sols:
        import itertools
        sols = []
        for vals in itertools.product(F, repeat=len(vars)):
            subst = dict(zip(vars, vals))
            if all(poly(**subst) == 0 for poly in polynomials):
                sols.append(subst)

    result_file = basisfile.replace(".txt", "_sols.txt")
    with open(result_file, "w") as out:
        out.write(f"# Solutions for {basisfile}\n")
        for sol in sols:
            # Display variable assignments in order
            out.write("{" + ", ".join(f"{v}: {sol[v]}" for v in vars) + "}\n")
    print(f"Solutions saved to {result_file}")
    print(f"Number of solutions: {len(sols)}")

    # Verify all solutions
    for sol in sols:
        strsol = {str(k): v for k, v in sol.items()}
        for poly in polynomials:
            assert poly(**strsol) == 0
    print("Verification completed.")

main()
