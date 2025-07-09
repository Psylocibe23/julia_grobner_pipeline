import sys
import time
import os

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

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main():
    if len(sys.argv) < 2:
        print("Usage: sage extract_solutions_from_lex.sage <basis_file>")
        sys.exit(1)
    basisfile = sys.argv[1]
    vars, p, polys = parse_field_vars_basis(basisfile)
    print(f"Inferred variables: {vars}")
    print(f"Field characteristic: {p}")

    # Set up results and logs paths
    base_name = os.path.splitext(os.path.basename(basisfile))[0]
    results_dir = "results"
    logs_dir = "logs"
    ensure_dir(results_dir)
    ensure_dir(logs_dir)
    result_file = os.path.join(results_dir, base_name + "_sols.txt")
    log_file = os.path.join(logs_dir, base_name + "_SOLUTIONS.log")

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

    t0 = time.time()
    try:
        sols = I.variety()
        method = "variety"
    except Exception as e:
        print("variety() failed, trying brute force...")
        import itertools
        sols = []
        for vals in itertools.product(F, repeat=len(vars)):
            subst = dict(zip(vars, vals))
            if all(poly(**subst) == 0 for poly in polynomials):
                sols.append(subst)
        method = "brute_force"
    t1 = time.time()
    solve_time = t1 - t0

    # Double-check brute force for all fields (robust for small/moderate fields)
    if not sols:
        import itertools
        sols = []
        for vals in itertools.product(F, repeat=len(vars)):
            subst = dict(zip(vars, vals))
            if all(poly(**subst) == 0 for poly in polynomials):
                sols.append(subst)
        method = "brute_force"

    # Write solutions file
    with open(result_file, "w") as out:
        out.write(f"# Solutions for {basisfile}\n")
        for sol in sols:
            # Display variable assignments in order
            out.write("{" + ", ".join(f"{v}: {sol[v]}" for v in vars) + "}\n")
    print(f"Solutions saved to {result_file}")
    print(f"Number of solutions: {len(sols)}")

    # Verify all solutions
    all_verified = True
    for sol in sols:
        strsol = {str(k): v for k, v in sol.items()}
        for poly in polynomials:
            if not poly(**strsol) == 0:
                all_verified = False
                print("Verification failed for solution:", sol)
                break

    print("Verification completed." if all_verified else "Some solutions did not verify.")

    # Save log
    with open(log_file, "w") as log:
        log.write(f"# Solution extraction log for {basisfile}\n")
        log.write(f"# Input LEX basis: {basisfile}\n")
        log.write(f"# Output solutions file: {result_file}\n")
        log.write(f"# Variables: {vars}\n")
        log.write(f"# Field characteristic: {p}\n")
        log.write(f"# Number of equations: {len(polys)}\n")
        log.write(f"# Solution method: {method}\n")
        log.write(f"# Time for solution extraction (sec): {solve_time:.6f}\n")
        log.write(f"# Number of solutions: {len(sols)}\n")
        log.write(f"# Verification passed: {all_verified}\n")
        # Optionally, print solutions if small
        if len(sols) <= 20:
            log.write("# Solutions:\n")
            for sol in sols:
                log.write("{" + ", ".join(f"{v}: {sol[v]}" for v in vars) + "}\n")
    print(f"Extraction log saved to {log_file}")

main()
