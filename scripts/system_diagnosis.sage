###############################################################################
# system_diagnosis.sage
#
# This script reads a polynomial system from file, constructs the corresponding
# polynomial ring and ideal, and prints a variety of diagnostic/statistical
# properties of the system. These include variable count, sparsity, degrees,
# homogeneity, Krull dimension, Hilbert polynomial, and more.
# It is written for use in algebraic cryptanalysis and system profiling.
###############################################################################

import sys

############################
# 1. Input parsing function
############################
def parse_input_system(filename):
    """
    Parse the input file describing the polynomial system.

    Input format:
        - First line: comma-separated variable names
        - Second line: field spec ("p" or "p^n" for GF(p^n))
        - Subsequent lines: one polynomial per line (as string), using the variable names

    Returns:
        - var_names: list of variable names (as strings)
        - F: the corresponding finite field (GF(p^n))
        - field_desc: string description of the field
        - polys: list of polynomial strings (for later parsing)
    """
    with open(filename) as f:
        lines = [l.strip() for l in f if l.strip() and not l.strip().startswith('#')]
    var_names = [v.strip() for v in lines[0].split(',')]
    field_line = lines[1]
    if "^" in field_line:
        p, n = [int(x.strip()) for x in field_line.split("^")]
        F = GF(p**n, 'a')
        field_desc = f"GF({p}^{n})"
    else:
        p = int(field_line)
        F = GF(p)
        field_desc = f"GF({p})"
    polys = lines[2:]
    return var_names, F, field_desc, polys

############################
# 2. Homogeneity check
############################
def is_homogeneous(poly):
    """
    Check if a polynomial is homogeneous.
    """
    degs = [m.degree() for m in poly.monomials()]
    return len(set(degs)) == 1

############################
# 3. Main diagnostic logic
############################
def main():
    if len(sys.argv) < 2:
        print("Usage: sage diagnose_system.sage <input_file.in>")
        sys.exit(1)
    infile = sys.argv[1]
    outlog = infile.replace(".in", "_diagnostics.log").replace("data/", "logs/")

    # Parse the input file
    vars, F, field_desc, polys_str = parse_input_system(infile)
    R = PolynomialRing(F, vars)
    polynomials = [R(s) for s in polys_str]
    I = R.ideal(polynomials)

    # --- Basic properties ---
    log = []
    log.append(f"System: {infile}")
    log.append(f"Field: {field_desc}")
    log.append(f"Variables: {vars}")
    log.append(f"Number of variables: {len(vars)}")
    log.append(f"Number of equations: {len(polynomials)}")

    # --- Degree, sparsity, and homogeneity ---
    degs = [p.total_degree() for p in polynomials]
    log.append(f"Degrees of input polynomials: {degs}")
    log.append(f"Max degree: {max(degs)}, Min degree: {min(degs)}")
    terms = [len(p.monomials()) for p in polynomials]
    log.append(f"Sparsity: Terms per polynomial: {terms}")
    homog = [is_homogeneous(p) for p in polynomials]
    log.append(f"Homogeneous polynomials: {homog}")
    all_homog = all(homog)
    log.append(f"All polynomials homogeneous? {all_homog}")

    # --- Krull dimension: geometric structure of solution set ---
    try:
        dim = I.dimension()
        log.append(f"Krull dimension of the ideal: {dim}")
        if dim == 0:
            log.append("System is zero-dimensional (finite number of solutions).")
        elif dim == -1:
            log.append("Ideal is the unit ideal (1 in the basis).")
        else:
            log.append("System is positive-dimensional (infinitely many solutions).")
    except Exception as e:
        log.append(f"Could not compute dimension: {e}")

    # --- Degree of the ideal (Hilbert polynomial) ---
    # Ideal degree (Hilbert polynomial at zero for zero-dim)
    try:
        hilb = I.hilbert_polynomial()
        # For zero-dimensional ideals, hilb(0) equals the number of solutions (counted with multiplicity)
        degI = hilb(0) if I.dimension() == 0 else hilb.degree()
        log.append(f"Degree of the ideal (zero-dim: number of points, else: degree of Hilbert polynomial): {degI}")
    except Exception as e:
        log.append(f"Could not compute ideal degree (Hilbert polynomial): {e}")

    # No shape position/lex basis attempted!
    log.append("Shape position check skipped (lex Groebner basis not computed).")

    # --- Special forms: Boolean and quadratic ---
    is_boolean = F.order() == 2 and all(p.total_degree() <= 2 for p in polynomials)
    log.append(f"System is Boolean quadratic? {'Yes' if is_boolean else 'No'}")
    is_quadratic = all(p.total_degree() <= 2 for p in polynomials)
    log.append(f"System is quadratic (all degree ≤2)? {'Yes' if is_quadratic else 'No'}")
    
    # Save log
    with open(outlog, "w") as f:
        for line in log:
            print(line)
            f.write(line + "\n")
    print(f"\nDiagnostics saved to {outlog}")

main()
