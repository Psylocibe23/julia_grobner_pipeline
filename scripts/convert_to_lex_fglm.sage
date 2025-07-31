###############################################################################
# convert_to_lex_fglm.sage
#
# Given a DRL Gröbner basis (e.g., output of F4 or F5), this script converts
# the basis to LEX order using the FGLM algorithm, checks for "shape position,"
# and logs key algebraic and computational statistics.
# Requirements: SAGE, DRL basis output as produced by pipeline
###############################################################################

import sys
import time
import os

############################
# 1. Parse the DRL basis output file
############################
def read_groebner_basis_file(result_file):
    """
    Parse an output file containing a DRL Gröbner basis.
    Returns:
        variables: list of variable names (str)
        p: characteristic of the base field
        polys: list of basis polynomials as strings
    """
    with open(result_file, 'r') as f:
        lines = f.readlines()
        variables = None
        p = None
        polys = []
        basis_start = None
        for i, line in enumerate(lines):
            if line.startswith("# Variables:"):
                variables = [v.strip() for v in line.split(":",1)[1].split(",")]
            elif line.startswith("# Field characteristic:"):
                p = int(line.split(":")[1].strip())
            elif line.startswith("# Field: GF("):
                # Support both "# Field characteristic:" and "# Field: GF(...)" lines
                p = int(line.split("GF(")[1].split(")")[0].strip())
            elif line.startswith("# --- Groebner basis ---"):
                basis_start = i + 1
                break
        else:
            raise ValueError("Could not find basis start in file.")
    # The rest are basis polynomials as strings
    for line in lines[basis_start:]:
        s = line.strip()
        if s and not s.startswith("#"):
            polys.append(s)
    return variables, p, polys

############################
# 2. Shape position check
############################
def is_shape_position(G_lex):
    """
    A system is in shape position if its LEX Gröbner basis is in 'triangular' form:
    each basis element involves only one new variable (i.e., univariate), which allows
    a "back-substitution" solution (like in triangular linear systems).
    This is a sufficient condition for efficiently extracting all solutions by root finding.
    """
    # Simple check: all polys are univariate (shape position for typical cryptanalytic systems)
    try:
        return all(len(g.lm().variables()) == 1 for g in G_lex)
    except Exception:
        return False

############################
# 3. Utility: ensure output directories exist
############################
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

############################
# 4. Main conversion and logging logic
############################
def main():
    if len(sys.argv) != 2:
        print("Usage: sage scripts/convert_to_lex_fglm.sage path/to/result_file.txt")
        sys.exit(1)
    result_file = sys.argv[1]

    # === Parse input basis file (with robust diagnostics) ===
    variables, p, polys = read_groebner_basis_file(result_file)
    print(f"Inferred variables: {variables}")
    print(f"Field characteristic: {p}")
    print(f"Number of input polynomials: {len(polys)}")

    # === Prepare output filenames in results/logs ===
    base_name = os.path.splitext(os.path.basename(result_file))[0]
    results_dir = "results"
    logs_dir = "logs"
    ensure_dir(results_dir)
    ensure_dir(logs_dir)
    lex_outfile = os.path.join(results_dir, base_name + "_LEX.txt")
    log_outfile = os.path.join(logs_dir, base_name + "_FGLM.log")

    # === Input DRL basis stats ===
    input_basis_size = len(polys)
    input_degrees = []
    try:
        # DRL ring for statistics only
        R_temp = PolynomialRing(GF(p), variables, order='deglex')
        input_degrees = [R_temp(s).total_degree() for s in polys]
        input_max_deg = max(input_degrees)
    except Exception as e:
        input_degrees = []
        input_max_deg = None
        print(f"Warning: could not compute DRL degree stats: {e}")

    # === Construct polynomial rings ===
    # We *always* reconstruct rings with the correct variables and order
    R_drl = PolynomialRing(GF(p), variables, order='deglex')
    G_drl = [R_drl(s) for s in polys]
    I_drl = R_drl.ideal(G_drl)

    # Now switch to LEX for FGLM
    R_lex = PolynomialRing(GF(p), variables, order='lex')
    # Map polynomials to the LEX ring (parse as strings)
    G_lex = [R_lex(str(p)) for p in G_drl]
    I_lex = R_lex.ideal(G_lex)

    # === FGLM computation ===
    t0 = time.time()
    G_lex_fglm = I_lex.groebner_basis(algorithm="singular:stdfglm")
    t1 = time.time()
    fglm_time = t1 - t0

    # === Output LEX basis stats ===
    output_basis_size = len(G_lex_fglm)
    output_degrees = [g.total_degree() for g in G_lex_fglm]
    output_max_deg = max(output_degrees)
    shape_pos = is_shape_position(G_lex_fglm)


    print("\nLEX Groebner basis via FGLM:")
    for g in G_lex_fglm:
        print(g)
    with open(lex_outfile, "w") as out:
        out.write(f"# Lex Groebner basis for {result_file}\n")
        out.write(f"# Variables: {', '.join(variables)}\n")
        out.write(f"# Field: GF({p})\n")
        for g in G_lex_fglm:
            out.write(str(g) + "\n")
    print(f"\nSaved to {lex_outfile}")

    # === Save all diagnostic and complexity information ===
    with open(log_outfile, "w") as log:
        log.write(f"# FGLM conversion log for {result_file}\n")
        log.write(f"# Input DRL Groebner basis: {input_basis_size} polys, max deg = {input_max_deg}, degrees = {input_degrees}\n")
        log.write(f"# Output LEX Groebner basis: {output_basis_size} polys, max deg = {output_max_deg}, degrees = {output_degrees}\n")
        log.write(f"# FGLM wall time: {fglm_time:.5f} seconds\n")
        log.write(f"# LEX basis shape position: {shape_pos}\n")
        log.write(f"# System zero-dimensional: {zero_dim}\n")
        log.write(f"# Number of solutions: {num_solutions}\n")
        log.write(f"# Input: {result_file}\n")
        log.write(f"# Output: {lex_outfile}\n")
    print(f"FGLM stats log saved to {log_outfile}")

if __name__ == "__main__":
    main()
