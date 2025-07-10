import sys
import time
import os

def read_groebner_basis_file(result_file):
    with open(result_file, 'r') as f:
        lines = f.readlines()
        variables = None
        p = None
        polys = []
        for i, line in enumerate(lines):
            if line.startswith("# Variables:"):
                variables = [v.strip() for v in line.split(":",1)[1].split(",")]
            elif line.startswith("# Field characteristic:"):
                p = int(line.split(":")[1].strip())
            elif line.startswith("# Field: GF("):
                p = int(line.split("GF(")[1].split(")")[0].strip())
            elif line.startswith("# --- Groebner basis ---"):
                basis_start = i + 1
                break
        else:
            raise ValueError("Could not find basis start in file.")
    # The rest are basis polynomials
    for line in lines[basis_start:]:
        s = line.strip()
        if s:
            polys.append(s)
    return variables, p, polys

def is_shape_position(G_lex):
    # Simple check: all polys are univariate (shape position for typical cryptanalytic systems)
    try:
        return all(len(g.lm().variables()) == 1 for g in G_lex)
    except Exception:
        return False

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main():
    if len(sys.argv) != 2:
        print("Usage: sage scripts/convert_to_lex_fglm.sage path/to/result_file.txt")
        sys.exit(1)
    result_file = sys.argv[1]
    variables, p, polys = read_groebner_basis_file(result_file)
    print(f"Inferred variables: {variables}")
    print(f"Field characteristic: {p}")

    # Prepare output filenames in results/logs
    base_name = os.path.splitext(os.path.basename(result_file))[0]
    results_dir = "results"
    logs_dir = "logs"
    ensure_dir(results_dir)
    ensure_dir(logs_dir)
    lex_outfile = os.path.join(results_dir, base_name + "_LEX.txt")
    log_outfile = os.path.join(logs_dir, base_name + "_FGLM.log")

    # Input DRL basis stats
    input_basis_size = len(polys)
    input_degrees = []
    try:
        R_temp = PolynomialRing(GF(p), variables)
        input_degrees = [R_temp(s).total_degree() for s in polys]
        input_max_deg = max(input_degrees)
    except Exception:
        input_degrees = []
        input_max_deg = None

    # Construct polynomial ring in LEX order
    R_lex = PolynomialRing(GF(p), variables, order='lex')
    G_lex = [R_lex(s) for s in polys]
    I_lex = Ideal(G_lex)

    # FGLM and stats
    t0 = time.time()
    G_lex_fglm = I_lex.groebner_basis(algorithm="singular:stdfglm")
    t1 = time.time()
    fglm_time = t1 - t0

    # Output LEX basis stats
    output_basis_size = len(G_lex_fglm)
    output_degrees = [g.total_degree() for g in G_lex_fglm]
    output_max_deg = max(output_degrees)
    shape_pos = is_shape_position(G_lex_fglm)

    # Attempt to compute number of solutions (for zero-dim ideals)
    num_solutions = None
    zero_dim = None
    try:
        zero_dim = (I_lex.dimension() == 0)
        if zero_dim:
            sols = I_lex.variety()
            num_solutions = len(sols)
        else:
            num_solutions = None
    except Exception:
        zero_dim = None
        num_solutions = None

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

