###############################################################################
# fast_extract_solutions_from_shape_lex.sage
#
# Efficiently extract all solutions from a LEX Gröbner basis in shape position.
# Requirements: SAGE, output from convert_to_lex_fglm.sage
###############################################################################

import sys, os, time

def parse_lex_basis(filename):
    # Returns variables (list), field characteristic (int), and polynomials (list of strings)
    variables, p, polys = None, None, []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("# Variables:"):
            variables = [v.strip() for v in line.split(":",1)[1].split(",")]
        elif line.startswith("# Field: GF("):
            p = int(line.split("GF(")[1].split(")")[0].strip())
        elif line and not line.startswith("#"):
            polys.append(line.strip())
    return variables, p, polys

def main():
    if len(sys.argv) < 2:
        print("Usage: sage fast_extract_solutions_from_shape_lex.sage <LEX_basis_file>")
        sys.exit(1)
    basisfile = sys.argv[1]
    variables, p, poly_strs = parse_lex_basis(basisfile)
    F = GF(p)
    R = PolynomialRing(F, variables)
    polys = [R(s) for s in poly_strs if s]
    n = len(variables)

    # Ensure shape position: each poly involves <=1 new variable than previous
    shapes = [set(f.variables()) for f in polys]
    is_shape = all(len(v) == 1 for v in shapes[:-1]) and len(shapes[-1]) == 1
    if not is_shape:
        print("System is not in shape position. Aborting fast extraction.")
        sys.exit(2)
    t0 = time.time()
    # Last polynomial is univariate: g(x_n)
    univar_poly = polys[-1]
    univar_var = list(univar_poly.variables())[0]
    roots = univar_poly.roots(multiplicities=False)
    solutions = []
    for alpha in roots:
        assign = {univar_var: alpha}
        for poly in reversed(polys[:-1]):
            var = list(poly.variables())[0]
            rhs = -(poly.subs(assign) - var)
            assign[var] = rhs
        # Return in order of variables
        solution = {str(v): assign[R(v)] for v in variables}
        solutions.append(solution)
    t1 = time.time()
    print(f"Found {len(solutions)} solutions in {t1-t0:.3f} seconds.")
    # Optionally write solutions to file
    out_file = os.path.splitext(basisfile)[0] + "_fastsols.txt"
    with open(out_file, "w") as out:
        for sol in solutions:
            out.write("{" + ", ".join(f"{k}: {v}" for k, v in sol.items()) + "}\n")
    print(f"Solutions saved to {out_file}")

if __name__ == "__main__":
    main()
