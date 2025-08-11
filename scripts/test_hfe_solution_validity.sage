###############################################################################
# test_hfe_solution_validity.sage
#
# PURPOSE
#   Given:
#     • a pipeline .in file (variables, field, public equations, field equations),
#     • a solutions file (one dict-like solution per line),
#     • optionally, a generation log (to recover the secret and S),
#   this script verifies:
#     (1) whether each candidate satisfies all *public* equations,
#     (2) whether each candidate satisfies the *field equations* x_i^2 + x_i,
#     (3) if a secret is available from the log, whether the candidate equals it,
#     (4) if S is available from the log, whether the candidate is S-equivalent
#         to the secret:  A_S * x + b_S == A_S * x* + b_S over GF(2).
#
# USAGE
#   sage scripts/test_hfe_solution_validity.sage <in_file> <solutions_file> [log_file]
#
# SOLUTIONS FORMAT
#   Expected one solution per line in a dict-like form, e.g.:
#       {x0: 0, x1: 1, x2: 0, x3: 1}
#   Keys may or may not be quoted; values must be 0/1 (GF(2)).
#
###############################################################################

import sys, re, ast
from sage.all import GF, PolynomialRing, Matrix, vector

# -----------------------------------------------------------------------------
# 1) Parse pipeline .in file
# -----------------------------------------------------------------------------
def parse_in_file(infile):
    """
    Returns:
        var_names: list[str]
        F: GF(p) (we expect p=2)
        R: polynomial ring GF(p)[x...]
        public_polys: list of R-polynomials (length n)
        field_polys: list of R-polynomials (length n) or None
    """
    with open(infile, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if len(lines) < 2:
        raise ValueError("Malformed .in file (need at least 2 header lines).")

    # Line 1: variable names "x0, x1, ..., x{n-1}"
    var_names = [v.strip() for v in lines[0].split(",")]
    n = len(var_names)
    if n == 0:
        raise ValueError("No variables found in first line of .in file.")

    # Line 2: field (accepts "2", "GF(2)", "GF(2^1)")
    fld = lines[1]
    if fld == "2":
        p = 2
    elif fld.upper().startswith("GF("):
        inside = fld.split("GF(", 1)[1].split(")", 1)[0].strip()
        if "^" in inside:
            base, _deg = inside.split("^", 1)
            p = int(base.strip())
        else:
            p = int(inside)
    else:
        p = int(fld)

    F = GF(p)
    R = PolynomialRing(F, var_names, order='lex')
    x = R.gens()

    # Remaining lines: first n are public equations; next n (if present) are field eqs
    eq_lines = lines[2:]
    if len(eq_lines) < n:
        raise ValueError(f"Expected at least {n} public equations; found {len(eq_lines)}.")

    public_polys = [R(s) for s in eq_lines[:n]]

    field_polys = None
    if len(eq_lines) >= 2*n:
        cand = [R(s) for s in eq_lines[n:2*n]]
        ok = all(cand[i] == x[i]**2 + x[i] for i in range(n))
        field_polys = cand if ok else None

    return var_names, F, R, public_polys, field_polys

# -----------------------------------------------------------------------------
# 2) Parse solutions file
# -----------------------------------------------------------------------------
def parse_solutions_file(solfile):
    """
    Parse lines like '{x0: 0, x1: 1, ...}' into dicts.
    Accepts optional quotes; ignores other non-solution lines.
    """
    sols = []
    with open(solfile) as f:
        for line in f:
            line = line.strip()
            if not line or "{" not in line or "}" not in line:
                continue
            # ensure JSON-like keys by quoting bare words before ':'
            fixed = re.sub(r'(\w+)\s*:', r'"\1":', line)
            try:
                d = ast.literal_eval(fixed)
            except Exception:
                continue
            # values to ints
            sol = {}
            for k, v in d.items():
                sol[str(k)] = int(v)
            sols.append(sol)
    return sols

# -----------------------------------------------------------------------------
# 3) Evaluate equations at a candidate solution (ordered by var_names)
# -----------------------------------------------------------------------------
def eval_polys(polys, R, var_names, sol_dict):
    """
    Evaluate each polynomial at the point specified by sol_dict,
    interpreting values in the ring's base field. Returns list[bool].
    """
    F = R.base_ring()
    # Build value tuple in the ring's variable order
    try:
        vals = tuple(F(int(sol_dict[name])) for name in var_names)
    except KeyError as e:
        missing = str(e).strip("'")
        raise KeyError(f"Solution is missing variable '{missing}'.")
    return [p(*vals) == 0 for p in polys]

# -----------------------------------------------------------------------------
# 4) Parse generator log to recover secret and (if present) A_S, b_S
# -----------------------------------------------------------------------------
def load_secret_and_S_from_log(logfile):
    """
    Try to recover:
        secret_vector : list[int] or None
        A_S : Matrix(GF(2)) or None
        b_S : vector(GF(2)) or None
    Supports several formats:
      • "Secret: [1,0,1,...]"  or  "Secret = [ ... ]"
      • "Affine input map S(x) = A_S x + b_S:"  followed by matrix + vector
      • "A_S =" block (matrix) followed by b_S line
    Returns (secret_vector_or_None, A_S_or_None, b_S_or_None).
    """
    txt = open(logfile, "r").read()

    # Secret (both ':' and '=' variants, optionally with 'x*')
    m = re.search(r"Secret(?:\s*x\*)?\s*[:=]\s*\[([01,\s]+)\]", txt, flags=re.IGNORECASE)
    secret = None
    if m:
        secret = [int(t) for t in re.split(r"[,\s]+", m.group(1).strip()) if t != ""]

    # Try to parse A_S, b_S in several formats
    A_S = None
    b_S = None

    # Helper to parse a Sage-style matrix block like:
    # [1 0 1]
    # [0 1 1]
    def parse_matrix_block(lines, start):
        rows = []
        i = start
        while i < len(lines) and lines[i].strip().startswith("["):
            row = [int(s) for s in lines[i].strip().replace("[","").replace("]","").split()]
            rows.append(row); i += 1
        return rows, i

    lines = txt.splitlines()

    # Pattern 1: explicitly labeled "A_S =" then rows, then a b_S line
    for i, ln in enumerate(lines):
        if re.match(r"^\s*A_S\s*=", ln):
            rows, j = parse_matrix_block(lines, i+1)
            if rows:
                A_S = Matrix(GF(2), rows)
                # next non-empty line try to parse as b_S
                while j < len(lines) and not lines[j].strip():
                    j += 1
                if j < len(lines):
                    b_line = lines[j].strip()
                    # allow "(1,0,...)" or "[1,0,...]" or space-separated
                    b_line = b_line.replace("(", "").replace(")", "").replace("[", "").replace("]", "")
                    b_vals = [int(t) for t in re.split(r"[,\s]+", b_line) if t != ""]
                    b_S = vector(GF(2), b_vals)
            break

    # Pattern 2: older logs "Affine input map S(x) = A_S x + b_S:" followed by matrix + vector
    if A_S is None:
        for i, ln in enumerate(lines):
            if "Affine input map S(x)" in ln and "A_S" in ln and "b_S" in ln:
                rows, j = parse_matrix_block(lines, i+1)
                if rows:
                    A_S = Matrix(GF(2), rows)
                    # parse vector line
                    while j < len(lines) and not lines[j].strip():
                        j += 1
                    if j < len(lines):
                        b_line = lines[j].strip()
                        b_line = b_line.replace("(", "").replace(")", "").replace("[", "").replace("]", "")
                        b_vals = [int(t) for t in re.split(r"[,\s]+", b_line) if t != ""]
                        b_S = vector(GF(2), b_vals)
                break

    return secret, A_S, b_S

# -----------------------------------------------------------------------------
# 5) S-equivalence test
# -----------------------------------------------------------------------------
def S_equivalent(var_names, sol_dict, secret_vector, A_S, b_S):
    """
    Return True iff A_S * sol + b_S == A_S * secret + b_S over GF(2).
    """
    F2 = GF(2)
    x_sol = vector(F2, [F2(int(sol_dict[v])) for v in var_names])
    x_sec = vector(F2, [F2(int(b))          for b in secret_vector])
    return (A_S * x_sol + b_S) == (A_S * x_sec + b_S)

# -----------------------------------------------------------------------------
# 6) Main
# -----------------------------------------------------------------------------
def main():
    if len(sys.argv) < 3:
        print("Usage: sage scripts/test_hfe_solution_validity.sage <in_file> <solutions_file> [log_file]")
        sys.exit(1)

    in_file  = sys.argv[1]
    sols_file = sys.argv[2]
    log_file = sys.argv[3] if len(sys.argv) >= 4 else None

    var_names, F, R, public_polys, field_polys = parse_in_file(in_file)
    n = len(var_names)
    solutions = parse_solutions_file(sols_file)

    # Optional: load secret and S from log (if provided)
    secret, A_S, b_S = (None, None, None)
    if log_file:
        try:
            secret, A_S, b_S = load_secret_and_S_from_log(log_file)
        except Exception as e:
            print(f"[WARN] Could not parse log '{log_file}': {e}")

    print("==============================================================")
    print(" HFE SOLUTION VALIDITY CHECK")
    print("==============================================================")
    print(f"System: {in_file}")
    print(f"Solutions file: {sols_file}")
    if log_file:
        print(f"Log file: {log_file}")
    print(f"Variables: {var_names}  (n={n})")
    print("--------------------------------------------------------------")
    print(f"Found {len(solutions)} solution(s) to test.\n")

    for i, sol in enumerate(solutions, 1):
        print(f"Solution #{i}: {sol}")

        # (1) Public equations
        ok_pub = eval_polys(public_polys, R, var_names, sol)
        print(f"  Satisfies public equations? {'YES' if all(ok_pub) else 'NO'}")
        if not all(ok_pub):
            bad = [k+1 for k,b in enumerate(ok_pub) if not b]
            print(f"    Failing eq indices: {bad}")

        # (2) Field equations (if present)
        if field_polys is not None:
            ok_field = eval_polys(field_polys, R, var_names, sol)
            print(f"  Satisfies field equations x_i^2+x_i? {'YES' if all(ok_field) else 'NO'}")
            if not all(ok_field):
                badf = [k+1 for k,b in enumerate(ok_field) if not b]
                print(f"    Failing field eq indices: {badf}")
        else:
            print("  Field equations not present (or not recognized) in .in — skipped.")

        # (3) Secret equality (if available)
        if secret is not None:
            vals = [int(sol[v]) for v in var_names]
            is_secret = (vals == [int(b) for b in secret])
            print(f"  Equals planted secret? {'YES' if is_secret else 'NO'}")
        else:
            print("  Secret not found in log — equality test skipped.")

        # (4) S-equivalence (if S and secret available)
        if secret is not None and A_S is not None and b_S is not None:
            try:
                eqS = S_equivalent(var_names, sol, secret, A_S, b_S)
                print(f"  S-equivalent to secret (A_S x + b_S)? {'YES' if eqS else 'NO'}")
            except Exception as e:
                print(f"  S-equivalence check failed: {e}")
        else:
            print("  S-equivalence not checked (missing S or secret).")

        print()

    print("--------------------------------------------------------------")
    print("Validation complete.")

if __name__ == "__main__":
    main()
