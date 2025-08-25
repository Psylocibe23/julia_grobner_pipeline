###############################################################################
# system_diagnosis.sage
#
# This script reads a polynomial system from file, constructs the corresponding
# polynomial ring and ideal, and prints a variety of diagnostic/statistical
# properties of the system. These include variable count, sparsity, degrees,
# homogeneity, Krull dimension, Hilbert polynomial, and more.
# It is written for use in algebraic cryptanalysis and system profiling.
#
# NEW in this version:
#   - Robust time budgets for expensive steps using a separate process:
#       * --dim=<time>     : time budget for I.dimension()
#       * --hilbert=<time> : time budget for I.hilbert_polynomial()
#     <time> may be seconds ("90"), or "60s", "2m", "1m30s"; a trailing '+' is
#     allowed and ignored ("60+" == "60"). 0 or omitted => skip that step.
###############################################################################

import sys, os, re, time, multiprocessing
from sage.all import *

############################
# 0. CLI helpers
############################

def parse_time_budget(tok):
    """
    Parse a time string into seconds.
    Accepts:
      - bare seconds: "90"
      - with units:   "60s", "2m", "1m30s", "1h", "1h2m3s"
      - a trailing '+' is allowed and ignored: "60+" -> 60
      - "0", "", "none", "off", "-" => 0 (skip)
    """
    if tok is None:
        return 0
    s = str(tok).strip().lower().rstrip('+')
    if s in ("", "0", "none", "off", "-", "skip"):
        return 0
    # bare integer seconds?
    if re.fullmatch(r"\d+", s):
        return int(s)
    # h/m/s pattern
    m = re.fullmatch(r"(?:(\d+)h)?(?:(\d+)m)?(?:(\d+)s)?", s)
    if not m:
        # fallback: try plain int
        try:
            return int(float(s))
        except Exception:
            return 0
    h = int(m.group(1) or 0)
    m_ = int(m.group(2) or 0)
    s_ = int(m.group(3) or 0)
    return h*3600 + m_*60 + s_

def get_flag(argv, name, default="0"):
    """
    Return the value of a --name=... flag if present, else default.
    """
    prefix = f"--{name}="
    for a in argv:
        if a.startswith(prefix):
            return a[len(prefix):]
    return default

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
# 2.5 Heavy computations in a child process (real timeouts)
############################

def _dim_worker(infile, outq):
    try:
        vars, F, _, polys_str = parse_input_system(infile)
        R = PolynomialRing(F, vars)
        polynomials = [R(s) for s in polys_str]
        I = R.ideal(polynomials)
        outq.put(("ok", int(I.dimension())))
    except Exception as e:
        outq.put(("err", str(e)))

def _hilbert_worker(infile, outq):
    try:
        vars, F, _, polys_str = parse_input_system(infile)
        R = PolynomialRing(F, vars)
        polynomials = [R(s) for s in polys_str]
        I = R.ideal(polynomials)
        hilb = I.hilbert_polynomial()
        # For zero-dimensional ideals, hilb(0) equals the number of solutions (with multiplicity)
        if I.dimension() == 0:
            degI = hilb(0)
            # Ensure it's an int if possible
            try:
                degI = int(degI)
            except Exception:
                pass
        else:
            degI = hilb.degree()
        outq.put(("ok", degI))
    except Exception as e:
        outq.put(("err", str(e)))

def run_with_timeout(target, args, timeout_sec):
    """
    Run 'target(*args)' in a separate process and enforce a wall-clock timeout.
    Returns (status, value_or_message), where status ∈ {"ok","timeout","err"}.
    """
    if timeout_sec <= 0:
        # treated as "skip"
        return ("skip", "skipped")

    # Prefer 'fork' on Linux; fall back to 'spawn' elsewhere.
    try:
        ctx = multiprocessing.get_context("fork")
    except Exception:
        ctx = multiprocessing.get_context("spawn")

    q = ctx.Queue()
    p = ctx.Process(target=target, args=(*args, q))
    p.start()
    p.join(timeout=timeout_sec)
    if p.is_alive():
        p.terminate()
        p.join()
        return ("timeout", f"timed out after {timeout_sec}s")
    # Read result
    try:
        status, value = q.get_nowait()
    except Exception:
        return ("err", "no result from worker")
    if status == "ok":
        return ("ok", value)
    else:
        return ("err", value)

############################
# 3. Main diagnostic logic
############################
def main():
    if len(sys.argv) < 2:
        print("Usage: sage system_diagnosis.sage <input_file.in> [--dim=SECONDS] [--hilbert=SECONDS]")
        sys.exit(1)

    # Positional: first non-flag argument is the file
    infile = None
    for a in sys.argv[1:]:
        if not a.startswith("--"):
            infile = a
            break
    if infile is None:
        print("Error: missing input file.")
        sys.exit(1)

    # Budgets
    dim_budget     = parse_time_budget(get_flag(sys.argv[1:], "dim", "0"))
    hilbert_budget = parse_time_budget(get_flag(sys.argv[1:], "hilbert", "0"))

    outlog = infile.replace(".in", "_diagnostics.log").replace("data/", "logs/")

    # Parse the input file (cheap) for all the quick stats
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
    if dim_budget > 0:
        status, value = run_with_timeout(_dim_worker, (infile,), dim_budget)
        if status == "ok":
            dim = int(value)
            log.append(f"Krull dimension of the ideal: {dim}")
            if dim == 0:
                log.append("System is zero-dimensional (finite number of solutions).")
            elif dim == -1:
                log.append("Ideal is the unit ideal (1 in the basis).")
            else:
                log.append("System is positive-dimensional (infinitely many solutions).")
        elif status == "timeout":
            log.append(f"Could not compute dimension: {value}")
        else:
            log.append(f"Could not compute dimension: {value}")
    else:
        log.append("Krull dimension computation skipped (no time budget).")

    # --- Degree of the ideal (Hilbert polynomial) ---
    if hilbert_budget > 0:
        status, value = run_with_timeout(_hilbert_worker, (infile,), hilbert_budget)
        if status == "ok":
            degI = value
            log.append(f"Degree of the ideal (zero-dim: number of points, else: degree of Hilbert polynomial): {degI}")
        elif status == "timeout":
            log.append(f"Could not compute ideal degree (Hilbert polynomial): {value}")
        else:
            log.append(f"Could not compute ideal degree (Hilbert polynomial): {value}")
    else:
        log.append("Degree (Hilbert polynomial) computation skipped (no time budget).")

    # No shape position/lex basis attempted!
    log.append("Shape position check skipped (lex Groebner basis not computed).")

    # --- Special forms: Boolean and quadratic ---
    is_boolean = F.order() == 2 and all(p.total_degree() <= 2 for p in polynomials)
    log.append(f"System is Boolean quadratic? {'Yes' if is_boolean else 'No'}")
    is_quadratic = all(p.total_degree() <= 2 for p in polynomials)
    log.append(f"System is quadratic (all degree ≤2)? {'Yes' if is_quadratic else 'No'}")

    # Save log (and echo to stdout)
    os.makedirs(os.path.dirname(outlog) or ".", exist_ok=True)
    with open(outlog, "w") as f:
        for line in log:
            print(line)
            f.write(line + "\n")
    print(f"\nDiagnostics saved to {outlog}")

main()
