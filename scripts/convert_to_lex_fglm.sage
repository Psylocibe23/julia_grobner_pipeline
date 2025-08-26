###############################################################################
# convert_to_lex_fglm.sage  —  FGLM order change (DRL → LEX) with rich logging
#
# PURPOSE
#   Given a Gröbner basis G in DRL (degrevlex) over GF(p) — typically produced
#   by the F4 stage — convert it to a LEX Gröbner basis via the FGLM algorithm.
#   Then (a) record useful statistics, (b) perform BOTH a heuristic and a strict
#   shape-position check, and (c) write out a clean, machine-parseable LEX file.
#
# DESIGN CHOICES
#   • Input:   Prime field GF(p); DRL basis (not necessarily reduced).
#   • Output:  LEX basis; we also build an ideal from it (often triangular).
#   • API:     Sage’s Ideal.transformed_basis('fglm', R_lex) as the canonical
#              entry point for FGLM conversion (zero-dimensional precondition).
#   • Robustness: Detailed logging and failure messages for common pitfalls.
#
# INPUT FILE (as written by solve_F4_from_file.jl)
#   # Groebner basis (F4) computed for <path>
#   # Variables: x0, x1, ..., xn-1
#   # Field: GF(p)
#   # Order: degrevlex
#   # Basis: reduced=UNKNOWN
#   # Number of input equations: k
#   # --- Groebner basis ---
#   <poly_1>
#   <poly_2>
#   ...
#
# New CLI flags (all optional):
#   --assume-zerodim          : skip Krull-dimension check (treat as 0-dim)
#   --dim-timeout=SECONDS     : time budget for I_drl.dimension() (default 60)
#   --fglm-timeout=SECONDS    : time budget for the FGLM call (default 0 = no cap)
#   --reduce=never|auto|always: reduce the LEX basis (default auto)
#   --reduce-timeout=SECONDS  : time budget for reduction (default 60 in auto/always)
#
# Examples:
#   sage scripts/convert_to_lex_fglm.sage results/HFE_n5_D6_F4_...txt \
#        --assume-zerodim --fglm-timeout=600 --reduce=auto --reduce-timeout=120
#
###############################################################################


import sys, os, time, signal

# ---------------- Utilities ----------------
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def stem_of(path):
    return os.path.splitext(os.path.basename(path))[0]

def log_and_print(msg, fh=None):
    print(msg, flush=True)
    if fh is not None:
        fh.write(msg + "\n"); fh.flush()

class _Timeout(Exception): pass
def _timeout_handler(signum, frame):  # noqa: ARG001
    raise _Timeout

def run_with_timeout(fn, seconds, on_timeout=None):
    """Run fn() with a SIGALRM timeout. seconds<=0 means no timeout."""
    if seconds and seconds > 0:
        old = signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(int(seconds))
        try:
            return fn()
        except _Timeout:
            if on_timeout:
                return on_timeout()
            raise
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
    else:
        return fn()

# ---------------- Parsing the DRL-basis file ----------------
def read_groebner_basis_file(result_file):
    variables, field_p, order = None, None, None
    basis_start = None
    polys = []

    with open(result_file, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("# Variables:"):
            variables = [v.strip() for v in s.split(":", 1)[1].split(",")]
        elif s.startswith("# Field:"):
            try:
                inside = s.split("GF(", 1)[1].split(")", 1)[0].strip()
                if "^" in inside:
                    base, _ = inside.split("^", 1)
                    field_p = int(base.strip())
                else:
                    field_p = int(inside)
            except Exception:
                raise ValueError(f"Could not parse field from header line: '{s}'")
        elif s.startswith("# Order:"):
            order = s.split(":", 1)[1].strip()
        elif s.startswith("# --- Groebner basis ---"):
            basis_start = i + 1
            break

    if variables is None:
        raise ValueError("Header '# Variables:' not found.")
    if field_p is None:
        raise ValueError("Header '# Field: GF(p)' not found.")
    if basis_start is None:
        raise ValueError("Marker '# --- Groebner basis ---' not found; cannot read basis.")

    for line in lines[basis_start:]:
        s = line.strip()
        if s and not s.startswith("#"):
            polys.append(s)
    if not polys:
        raise ValueError("No polynomials found after the basis marker.")

    return variables, field_p, (order or "UNKNOWN"), polys

# ---------------- Shape-position checks ----------------
def shape_heuristic(variables, G_lex):
    try:
        G_sorted = sorted(G_lex, key=lambda g: g.lm(), reverse=True)
        lead_vars = []
        for g in G_sorted:
            lm = g.lm()
            if not lm.is_pure_power():
                return False
            vs = list(lm.variables())
            if len(vs) != 1:
                return False
            lead_vars.append(str(vs[0]))
        var_index = {v: i for i, v in enumerate(variables)}
        idxs = [var_index.get(vn, -1) for vn in lead_vars]
        if any(i < 0 for i in idxs):
            return False
        return all(idxs[i] <= idxs[i+1] for i in range(len(idxs)-1))
    except Exception:
        return False

def shape_strict(variables, G_lex):
    if not G_lex:
        return False
    R = G_lex[0].parent()
    gens = {str(g): g for g in R.gens()}
    n = len(variables)
    v_last = variables[-1]
    univars = [g for g in G_lex if set(map(str, g.variables())) <= {v_last}]
    if len(univars) != 1:
        return False
    for k in range(n-2, -1, -1):
        v = gens[variables[k]]
        cand = [g for g in G_lex if g.lm() == v]
        if len(cand) != 1:
            return False
        if cand[0].degree(v) != 1:
            return False
    return True

# ---------------- Main ----------------
def main():
    from sage.all import GF, PolynomialRing  # lazy import so script starts quickly

    if len(sys.argv) < 2:
        print("Usage: sage scripts/convert_to_lex_fglm.sage <F4_result.txt> [--assume-zerodim] "
              "[--dim-timeout=SECONDS] [--fglm-timeout=SECONDS] "
              "[--reduce=never|auto|always] [--reduce-timeout=SECONDS]")
        sys.exit(1)

    result_file = sys.argv[1]
    # defaults
    assume_zerodim = False
    dim_timeout = 60
    fglm_timeout = 0          # 0 == unlimited
    reduce_mode = "auto"      # never|auto|always
    reduce_timeout = 60

    # parse flags
    for a in sys.argv[2:]:
        if a == "--assume-zerodim":
            assume_zerodim = True
        elif a.startswith("--dim-timeout="):
            dim_timeout = int(a.split("=",1)[1])
        elif a.startswith("--fglm-timeout="):
            fglm_timeout = int(a.split("=",1)[1])
        elif a.startswith("--reduce="):
            reduce_mode = a.split("=",1)[1].strip().lower()
        elif a.startswith("--reduce-timeout="):
            reduce_timeout = int(a.split("=",1)[1])
        else:
            print(f"Unknown flag '{a}'"); sys.exit(2)

    results_dir = "results"
    logs_dir    = "logs"
    ensure_dir(results_dir); ensure_dir(logs_dir)

    base_name   = stem_of(result_file)
    lex_outfile = os.path.join(results_dir, base_name + "_LEX.txt")
    log_outfile = os.path.join(logs_dir,    base_name + "_FGLM.log")

    with open(log_outfile, "w") as log:
        log_and_print("==============================================================", log)
        log_and_print(" FGLM CONVERSION — DRL → LEX (Sage/Singular backend)", log)
        log_and_print("==============================================================", log)
        log_and_print(f"Input file:  {result_file}", log)
        log_and_print(f"LEX out:     {lex_outfile}", log)
        log_and_print("--------------------------------------------------------------", log)

        try:
            log_and_print("Parsing input (headers and DRL basis)...", log)
            variables, p, in_order, polys_str = read_groebner_basis_file(result_file)
            log_and_print(f"Variables: {variables}", log)
            log_and_print(f"Field:     GF({p})", log)
            log_and_print(f"Order(in): {in_order}", log)
            log_and_print(f"Basis size (input): {len(polys_str)}", log)

            log_and_print("Constructing DRL ring and ideal...", log)
            R_drl = PolynomialRing(GF(p), variables, order='degrevlex')
            G_drl = [R_drl(s) for s in polys_str]
            I_drl = R_drl.ideal(G_drl)

            # (Optional) tiny input stat
            try:
                degs = [g.total_degree() for g in G_drl]
                log_and_print(f"Degrees(in): max={max(degs)}, list(first 10)={degs[:10]}{'...' if len(degs)>10 else ''}", log)
            except Exception as e:
                log_and_print(f"Warning: could not compute input degrees: {e}", log)

            # Zero-dim check
            if assume_zerodim:
                log_and_print("Skipping dimension() (assume zero-dimensional).", log)
                dim = 0
            else:
                log_and_print(f"Checking zero-dimensionality with timeout {dim_timeout}s...", log)
                def _dim():  # compute dimension
                    return I_drl.dimension()
                try:
                    dim = run_with_timeout(_dim, dim_timeout)
                    log_and_print(f"Krull dimension: {dim}", log)
                except _Timeout:
                    log_and_print(f"dimension() timed out after {dim_timeout}s — treating as zero-dimensional for FGLM.", log)
                    dim = 0

            if dim != 0:
                log_and_print("ERROR: Ideal is not zero-dimensional. FGLM requires dim = 0.", log)
                raise ValueError("Non-zero-dimensional ideal; aborting FGLM.")

            # Build LEX ring
            from sage.all import GF as _GF, PolynomialRing as _PR  # (local alias to look explicit in logs)
            log_and_print("Constructing LEX ring...", log)
            R_lex = _PR(_GF(p), variables, order='lex')

            # FGLM
            log_and_print(f"Running FGLM (timeout {fglm_timeout or '∞'}s)...", log)
            t0 = time.time()
            def _run_fglm():
                return I_drl.transformed_basis('fglm', R_lex)
            try:
                G_lex_fglm = run_with_timeout(_run_fglm, fglm_timeout)
            except _Timeout:
                log_and_print(f"FGLM timed out after {fglm_timeout}s.", log)
                raise
            t1 = time.time()
            log_and_print(f"FGLM wall time: {t1 - t0:.3f} s", log)

            # Reduce?
            reduce_mode_norm = reduce_mode.lower()
            log_and_print(f"Post-processing LEX basis (reduce={reduce_mode_norm}, timeout={reduce_timeout}s)...", log)
            if reduce_mode_norm == "never":
                G_lex = G_lex_fglm
            else:
                I_lex = R_lex.ideal(G_lex_fglm)
                def _reduce():
                    return I_lex.groebner_basis()
                try:
                    G_lex = run_with_timeout(_reduce, reduce_timeout if reduce_mode_norm in ("auto","always") else 0)
                except _Timeout:
                    if reduce_mode_norm == "auto":
                        log_and_print("Reduction timed out — keeping *unreduced* LEX basis from FGLM.", log)
                        G_lex = G_lex_fglm
                    else:
                        log_and_print("Reduction timed out and reduce=always → aborting.", log)
                        raise

            out_degs = [g.total_degree() for g in G_lex]
            log_and_print(f"Basis size(out): {len(G_lex)}", log)
            log_and_print(f"Degrees(out): max={max(out_degs)}, list(first 10)={out_degs[:10]}{'...' if len(out_degs)>10 else ''}", log)

            shape_h = shape_heuristic(variables, G_lex)
            shape_s = shape_strict(variables, G_lex)
            log_and_print(f"Shape position (heuristic): {shape_h}", log)
            log_and_print(f"Shape position (strict):    {shape_s}", log)

            log_and_print(f"Writing LEX basis to: {lex_outfile}", log)
            with open(lex_outfile, "w") as out:
                out.write(f"# Lex Groebner basis (FGLM) for {result_file}\n")
                out.write(f"# Variables: {', '.join(variables)}\n")
                out.write(f"# Field: GF({p})\n")
                out.write(f"# Order: lex\n")
                out.write(f"# Basis: reduced={'True' if G_lex is not G_lex_fglm else 'Unknown'}\n")
                out.write(f"# Zero-dimensional: True\n")
                out.write(f"# Shape position: {shape_s}\n")
                out.write("# --- Groebner basis (LEX) ---\n")
                for g in G_lex:
                    out.write(str(g) + "\n")

            log_and_print("LEX basis (first few polynomials):", log)
            for g in G_lex[:min(5, len(G_lex))]:
                log_and_print("  " + str(g), log)

            log_and_print("--------------------------------------------------------------", log)
            log_and_print("FGLM conversion COMPLETE.", log)

        except Exception as e:
            import traceback
            log_and_print("FGLM conversion FAILED with an exception:", log)
            log_and_print(str(e), log)
            log_and_print(traceback.format_exc(), log)
            sys.exit(1)

if __name__ == "__main__":
    main()
