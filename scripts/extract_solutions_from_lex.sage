###############################################################################
# extract_solutions_from_lex.sage  —  Solution extraction from a LEX GB
#
# PURPOSE
#   Given a Gröbner basis G in LEX order over GF(p) extract all solutions in the base field GF(p).
#
# MAIN
#   1) Parse headers and the LEX basis polynomials from the input file.
#   2) Build R = GF(p)[variables] with 'lex' order and form the ideal I.
#   3) Fast route for scale:
#        • Look for a univariate in the last variable; if found, solve it in F.
#        • Else (default) just try all values of the last variable in F
#          (for GF(2), {0,1}), and branch upward with pruning.
#        • At each step:
#            - substitute assigned vars,
#            - if any polynomial becomes a nonzero constant then prune branch,
#            - otherwise solve next variable.
#   4) Verify each candidate against the full LEX basis and deduplicate.
#   5) Optional fallback (SMALL systems only): I.variety(GF(p)).
#
# CLI
#   sage scripts/extract_solutions_from_lex.sage results/<stem>_LEX.txt \
#        [--first-only] [--max-solutions N] [--node-limit N] \
#        [--no-fallback] [--check-dim] [--elim]
#
#   --check-dim : force computing I.dimension() (SLOW / can overflow; default OFF)
#   --elim : allow elimination-ideal fallback for last variable (default OFF)
#
# OUTPUT
#   results/<stem>_LEX_sols.txt : one solution per line, e.g. {x0: 1, x1: 0}
#   logs/<stem>_LEX_extract.log : verbose run log
#
###############################################################################

import sys, os, time
from sage.all import GF, PolynomialRing

# ================= Utilities ================
def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def stem_of(path):
    return os.path.splitext(os.path.basename(path))[0]

def log_and_print(msg, fh=None):
    print(msg, flush=True)
    if fh is not None:
        fh.write(msg + "\n")
        fh.flush()
# ================ 1: parse file ================
def read_lex_basis_file(lex_file):
    """
    Read headers and polynomials from a LEX GB file produced by fglm_adjusted.sage.

    Returns: (variables, p, order, polys_as_strings, shape_hdr_or_None, zerodim_hdr_or_None)
    """
    variables, p, order, shape_hdr, zerodim_hdr = None, None, None, None, None
    polys_str = []
    basis_start = None

    with open(lex_file, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("# Variables:"):
            variables = [v.strip() for v in s.split(":", 1)[1].split(",")]
        elif s.startswith("# Field:"):
            # Expect "# Field: GF(p)" (accepts GF(p^1) but takes base p)
            try:
                inside = s.split("GF(", 1)[1].split(")", 1)[0].strip()
                if "^" in inside:
                    base, _ = inside.split("^", 1)
                    p = int(base.strip())
                else:
                    p = int(inside)
            except Exception:
                raise ValueError(f"Could not parse field from header line: '{s}'")
        elif s.startswith("# Order:"):
            order = s.split(":", 1)[1].strip()
        elif s.startswith("# Shape position:"):
            txt = s.split(":", 1)[1].strip()
            shape_hdr = (txt.lower() == "true")
        elif s.startswith("# Zero-dimensional:"):
            txt = s.split(":", 1)[1].strip()
            zerodim_hdr = (txt.lower() == "true")
        elif s.startswith("# --- Groebner basis"):
            basis_start = i + 1
            break

    if variables is None:
        raise ValueError("Header '# Variables:' not found.")
    if p is None:
        raise ValueError("Header '# Field: GF(p)' not found.")
    if order is None:
        order = "UNKNOWN"
    if basis_start is None:
        raise ValueError("Marker '# --- Groebner basis (LEX) ---' not found.")

    for line in lines[basis_start:]:
        s = line.strip()
        if s and not s.startswith("#"):
            polys_str.append(s)
    if not polys_str:
        raise ValueError("No polynomials found after the basis marker.")

    return variables, p, order, polys_str, shape_hdr, zerodim_hdr

# ================ 2: shape position check ================
def lex_shape_position(variables, G_lex):
    """
    Heuristic 'shape' test: leading monomials should be pure powers, and
    the sequence of leading variables should be nondecreasing in lex order.
    (Sufficient but not necessary.)
    """
    try:
        G_sorted = sorted(G_lex, key=lambda g: g.lm(), reverse=True)
        lead_vars = []
        for g in G_sorted:
            lm = g.lm()
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

# ================ 3: root finding ================
def to_true_univariate(poly_mv, var, R, F):
    """
    Convert a multivariate polynomial `poly_mv` into a true univariate f(T) in F[T] by reading its
    exponent tuples and coefficients.
    """
    gens = R.gens()
    try:
        idx = list(gens).index(var)
    except ValueError:
        raise TypeError("Requested variable is not a generator of the ring.")

    exps = poly_mv.exponents()       # list of tuples (e_0, ..., e_{n-1})
    coeffs = poly_mv.coefficients()  # list of coefficients, same order as exps
    if len(exps) != len(coeffs):
        raise RuntimeError("Inconsistent monomial/coeff arrays in polynomial.")

    maxdeg = 0
    for e in exps:
        other = sum(e[j] for j in range(len(e)) if j != idx)
        if other != 0:
            raise TypeError("Polynomial is not univariate in the requested variable.")
        if e[idx] > maxdeg:
            maxdeg = e[idx]

    U.<T> = PolynomialRing(F)
    arr = [F(0)] * (maxdeg + 1)
    for e, c in zip(exps, coeffs):
        arr[e[idx]] += F(c)
    return U(arr)

def roots_over_field_univariate(poly_mv, var, R, F):
    """
    Find roots in the base field F for a polynomial that is univariate in `var`.
    """
    f_uni = to_true_univariate(poly_mv, var, R, F)      # f ∈ F[T]
    return [a for (a, mult) in f_uni.roots(multiplicities=True, ring=F)]

# ================ 4: branch triangular solver ================
def solve_shape_lex_branch(R, variables, G_lex, F, log=None,
                           allow_elimination=False,
                           first_only=False, max_solutions=None, max_nodes=10_000_000):
    """
    Branching solver for triangular LEX bases.

    • Preferred start: pick a univariate in the last variable x_{n-1} if present.
      If none and allow_elimination=True, try an elimination ideal (may be heavy).
      Otherwise, simply iterate over all a in F for x_{n-1} (cheap for GF(2)).
    • pruning:
        - if any g in G_lex reduces to a nonzero constant, prune.
        - if step polynomial for vk is univariate: solve over F (linear shortcut if deg=1).
        - if not univariate: for GF(2) try {0,1}; for small base fields iterate small set.

    Returns a list of dicts {R.gen(i): value} representing solutions.
    """
    n = len(variables)
    gens = R.gens()
    name_to_var = {str(gens[i]): gens[i] for i in range(n)}

    v_last = name_to_var[variables[-1]]
    # Prefer a univariate already present
    univariates = [g for g in G_lex if set(map(str, g.variables())) <= {variables[-1]}]

    roots_last = None
    if univariates:
        g_last = max(univariates, key=lambda g: g.degree(v_last))
        roots_last = roots_over_field_univariate(g_last, v_last, R, F)
        log_and_print(f"Eliminant in {variables[-1]} has {len(roots_last)} root(s): {roots_last}", log)
    elif allow_elimination:
        log_and_print("No univariate in last variable; attempting elimination ideal (may be heavy)...", log)
        try:
            El = R.ideal(G_lex).elimination_ideal(gens[:-1])
            cand = [g for g in El.gens() if set(map(str, g.variables())) <= {variables[-1]}]
            if cand:
                g_last = max(cand, key=lambda g: g.degree(v_last))
                roots_last = roots_over_field_univariate(g_last, v_last, R, F)
                log_and_print(f"[elim] Eliminant in {variables[-1]} has {len(roots_last)} root(s): {roots_last}", log)
        except Exception as e:
            log_and_print(f"[elim] Elimination failed ({e}); falling back to brute-force over base field.", log)

    if roots_last is None:
        # Cheap & robust fallback: try all values in the base field
        # (for GF(2) this is just {0,1})
        roots_last = list(F)
        log_and_print(f"No univariate eliminant available. Trying all {len(roots_last)} value(s) in GF({F.order()}) for {variables[-1]}.", log)

    # Representative polynomial per leading variable 
    lead_poly = {}
    for g in G_lex:
        lm = g.lm()
        vs = list(lm.variables())
        if len(vs) == 1:
            vn = str(vs[0])
            lead_poly.setdefault(vn, g)

    # Pruning: if any g becomes a nonzero constant under current partial assignment then prune
    def prunable(assign_dict):
        Fbase = R.base_ring()
        for g in G_lex:
            h = g.subs(assign_dict)
            # Polynomial case
            if hasattr(h, "is_constant"):
                if h.is_constant() and h != 0:
                    return True
            else:
                # Base-field scalar (or python int)
                if h != 0:
                    return True
        return False

    solutions = []
    nodes = 0

    for r in roots_last:
        start = {v_last: r}
        if prunable(start):
            continue
        stack = [(n-2, start)]

        while stack:
            if max_nodes is not None and nodes >= max_nodes:
                log_and_print(f"[WARN] node limit {max_nodes} reached; stopping search.", log)
                return solutions
            k, assign = stack.pop()
            nodes += 1

            if k < 0:
                if all(g.subs(assign) == 0 for g in G_lex):
                    solutions.append(dict(assign))
                    if first_only:
                        return solutions
                    if max_solutions is not None and len(solutions) >= max_solutions:
                        return solutions
                continue

            vk_name = variables[k]
            vk = name_to_var[vk_name]

            poly = lead_poly.get(vk_name, None)
            if poly is None:
                cands = [g for g in G_lex if vk in g.variables()]
                if not cands:
                    vals = list(F) if F.order() <= 3 else [F(0), F(1)]
                    for val in vals:
                        nxt = dict(assign); nxt[vk] = val
                        if not prunable(nxt):
                            stack.append((k-1, nxt))
                    continue
                poly = cands[0]

            pk = poly.subs(assign)
            if pk == 0:
                stack.append((k-1, assign))
                continue

            if hasattr(pk, "is_constant"):
                if pk.is_constant():
                    stack.append((k-1, assign))
                    continue
            else:
                stack.append((k-1, assign))
                continue

            # Try univariate solve in vk
            cand_vals = None
            try:
                fk = to_true_univariate(pk, vk, R, F)
            except (TypeError, ValueError):
                # Not univariate: iterate small set
                cand_vals = list(F) if F.order() <= 3 else [F(0), F(1)]
            else:
                if fk.degree() == 1:
                    L = fk.list()
                    b = L[0] if len(L) > 0 else F(0)
                    a = L[1] if len(L) > 1 else F(0)
                    if a != 0:
                        cand_vals = [(-b)/a]
                    else:
                        cand_vals = [a for (a, _) in fk.roots(multiplicities=True, ring=F)]
                else:
                    cand_vals = [a for (a, _) in fk.roots(multiplicities=True, ring=F)]

            if not cand_vals:
                cand_vals = list(F) if F.order() <= 3 else [F(0), F(1)]

            for val in cand_vals:
                nxt = dict(assign); nxt[vk] = val
                if not prunable(nxt):
                    stack.append((k-1, nxt))

    # Deduplicate by tuple of ints
    uniq, seen = [], set()
    name_to_var = {str(g): g for g in gens}
    for sol in solutions:
        key = tuple(int(sol[name_to_var[v]]) for v in variables)
        if key not in seen:
            seen.add(key)
            uniq.append(sol)
    return uniq

def pretty_value(a):
    try:
        return str(int(a))
    except Exception:
        return str(a)

def write_solutions(outfile, variables, solutions, name_to_var):
    with open(outfile, "w") as out:
        out.write(f"# Solutions extracted from LEX basis\n")
        out.write(f"# Variables: {', '.join(variables)}\n")
        out.write(f"# Count: {len(solutions)}\n")
        for sol in solutions:
            parts = []
            for vname in variables:
                rv = name_to_var[vname]
                parts.append(f"{vname}: {pretty_value(sol[rv])}")
            out.write("{" + ", ".join(parts) + "}\n")

# ================ Main ================
def main():
    # Args
    if len(sys.argv) < 2:
        print("Usage: sage scripts/extract_solutions_from_lex.sage results/<stem>_LEX.txt "
              "[--first-only] [--max-solutions N] [--node-limit N] "
              "[--no-fallback] [--check-dim] [--elim]")
        sys.exit(1)

    # Positional
    lex_file = sys.argv[1]
    if not os.path.isfile(lex_file):
        print(f"Input file not found: {lex_file}")
        sys.exit(2)

    # Flags
    first_only = False
    max_solutions = None
    node_limit = 10_000_000
    no_fallback = False
    check_dim = False # OFF by default (prevents stdhilb overflow)
    allow_elim = False # OFF by default (avoid heavy elimination GB)

    i = 2
    while i < len(sys.argv):
        a = sys.argv[i]
        if a == "--first-only":
            first_only = True; i += 1
        elif a.startswith("--max-solutions"):
            if "=" in a:  max_solutions = int(a.split("=",1)[1]); i += 1
            else:         max_solutions = int(sys.argv[i+1]); i += 2
        elif a.startswith("--node-limit"):
            if "=" in a:  node_limit = int(a.split("=",1)[1]); i += 1
            else:         node_limit = int(sys.argv[i+1]); i += 2
        elif a == "--no-fallback":
            no_fallback = True; i += 1
        elif a == "--check-dim":
            check_dim = True; i += 1
        elif a == "--elim":
            allow_elim = True; i += 1
        else:
            print(f"Unrecognized argument: {a}")
            sys.exit(3)

    # Paths
    results_dir = "results"
    logs_dir = "logs"
    ensure_dir(results_dir)
    ensure_dir(logs_dir)

    base_name = stem_of(lex_file)
    sols_out = os.path.join(results_dir, base_name + "_sols.txt")
    log_out = os.path.join(logs_dir,    base_name + "_LEX_extract.log")

    with open(log_out, "w") as log:
        log_and_print("==============================================================", log)
        log_and_print(" LEX SOLUTION EXTRACTION — GF(p) (Sage)", log)
        log_and_print("==============================================================", log)
        log_and_print(f"Input LEX basis: {lex_file}", log)
        log_and_print(f"Output solutions: {sols_out}", log)
        log_and_print("--------------------------------------------------------------", log)

        try:
            # Parse input
            variables, p, order, polys_str, shape_hdr, zerodim_hdr = read_lex_basis_file(lex_file)
            log_and_print(f"Variables: {variables}", log)
            log_and_print(f"Field: GF({p})", log)
            log_and_print(f"Order: {order}", log)
            log_and_print(f"Basis size: {len(polys_str)}", log)
            if shape_hdr is not None:
                log_and_print(f"Header says shape position: {shape_hdr}", log)
            if zerodim_hdr is not None:
                log_and_print(f"Header says zero-dimensional: {zerodim_hdr}", log)

            # Build LEX ring/ideal
            F = GF(p)
            R = PolynomialRing(F, variables, order='lex')
            G_lex = [R(s) for s in polys_str]
            I_lex = R.ideal(G_lex)
            gens = R.gens()
            name_to_var = {str(gens[i]): gens[i] for i in range(len(gens))}

            # (Optional) zero-dimensionality check (OFF by default for scale)
            if check_dim:
                try:
                    dim = I_lex.dimension()
                    log_and_print(f"Krull dimension (requested): {dim}", log)
                except Exception as e:
                    log_and_print(f"[WARN] dimension() failed with: {e}. Proceeding without it.", log)
            else:
                log_and_print("Skipping dimension() (default). Use --check-dim to force it (may be slow/unstable).", log)

            # Heuristic shape info
            in_shape = lex_shape_position(variables, G_lex)
            log_and_print(f"Shape position (heuristic): {in_shape}", log)

            # Branching triangular solver
            t0 = time.time()
            solutions = solve_shape_lex_branch(
                R, variables, G_lex, F, log=log,
                allow_elimination=allow_elim,
                first_only=first_only, max_solutions=max_solutions, max_nodes=node_limit
            )
            t1 = time.time()
            log_and_print(f"Branching triangular solver produced {len(solutions)} solution(s) in {t1 - t0:.6f} s.", log)

            # fallback: full enumeration (SMALL systems only)
            if not solutions and not no_fallback:
                if len(variables) <= 16:
                    log_and_print("Falling back to I_lex.variety(GF(p)) enumeration (n ≤ 16)...", log)
                    t1 = time.time()
                    sols_dicts = I_lex.variety(F)
                    t2 = time.time()
                    log_and_print(f"variety(GF({p})) returned {len(sols_dicts)} solution(s) in {t2 - t1:.6f} s.", log)
                    solutions = []
                    for d in sols_dicts:
                        sol = {}
                        for vname in variables:
                            rv = name_to_var[vname]
                            sol[rv] = d[rv]
                        solutions.append(sol)
                else:
                    log_and_print("No solutions from branching solver; skipping variety() (n is large). Use --no-fallback to silence.", log)
            elif not solutions and no_fallback:
                log_and_print("No solutions found by branching solver; fallback disabled (--no-fallback).", log)

            # Final verification
            log_and_print("Verifying all solutions against the LEX basis...", log)
            verified = []
            for sol in solutions:
                if all(g.subs(sol) == 0 for g in G_lex):
                    verified.append(sol)
            if len(verified) != len(solutions):
                log_and_print(f"Discarded {len(solutions) - len(verified)} non-solutions after verification.", log)
            solutions = verified

            # Deduplicate and write
            uniq, seen = [], set()
            for sol in solutions:
                key = tuple(int(sol[name_to_var[v]]) for v in variables)
                if key not in seen:
                    seen.add(key)
                    uniq.append(sol)
            solutions = uniq

            write_solutions(sols_out, variables, solutions, name_to_var)
            log_and_print(f"Wrote {len(solutions)} solution(s) to {sols_out}", log)

            for i, sol in enumerate(solutions[:min(5, len(solutions))]):
                pretty = ", ".join(f"{v}: {pretty_value(sol[name_to_var[v]])}" for v in variables)
                log_and_print(f"Solution #{i+1}: {{{pretty}}}", log)

            log_and_print("--------------------------------------------------------------", log)
            log_and_print("LEX solution extraction COMPLETE.", log)

        except Exception as e:
            import traceback
            log_and_print("Extraction FAILED with an exception:", log)
            log_and_print(str(e), log)
            log_and_print(traceback.format_exc(), log)
            sys.exit(1)

if __name__ == "__main__":
    main()
