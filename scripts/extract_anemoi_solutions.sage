###############################################################################
# extract_anemoi_solutions.sage
#
# Given a LEX Gröbner basis G over GF(p) (from FGLM), this script:
#   • Finds the univariate "eliminant" in the last variable.
#   • Reports counts of solutions over the algebraic closure
#       - total (with multiplicity): deg(f_last)
#       - distinct (set-wise): deg(squarefree(f_last))
#       - also breaks down by extension degrees via factorization over F_p
#   • Optionally enumerates actual solutions in minimal extensions GF(p^d)
#     for small factor degrees 
#
# Output:
#   results/<stem>_AC_report.txt (counts + factor breakdown)
#   results/<stem>_AC_solutions.txt (optional: enumerated solutions)
#   logs/<stem>_AC_extract.log (detailed log)
###############################################################################

import sys, os, time

# ================ Helpers ================
def ensure_dir(d): os.path.isdir(d) or os.makedirs(d)
def stem_of(p): return os.path.splitext(os.path.basename(p))[0]
def logprint(msg, fh=None):
    print(msg, flush=True)
    if fh: fh.write(msg + "\n"); fh.flush()

# ================ Parsing LEX basis file ================
def read_lex_basis_file(path):
    """
    Expect the header your FGLM scripts write:
      # Variables: v0, v1, ..., v_{n-1}
      # Field: GF(p)
      # Order: lex
      # (optional) Shape position: True|False
      # (optional) Zero-dimensional: True|False
      # --- Groebner basis (LEX) ---
      <poly_1>
      ...
    """
    variables, p, order = None, None, None
    shape_hdr, zerodim_hdr = None, None
    polys_str, start = [], None

    with open(path, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("# Variables:"):
            variables = [v.strip() for v in s.split(":",1)[1].split(",")]
        elif s.startswith("# Field:"):
            inside = s.split("GF(",1)[1].split(")",1)[0].strip()
            p = int(inside.split("^",1)[0])
        elif s.startswith("# Order:"):
            order = s.split(":",1)[1].strip()
        elif s.startswith("# Shape position:"):
            shape_hdr = (s.split(":",1)[1].strip().lower()=="true")
        elif s.startswith("# Zero-dimensional:"):
            zerodim_hdr = (s.split(":",1)[1].strip().lower()=="true")
        elif s.startswith("# --- Groebner basis"):
            start = i+1
            break

    if variables is None: raise ValueError("Missing '# Variables:' header.")
    if p is None: raise ValueError("Missing '# Field: GF(p)' header.")
    if start is None: raise ValueError("Missing basis marker.")

    for line in lines[start:]:
        s = line.strip()
        if s and not s.startswith("#"):
            polys_str.append(s)
    if not polys_str: raise ValueError("No polynomials after basis marker.")

    return variables, p, (order or "UNKNOWN"), polys_str, shape_hdr, zerodim_hdr

# ================ Shape position tests ================
def shape_heuristic(variables, G):
    try:
        Gs = sorted(G, key=lambda g: g.lm(), reverse=True)
        idx = {v:i for i,v in enumerate(variables)}
        prev = -1
        for g in Gs:
            lm = g.lm()
            vs = list(lm.variables())
            if len(vs)!=1: return False
            cur = idx[str(vs[0])]
            if cur < prev: return False
            prev = cur
        return True
    except Exception:
        return False

def shape_strict(variables, G):
    if not G: return False
    R = G[0].parent()
    gens = {str(g): g for g in R.gens()}
    n = len(variables)
    last = variables[-1]
    # exactly one univariate in the last variable
    uni = [g for g in G if set(map(str,g.variables())) <= {last}]
    if len(uni) != 1: return False
    # for k=n-2..0, we need lm == v_k and degree in v_k equals 1
    for k in range(n-2, -1, -1):
        v = gens[variables[k]]
        cand = [g for g in G if g.lm() == v]
        if len(cand) != 1: return False
        if cand[0].degree(v) != 1: return False
    return True

# ---------- Univariate conversion ----------
def to_true_univariate(poly_mv, var, R, F):
    gens = list(R.gens())
    i = gens.index(var)
    exps = poly_mv.exponents()
    coeffs = poly_mv.coefficients()
    maxd = 0
    for e in exps:
        if sum(e[j] for j in range(len(e)) if j != i) != 0:
            raise TypeError("Not univariate in the requested variable.")
        if e[i] > maxd: maxd = e[i]
    U.<T> = PolynomialRing(F)
    arr = [F(0)]*(maxd+1)
    for e,c in zip(exps, coeffs):
        arr[e[i]] += F(c)
    return U(arr)

# ---------- Triangular back-substitution over a field K ----------
def solve_triangular_over_field(variables, G_lex, K, last_roots,
                                first_only=False, max_solutions=None, log=None):
    """
    Assumes near-triangular LEX basis; starts from roots for last var,
    climbs up using the polynomials with lm == v_k when available.
    Verifies all candidates against G_lex.
    """
    Rk = PolynomialRing(K, variables, order='lex')
    # lift basis to K
    Gk = [Rk(g) for g in G_lex]
    gens = {str(g): g for g in Rk.gens()}
    n = len(variables)
    v_last = gens[variables[-1]]

    # map leading variable -> representative poly
    lead = {}
    for g in Gk:
        lm = g.lm()
        vs = list(lm.variables())
        if len(vs)==1:
            vn = str(vs[0])
            # keep the one with lm == vn if linear; otherwise prefer by degree
            if vn not in lead: lead[vn] = g
            else:
                if g.degree(gens[vn]) < lead[vn].degree(gens[vn]):
                    lead[vn] = g

    sols, count = [], 0
    for r in last_roots:
        assign = {v_last: K(r)}
        # prune quickly
        if any(Rk(g).subs(assign)==K(0) for g in []): pass
        valid = True
        # climb variables n-2..0
        for k in range(n-2, -1, -1):
            vk = gens[variables[k]]
            gk = lead.get(variables[k], None)
            # If no dedicated leading poly, pick any involving vk
            if gk is None:
                cands = [gg for gg in Gk if vk in gg.variables()]
                gk = cands[0] if cands else None
            if gk is None:
                # unconstrained -> set 0
                assign[vk] = K(0); continue
            hk = gk.subs(assign)
            if hk == 0:
                # no constraint at this step; set 0
                assign[vk] = K(0); continue
            # try to solve for vk
            try:
                Uk.<T> = PolynomialRing(K)
                # coerce to true univariate if possible
                fk = to_true_univariate(hk, vk, Rk, K)
            except Exception:
                # attempt linear solve via coefficient extraction
                deg1 = hk.degree(vk)
                if deg1 == 0:
                    # constant≠0 -> infeasible
                    if hk != 0: valid = False; break
                    assign[vk] = K(0); continue
                if deg1 == 1:
                    a = hk.coefficient({vk:1})
                    b = hk.subs({vk:K(0)})
                    if a == 0: valid = False; break
                    assign[vk] = (-b)/a
                    continue
                # high-degree multi-variate: give up on this branch
                valid = False; break
            else:
                roots = [a for (a, m) in fk.roots(multiplicities=True, ring=K)]
                if not roots:
                    valid = False; break
                # pick all alternatives? here we choose the first to keep cheap.
                assign[vk] = roots[0]

        if not valid: continue
        # verification
        if all(gg.subs(assign) == K(0) for gg in Gk):
            sols.append(dict(assign))
            count += 1
            if first_only: break
            if max_solutions is not None and count >= max_solutions: break

    return sols

# ---------- Main driver ----------
def main():
    if len(sys.argv) < 2:
        print("Usage: sage scripts/extract_anemoi_solutions.sage results/<stem>_LEX.txt "
              "[--enum-maxdeg k] [--first-only] [--max-solutions N] [--elim] [--skip-enum]")
        sys.exit(1)

    lex_file = sys.argv[1]
    if not os.path.isfile(lex_file):
        print(f"File not found: {lex_file}"); sys.exit(2)

    # Flags
    enum_maxdeg   = 4
    first_only    = False
    max_solutions = None
    try_elim      = False
    skip_enum     = False

    i = 2
    while i < len(sys.argv):
        a = sys.argv[i]
        if a.startswith("--enum-maxdeg"):
            if "=" in a: enum_maxdeg = int(a.split("=",1)[1]); i += 1
            else:        enum_maxdeg = int(sys.argv[i+1]); i += 2
        elif a == "--first-only":
            first_only = True; i += 1
        elif a.startswith("--max-solutions"):
            if "=" in a: max_solutions = int(a.split("=",1)[1]); i += 1
            else:        max_solutions = int(sys.argv[i+1]); i += 2
        elif a == "--elim":
            try_elim = True; i += 1
        elif a == "--skip-enum":
            skip_enum = True; i += 1
        else:
            print(f"Unrecognized flag: {a}"); sys.exit(3)

    # Paths
    results_dir = "results"; logs_dir = "logs"
    ensure_dir(results_dir); ensure_dir(logs_dir)
    base = stem_of(lex_file)
    report_out = os.path.join(results_dir, base + "_AC_report.txt")
    sols_out   = os.path.join(results_dir, base + "_AC_solutions.txt")
    log_out    = os.path.join(logs_dir,    base + "_AC_extract.log")

    with open(log_out, "w") as log:
        logprint("==============================================================", log)
        logprint(" ALGEBRAIC-CLOSURE SOLUTIONS — LEX GB (Sage)", log)
        logprint("==============================================================", log)
        logprint(f"Input LEX basis: {lex_file}", log)
        logprint(f"enum_maxdeg={enum_maxdeg}, first_only={first_only}, max_solutions={max_solutions}, elim={try_elim}, skip_enum={skip_enum}", log)
        logprint("--------------------------------------------------------------", log)

        # Parse + build ring/ideal
        variables, p, order, polys_s, shape_hdr, zdim_hdr = read_lex_basis_file(lex_file)
        from sage.all import GF, PolynomialRing
        F = GF(p)
        R = PolynomialRing(F, variables, order='lex')
        G_lex = [R(s) for s in polys_s]
        I = R.ideal(G_lex)

        logprint(f"Variables: {variables}", log)
        logprint(f"Field:     GF({p})", log)
        logprint(f"Order:     {order}", log)
        logprint(f"Basis size: {len(G_lex)}", log)
        if shape_hdr is not None:  logprint(f"Header says shape position: {shape_hdr}", log)
        if zdim_hdr  is not None:  logprint(f"Header says zero-dimensional: {zdim_hdr}", log)

        # Shape check
        shp_h = shape_heuristic(variables, G_lex)
        shp_s = shape_strict(variables, G_lex)
        logprint(f"Shape position (heuristic): {shp_h}", log)
        logprint(f"Shape position (strict):    {shp_s}", log)

        # Find the univariate eliminant in the last variable
        v_last_name = variables[-1]
        v_last = R.gens()[-1]
        uni = [g for g in G_lex if set(map(str, g.variables())) <= {v_last_name}]
        eliminant = None
        if uni:
            # choose the one with highest degree in v_last
            eliminant = max(uni, key=lambda g: g.degree(v_last))
        elif try_elim:
            logprint("No univariate found in basis; trying elimination ideal (may be heavy)...", log)
            try:
                El = I.elimination_ideal(R.gens()[:-1])
                cand = [g for g in El.gens() if set(map(str,g.variables())) <= {v_last_name}]
                if cand:
                    eliminant = max(cand, key=lambda g: g.degree(v_last))
            except Exception as e:
                logprint(f"[elim] failed: {e}", log)

        if eliminant is None:
            logprint("ERROR: Could not obtain an eliminant in the last variable. Aborting.", log)
            sys.exit(4)

        # Convert to true univariate f(T) over F
        f_last = to_true_univariate(eliminant, v_last, R, F)
        deg_f  = f_last.degree()
        sf     = f_last.squarefree_part()
        deg_sf = sf.degree()

        # Factorization over F_p (to see required extensions & multiplicities)
        fac = f_last.factor()  # [(h,e), ...]
        # Build factor summary
        breakdown = {}  # degree d -> (sum multiplicity-degree, num factors, total roots over K)
        for (h,e) in fac:
            d = h.degree()
            # total roots with multiplicity from this block over algebraic closure: d*e
            if d not in breakdown:
                breakdown[d] = {"factors":0, "total_mult_roots":0, "distinct_roots":0}
            breakdown[d]["factors"] += 1
            breakdown[d]["total_mult_roots"] += d*e
            breakdown[d]["distinct_roots"]   += d  # one per root (each distinct over \bar{F}_p)

        # Write report (counts over algebraic closure)
        with open(report_out, "w") as rep:
            rep.write("# Algebraic-closure count report\n")
            rep.write(f"# Source: {lex_file}\n")
            rep.write(f"# Field: GF({p})\n")
            rep.write(f"# Variables: {', '.join(variables)}\n")
            rep.write(f"Eliminant variable: {v_last_name}\n")
            rep.write(f"deg(f_last) (total roots, with multiplicity): {deg_f}\n")
            rep.write(f"deg(squarefree(f_last)) (distinct roots):     {deg_sf}\n")
            rep.write("Factorization over GF(p):\n")
            for d in sorted(breakdown):
                info = breakdown[d]
                rep.write(f"  - degree {d}: factors={info['factors']}, "
                          f"distinct_roots={info['distinct_roots']}, "
                          f"mult_roots={info['total_mult_roots']}\n")
        logprint(f"Wrote algebraic-closure count report: {report_out}", log)

        # Optional enumeration in small extensions
        all_solutions = []
        if not skip_enum:
            logprint("Enumerating solutions in minimal extensions for small factor degrees...", log)
            t0 = time.time()
            enumerated = 0
            for (h,e) in fac:
                d = h.degree()
                if d == 0: continue
                if d > enum_maxdeg:
                    logprint(f"  [skip] factor deg {d} > enum_maxdeg {enum_maxdeg}", log)
                    continue
                from sage.all import GF as GFq
                K = GFq(p**d, name='a')
                RK.<T> = PolynomialRing(K)
                hk = RK(h)
                roots = hk.roots(multiplicities=True, ring=K)
                logprint(f"  [deg {d}] roots in GF(p^{d}): {len(roots)} (with multiplicities)", log)

                # back-substitute each root to a full solution
                sols_k = solve_triangular_over_field(
                    variables, G_lex, K, [r for (r,m) in roots],
                    first_only=first_only, max_solutions=max_solutions, log=log
                )
                enumerated += len(sols_k)
                all_solutions.extend([(K, sol) for sol in sols_k])
                if first_only or (max_solutions is not None and enumerated >= max_solutions):
                    break
            t1 = time.time()
            logprint(f"Enumeration done in {t1 - t0:.3f} s — {len(all_solutions)} solution(s) found.", log)

            # Write solutions (pretty)
            if all_solutions:
                with open(sols_out, "w") as out:
                    out.write(f"# Solutions (possibly in extensions) for {lex_file}\n")
                    out.write(f"# Each solution is printed as {variables} over its field.\n")
                    for (K, sol) in all_solutions:
                        Rk = PolynomialRing(K, variables, order='lex')
                        gens = {str(g): g for g in Rk.gens()}
                        vals = []
                        for v in variables:
                            vals.append(str(sol[gens[v]]))
                        out.write(f"GF({K.order()}): " + "{" + ", ".join(f"{v}: {vals[i]}" for i,v in enumerate(variables)) + "}\n")
                logprint(f"Wrote {len(all_solutions)} enumerated solution(s) to: {sols_out}", log)
            else:
                logprint("No solutions were enumerated under current limits (try raising --enum-maxdeg).", log)
        else:
            logprint("Enumeration skipped per --skip-enum.", log)

        logprint("--------------------------------------------------------------", log)
        logprint("Algebraic-closure extraction COMPLETE.", log)

if __name__ == "__main__":
    main()
