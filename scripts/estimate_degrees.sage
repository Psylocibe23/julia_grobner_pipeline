#!/usr/bin/env sage
# -*- coding: utf-8 -*-
###############################################################################
# estimate_degrees.sage  —  fast, pre-GB degree/toughness estimates for MQ over GF(2)
#
# Usage:
#   sage scripts/estimate_degrees.sage data/HFE_n80_D96.in [--trials 64] [--skip-m3] > logs/HFE_n6_D6_estdeg.txt
###############################################################################

import sys, re, time
from random import randrange
from sage.all import GF, PolynomialRing, matrix, random_vector

# ================= Utilities =================
def parse_in_file(infile):
    with open(infile, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 2:
        raise ValueError("Malformed .in (need ≥ 2 lines).")

    var_names = [v.strip() for v in lines[0].split(",")]
    n = len(var_names)
    if n == 0:
        raise ValueError("No variables on line 1.")

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

# ================ Macaulay bound for d_solv (DRL setting) =================
def macaulay_bound_top(n, degrees):
    if len(degrees) < n+1:
        return None, f"Not enough equations for the bound (need ≥ n+1, have {len(degrees)})."
    dsorted = sorted(degrees, reverse=True)
    top = dsorted[:n+1]
    bound = sum(top) - n
    return bound, f"Top (n+1) degrees = {top} ⇒ bound = {bound}"

# ================= Build degree-2 / degree-3 monomial bases (with rep) =================
def deg2_monomials(R):
    x = R.gens()
    n = len(x)
    monos = []
    for i in range(n):
        xi = x[i]
        for j in range(i, n):
            monos.append(xi * x[j])
    return monos

def deg3_monomials(R):
    x = R.gens()
    n = len(x)
    monos = []
    for i in range(n):
        xi = x[i]
        for j in range(i, n):
            xij = xi * x[j]
            for k in range(j, n):
                monos.append(xij * x[k])
    return monos

# ================= Boolean Macaulay M2: rows=m_public, cols=deg-2 monomials =================
def boolean_macaulay_M2(public_polys, R):
    F2 = GF(2)
    deg2 = deg2_monomials(R)
    col_index = {deg2[c]: c for c in range(len(deg2))}
    ncols = len(deg2)
    nrows = len(public_polys)

    entries = {}
    for i, f in enumerate(public_polys):
        for m, a in zip(f.monomials(), f.coefficients()):
            if m.total_degree() == 2 and (int(a) & 1):
                j = col_index.get(m, None)
                if j is not None:
                    entries[(i, j)] = 1

    M = matrix(F2, nrows, ncols, entries, sparse=True)
    return M, deg2

# ================= Boolean Macaulay M3: rows=m_public*n (x_j * f_i), cols=deg-3 =================
def boolean_macaulay_M3(public_polys, R):
    F2 = GF(2)
    x = R.gens()
    n = len(x)

    deg3 = deg3_monomials(R)
    col_index = {deg3[c]: c for c in range(len(deg3))}
    ncols = len(deg3)
    nrows = len(public_polys) * n

    entries = {}
    row = 0
    for f in public_polys:
        for xj in x:
            g = xj * f
            for m, a in zip(g.monomials(), g.coefficients()):
                if m.total_degree() == 3 and (int(a) & 1):
                    c = col_index.get(m, None)
                    if c is not None:
                        entries[(row, c)] = 1
            row += 1

    M = matrix(F2, nrows, ncols, entries, sparse=True)
    return M, deg3

# ================= First-fall heuristic at degree 3 (randomized, very fast) =================
def first_fall_degree3_heuristic(public_polys, R, trials=64, sample_per_trial=None, rng=None):
    if rng is None:
        from random import Random
        rng = Random()

    x = R.gens()
    n = len(x)
    m = len(public_polys)

    all_pairs = [(i, j) for i in range(m) for j in range(n)]
    total_pairs = len(all_pairs)
    if sample_per_trial is None:
        sample_per_trial = min(50, total_pairs)

    for t in range(trials):
        if sample_per_trial >= total_pairs:
            chosen = all_pairs
        else:
            chosen_idx = set()
            while len(chosen_idx) < sample_per_trial:
                chosen_idx.add(rng.randrange(total_pairs))
            chosen = [all_pairs[k] for k in chosen_idx]

        H = R(0)
        for (i, j) in chosen:
            H += x[j] * public_polys[i]

        if H == 0:
            return True, t+1

        deg3_part = R(0)
        for mmono, coeff in zip(H.monomials(), H.coefficients()):
            if mmono.total_degree() == 3 and (int(coeff) & 1):
                deg3_part += coeff * mmono

        if deg3_part == 0 and H != 0:
            return True, t+1

    return False, trials

def approx_nnz(M):
    """
    Try to count nonzeros; if not available, estimate via density.
    """
    # Try a reliable path first
    try:
        return len(list(M.nonzero_positions()))
    except Exception:
        pass
    # Fall back to density-based estimate
    try:
        return int(round(M.nrows() * M.ncols() * float(M.density())))
    except Exception:
        return None

# ================= Main =================
def main():
    if len(sys.argv) < 2:
        print("Usage: sage scripts/estimate_degrees.sage <input.in> [--trials N] [--skip-m3]")
        sys.exit(1)

    infile = sys.argv[1]
    trials = 64
    do_m3 = True

    args = sys.argv[2:]
    i = 0
    while i < len(args):
        if args[i] == "--trials" and i+1 < len(args):
            trials = int(args[i+1]); i += 2
        elif args[i] == "--skip-m3":
            do_m3 = False; i += 1
        else:
            print(f"Unrecognized argument: {args[i]}")
            sys.exit(2)

    var_names, F, R, public_polys, field_polys = parse_in_file(infile)
    n = len(var_names)
    m_pub = len(public_polys)

    print("==============================================================")
    print(" FAST DEGREE ESTIMATES (pre-GB) — Boolean MQ over GF(2)")
    print("==============================================================")
    print(f"System: {infile}")
    print(f"Variables: n = {n}   ({var_names})")
    print(f"Public eqs: m = {m_pub}")
    print(f"Field eqs: {'present' if field_polys is not None else 'absent / unrecognized'}")
    print("--------------------------------------------------------------")

    degs_pub = [f.total_degree() for f in public_polys]
    print(f"Public degrees: {degs_pub}  (max={max(degs_pub) if degs_pub else 'NA'})")

    bound, detail = macaulay_bound_top(n, degs_pub)
    if bound is None:
        print("Macaulay bound: N/A  (need at least n+1 public equations).")
    else:
        print(f"Macaulay bound on d_solv: ≤ {bound}")
        print(f"  Details: {detail}")

    # Degree 2
    t0 = time.time()
    M2, deg2cols = boolean_macaulay_M2(public_polys, R)
    r2 = M2.rank()
    t1 = time.time()
    nnz2 = approx_nnz(M2)
    nnz2_str = f"{nnz2}" if nnz2 is not None else "N/A"
    print("--------------------------------------------------------------")
    print("Boolean Macaulay @ degree 2:")
    print(f"  size: rows={M2.nrows()}  cols={M2.ncols()}   nnz≈{nnz2_str}")
    print(f"  rank over GF(2): {r2}")
    print(f"  time: {t1 - t0:.3f} s")
    if r2 == M2.ncols():
        print("  -> Full column rank at degree 2 -> predict d_solv ≤ 3")

    # Degree 3
    if do_m3:
        t2 = time.time()
        M3, deg3cols = boolean_macaulay_M3(public_polys, R)
        r3 = M3.rank()
        t3 = time.time()
        nnz3 = approx_nnz(M3)
        nnz3_str = f"{nnz3}" if nnz3 is not None else "N/A"
        print("--------------------------------------------------------------")
        print("Boolean Macaulay @ degree 3:")
        print(f"  size: rows={M3.nrows()}  cols={M3.ncols()}   nnz≈{nnz3_str}")
        print(f"  rank over GF(2): {r3}")
        print(f"  time: {t3 - t2:.3f} s")
        if r2 != M2.ncols() and r3 == M3.ncols():
            print("  -> Full column rank at degree 3 -> predict d_solv ≤ 4")
    else:
        print("--------------------------------------------------------------")
        print("Boolean Macaulay @ degree 3: skipped (use --skip-m3 to skip / default runs it).")

    # First-fall heuristic at degree 3
    ff_found, used = first_fall_degree3_heuristic(public_polys, R, trials=trials)
    print("--------------------------------------------------------------")
    if ff_found:
        print(f"First-fall heuristic @ degree 3: FOUND in ≤ {used} trial(s) ⇒ predict d_ff ≤ 3")
    else:
        print(f"First-fall heuristic @ degree 3: none found in {trials} trials")

    print("--------------------------------------------------------------")
    print("Done.")

if __name__ == "__main__":
    main()
