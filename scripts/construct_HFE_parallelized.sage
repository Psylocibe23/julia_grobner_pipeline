###############################################################################
# construct_HFE_system.sage  —  Classical HFE instance generator for GF(2)
#
# OUTPUT (pipeline .in, no comments):
#   line 1: x0, x1, ..., x{n-1}
#   line 2: 2
#   next n lines: public equations over GF(2)
#   last n lines: field equations  x_i^2 + x_i
#
# CLASSICAL HFE (char 2):
#   F(X) uses only exponents 2^i (linearized) and 2^i+2^j with i<j (true HFE quad),
#   subject to 2^i ≤ D and 2^i + 2^j ≤ D.
#
# LOGS:
#   • Regular generation log (progress + summary) — includes the secret vector:
#       logs/HFE_n{n}_D{D}_genlog.txt
#   • Secret log (parse-friendly details to verify attacks):
#       logs/HFE_n{n}_D{D}_secret.txt
#     contains Field/modulus, F(X), A_S, b_S, A_T, b_T, Secret.
#
# PERFORMANCE UPGRADES IN THIS VERSION
#   1) Precompute all ∂p_i/∂x_j symbolically ONCE per (S,T) map (big speed-up).
#   2) Probe secrets in PARALLEL across CPU workers until rank ≥ rank_min is found.
#      - Workers default to $SLURM_CPUS_PER_TASK or os.cpu_count().
#      - Optional CLI override: add a [workers] arg at the very end.
###############################################################################

import sys, os, time, random, multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from sage.all import *

# ============================================================================
# Utilities
# ============================================================================

def ensure_dir_for(path):
    """Create parent directory for `path`, if needed (idempotent)."""
    d = os.path.dirname(path)
    if d and not os.path.exists(d):
        os.makedirs(d)

def now_str():
    """Wall-clock time string for user-facing messages."""
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def log_write(L, s):
    """Write a line to an open log file (if any)."""
    if L is not None:
        L.write(s + "\n"); L.flush()

def pick_workers(cli_workers=None):
    """
    Decide worker count in this order:
      1) explicit CLI override
      2) SLURM_CPUS_PER_TASK
      3) os.cpu_count()
      4) fallback = 1
    """
    if cli_workers is not None:
        try:
            w = int(cli_workers)
            return max(1, w)
        except Exception:
            pass
    env = os.environ.get("SLURM_CPUS_PER_TASK")
    if env:
        try:
            return max(1, int(env))
        except Exception:
            pass
    try:
        return max(1, os.cpu_count() or 1)
    except Exception:
        return 1

# ============================================================================
# HFE univariate over K = GF(2^n)
# ============================================================================

def random_hfe_polynomial(K, n, D, prob_quad=0.6, prob_lin=0.6,
                          must_have_quad=True, must_have_lin=True):
    """
    Build F(X) ∈ K[X] with HFE-degree ≤ D using:
      • linearized terms: X^(2^i) when 2^i ≤ D;
      • TRUE quadratic terms: X^(2^i+2^j) with i<j and 2^i+2^j ≤ D.
    (We exclude i=j because X^(2^i+2^i) = X^(2^{i+1}) in char 2 → linearized.)
    """
    R.<X> = PolynomialRing(K)
    F = R(0)

    # admissible 2^i
    lin_is = []
    e = 1; i = 0
    while i < n and e <= D:
        lin_is.append(i)
        i += 1; e <<= 1

    # admissible 2^i+2^j with i<j
    quad_pairs = []
    for i in range(n):
        e_i = 1 << i
        if e_i > D: break
        for j in range(i+1, n):           # STRICT i<j
            e_ij = e_i + (1 << j)
            if e_ij > D: break
            quad_pairs.append((i, j))

    have_quad = False; have_lin = False

    # quadratic terms
    for (i, j) in quad_pairs:
        if random.random() < prob_quad:
            c = K.random_element()
            if c != 0:
                F += c * X**((1 << i) + (1 << j))
                have_quad = True

    # linearized terms
    for i in lin_is:
        if random.random() < prob_lin:
            c = K.random_element()
            if c != 0:
                F += c * X**(1 << i)
                have_lin = True

    # constant term (optional)
    c0 = K.random_element()
    if c0 != 0:
        F += c0

    # enforce presence if admissible
    if must_have_quad and not have_quad and quad_pairs:
        i, j = quad_pairs[0]
        F += K(1) * X**((1 << i) + (1 << j))
        have_quad = True
    if must_have_lin and not have_lin and lin_is:
        i = lin_is[0]
        F += K(1) * X**(1 << i)
        have_lin = True

    return F, have_quad, have_lin

# ============================================================================
# Project K[x] → n coordinates in F2[x]
# ============================================================================

def coords_over_F2(poly_Kx, K, a, n, R_F2):
    """
    Given poly_Kx ∈ K[x0,...,x{n-1}], expand coefficients on the power basis
    {1, a, ..., a^{n-1}} and return [g_0,...,g_{n-1}] with g_t ∈ F2[x].
    """
    coords = [R_F2(0) for _ in range(n)]
    for mon, coeff in poly_Kx.dict().items():
        cvec = K(coeff)._vector_()     # length n over GF(2)
        mon_R = R_F2({mon: 1})
        for t in range(n):
            if cvec[t] != 0:
                coords[t] += mon_R
    return coords

def boolean_reduce(poly, R):
    """
    Boolean reduction modulo <x_i^2 - x_i>:
      for each variable, any exponent e >= 1 collapses to 1 (so x^k = x for k≥1).
    """
    terms = {}
    for mon, coeff in poly.dict().items():
        if coeff == 0:
            continue
        red_mon = tuple(1 if e > 0 else 0 for e in mon)
        terms[red_mon] = terms.get(red_mon, R.base_ring()(0)) + coeff
    return R(terms)

# ============================================================================
# Jacobian helpers: PRECOMPUTE once per (S,T) + fast evaluation
# ============================================================================

def precompute_jacobian_derivatives(polys, R):
    """
    Precompute all partial derivatives d p_i / d x_j symbolically ONCE per map.
    Returns derivs[i][j] ∈ R.
    """
    xs = R.gens()
    n = len(xs)
    return [[polys[i].derivative(xs[j]) for j in range(n)] for i in range(n)]

def jacobian_rank_at_from_derivs(derivs, R, x_star):
    """
    Evaluate J(x_star) from precomputed 'derivs' and return GF(2)-rank.
    """
    F2 = R.base_ring(); xs = R.gens(); n = len(xs)
    subst = {xs[j]: F2(int(x_star[j])) for j in range(n)}
    J = Matrix(F2, n, n)
    for i in range(n):
        for j in range(n):
            J[i, j] = F2(derivs[i][j].subs(subst))
    return J.rank()

# ============================================================================
# Secret log writer
# ============================================================================

def write_secret_log(secret_logfile, n, K, F_univar, A_S, b_S, A_T, b_T, secret_vec):
    """Emit a parse-friendly secret log with F, S, T and the secret."""
    ensure_dir_for(secret_logfile)
    with open(secret_logfile, "w") as S:
        S.write(f"Secret log — {now_str()}\n")
        S.write(f"Field: GF(2^{n})\n")
        S.write(f"Modulus polynomial: {K.modulus()}\n")
        S.write(f"F(X) = {F_univar}\n")
        S.write("A_S =\n")
        for i in range(A_S.nrows()):
            S.write("[" + " ".join(str(int(A_S[i,j])) for j in range(A_S.ncols())) + "]\n")
        S.write("b_S = (" + ", ".join(str(int(b)) for b in b_S) + ")\n")
        S.write("A_T =\n")
        for i in range(A_T.nrows()):
            S.write("[" + " ".join(str(int(A_T[i,j])) for j in range(A_T.ncols())) + "]\n")
        S.write("b_T = (" + ", ".join(str(int(b)) for b in b_T) + ")\n")
        S.write("Secret = [" + ", ".join(str(int(b)) for b in secret_vec) + "]\n")

# ============================================================================
# PARALLEL secret probing infrastructure
#   We avoid passing Sage objects to workers (pickling issues) by relying on
#   Linux 'fork': children inherit globals set by the parent. We install the
#   global derivs/R/F2/n into the child at fork-time and only return small
#   Python tuples.
# ============================================================================

# Globals that workers will see after fork (DO NOT MODIFY in workers).
_GLOBALS = {
    "derivs": None,   # matrix of R-polynomials (∂p_i/∂x_j)
    "R":      None,   # polynomial ring R = GF(2)[x...]
    "n":      None,   # number of variables
}

def _install_worker_state(derivs, R, n):
    """Set shared, read-only state for worker processes."""
    _GLOBALS["derivs"] = derivs
    _GLOBALS["R"]      = R
    _GLOBALS["n"]      = n

def _probe_batch(batch_size, rank_min, seed=None):
    """
    Worker task: try 'batch_size' random secrets; return the first one with
    rank ≥ rank_min if found, otherwise return the best rank seen and its point.
    Result is a plain dict with small Python data (no Sage objects).
    """
    if seed is not None:
        random.seed(seed)

    derivs = _GLOBALS["derivs"]; R = _GLOBALS["R"]; n = _GLOBALS["n"]
    if derivs is None or R is None or n is None:
        # Defensive check; in normal use, parent installed globals before submit.
        return {"found": False, "best_rank": -1, "best_x": None}

    best_rank = -1
    best_x    = None

    # Try 'batch_size' independent random secrets
    F2 = R.base_ring()
    for _ in range(batch_size):
        x_tuple = tuple(random.getrandbits(1) for _ in range(n))   # 0/1 ints
        # Evaluate rank quickly using precomputed derivs
        rankJ = jacobian_rank_at_from_derivs(derivs, R, x_tuple)

        if rankJ > best_rank:
            best_rank, best_x = rankJ, x_tuple

        if rankJ >= rank_min:
            return {"found": True, "x": x_tuple, "rank": rankJ,
                    "best_rank": best_rank, "best_x": best_x}

    return {"found": False, "best_rank": best_rank, "best_x": best_x}

# ============================================================================
# Builder: compute in K[x], then project to F2[x], probe in parallel
# ============================================================================

def build_and_export_instance(n, D, out_infile, seed=None,
                              prob_quad=0.70, prob_lin=0.60,
                              max_maps=20, max_secret_tries=2048,
                              rank_min=None, allow_fallback=False,
                              logfile=None, secret_logfile=None,
                              verbose=False, workers=None, batch_size=256):
    """
    Conditions:
     (i)  F has linearized and TRUE quadratic HFE terms (when admissible).
     (ii) Public system has degree ≥ 2 after Boolean reduction.
     (iii)∃ x* with Jacobian rank ≥ rank_min (default n−1).
          If not found but allow_fallback=True, accept the best-rank seen.

    Parallel probing:
      • For each (S,T), precompute Jacobian derivatives once.
      • Split 'max_secret_tries' into batches of size 'batch_size' and farm
        them out to 'workers' processes. Stop early if a full-rank secret is found.

    Writes:
      - .in file for the pipeline
      - generation log (logfile)
      - secret log (secret_logfile) with F, S, T, secret
    """
    # ---- logging files ----
    L = None
    if logfile:
        ensure_dir_for(logfile)
        L = open(logfile, "w")
        log_write(L, f"HFE generation start — {now_str()}")
        log_write(L, f"Params: n={n}, D={D}, prob_quad={prob_quad}, prob_lin={prob_lin}, "
                     f"max_maps={max_maps}, max_secret_tries={max_secret_tries}, "
                     f"rank_min={'auto' if rank_min is None else rank_min}, "
                     f"allow_fallback={allow_fallback}")

    try:
        if seed is not None:
            set_random_seed(int(seed))
            random.seed(int(seed))
            log_write(L, f"Seed set to {int(seed)}")

        if rank_min is None:
            rank_min = max(0, n - 1)

        # Base structures
        F2 = GF(2)
        K.<a> = GF(2**n)                       # field for HFE
        modulus = K.modulus()
        names = tuple(f"x{i}" for i in range(n))
        R = PolynomialRing(F2, n, names=names) # R = F2[x]
        XK = PolynomialRing(K, n, names=names) # XK = K[x]
        x_R  = R.gens()
        x_K  = XK.gens()

        # Draw F(X) in K[X]
        F_univar, have_quad, have_lin = random_hfe_polynomial(
            K, n, D, prob_quad=prob_quad, prob_lin=prob_lin,
            must_have_quad=True, must_have_lin=True
        )
        log_write(L, f"F has true quadratic? {have_quad}; linearized? {have_lin}")

        # Helper: random invertible over GF(2)
        def rnd_inv(n, F2):
            while True:
                M = random_matrix(F2, n, n)
                if M.is_invertible():
                    return M

        # For fallback bookkeeping
        best = {"rank": -1, "bundle": None}  # (A_S,b_S,A_T,z0_R,x_star,degs)

        # Worker count
        W = pick_workers(workers)
        if L: log_write(L, f"Workers: {W}   (batch_size per task: {batch_size})")

        # Try up to max_maps pairs (S,T)
        for amap in range(1, max_maps+1):
            A_S = rnd_inv(n, F2)
            b_S = vector(F2, [F2.random_element() for _ in range(n)])
            A_T = rnd_inv(n, F2)

            # --- Compute S(x) in K[x] cleanly ---
            s_vec_K = []
            for k in range(n):
                expr = XK(0)
                for j in range(n):
                    if int(A_S[k, j]) != 0:
                        expr += x_K[j]
                if int(b_S[k]) != 0:
                    expr += 1
                s_vec_K.append(expr)

            # Embed into K via the fixed basis: sK(x) = Σ s_k(x) a^k ∈ K[x]
            sK = sum(s_vec_K[k] * (a**k) for k in range(n))

            # Evaluate the HFE univariate: yK(x) = F( sK(x) ) ∈ K[x]
            yK = F_univar(sK)

            # Project K[x] → n polynomials over F2[x]
            y_vec_R = coords_over_F2(yK, K, a, n, R)

            # Apply output linear map (no constant yet): z0 = A_T * y
            z0_R = []
            for i in range(n):
                acc = R(0)
                for j in range(n):
                    if int(A_T[i, j]) != 0:
                        acc += y_vec_R[j]
                z0_R.append(acc)

            # Boolean reduction
            z0_R = [boolean_reduce(p, R) for p in z0_R]
            degs = [p.total_degree() for p in z0_R]
            maxdeg = max(degs) if degs else -1
            log_write(L, f"[map {amap}] public max degree after Boolean: {maxdeg}")

            # Condition (ii): we need degree ≥ 2 publicly
            if maxdeg < 2:
                log_write(L, f"[map {amap}] skipped: maxdeg < 2")
                continue

            # >>> PRECOMPUTE Jacobian derivatives ONCE for this map
            derivs = precompute_jacobian_derivatives(z0_R, R)

            # Install globals for workers (children inherit by fork)
            _install_worker_state(derivs, R, n)

            # --- PARALLEL secret probing until rank ≥ rank_min or tries exhausted ---
            remaining = max_secret_tries
            success   = None         # (x_tuple, rankJ) on success
            best_rank_seen = -1
            best_point     = None

            if W <= 1:
                # Single-process fallback (no executor)
                t = 0
                while remaining > 0:
                    bs = min(batch_size, remaining)
                    # Inline batch work
                    res = _probe_batch(bs, rank_min, seed=None)
                    remaining -= bs; t += bs
                    if res.get("found"):
                        success = (tuple(res["x"]), int(res["rank"]))
                        break
                    # Track best
                    if res["best_rank"] > best_rank_seen:
                        best_rank_seen, best_point = int(res["best_rank"]), tuple(res["best_x"])
                    # Heartbeat every ~1024 probes
                    if (t & 1023) == 0:
                        log_write(L, f"[map {amap}] probe {t}/{max_secret_tries} — best rank so far {best_rank_seen}")
            else:
                # Multi-process probing
                ctx = multiprocessing.get_context("fork")
                with ProcessPoolExecutor(max_workers=W, mp_context=ctx) as ex:
                    futures = set()
                    # Seed offset so different workers get different RNG streams (not critical)
                    seed_base = int(time.time()) ^ random.getrandbits(32)

                    def submit_some(limit):
                        nonlocal remaining, futures
                        while remaining > 0 and len(futures) < W:
                            bs = min(batch_size, remaining)
                            seed = seed_base + remaining
                            fut = ex.submit(_probe_batch, bs, rank_min, seed)
                            futures.add(fut)
                            remaining -= bs

                    submit_some(remaining)

                    while futures and success is None:
                        for fut in as_completed(futures, timeout=None):
                            futures.remove(fut)
                            try:
                                res = fut.result()
                            except Exception as e:
                                # If a worker crashes, continue; re-submit a batch
                                res = {"found": False, "best_rank": -1, "best_x": None}

                            if res.get("found"):
                                success = (tuple(res["x"]), int(res["rank"]))
                                break

                            if res["best_rank"] > best_rank_seen:
                                best_rank_seen, best_point = int(res["best_rank"]), tuple(res["best_x"])

                            # Submit next chunk if any left
                            if remaining > 0:
                                bs = min(batch_size, remaining)
                                seed = seed_base + remaining
                                fut2 = ex.submit(_probe_batch, bs, rank_min, seed)
                                futures.add(fut2)
                                remaining -= bs

                        # Heartbeat: every time we drain a wave, print best so far
                        probed = max_secret_tries - remaining
                        if (probed & 1023) == 0:
                            log_write(L, f"[map {amap}] probe {probed}/{max_secret_tries} — best rank so far {best_rank_seen}")

            # Record best in case we need fallback
            if best_rank_seen > best["rank"]:
                # When success is None, we keep z0_R for fallback; when success found we return, so this is fine.
                best["rank"]   = best_rank_seen
                best["bundle"] = (A_S, b_S, A_T, z0_R, vector(F2, best_point) if best_point else None, degs)

            # If we found a good secret, finish the instance
            if success is not None:
                x_tuple, rankJ = success
                x_star = vector(F2, x_tuple)

                # Choose b_T so P(x*) = 0 (compute at field level)
                s_star  = A_S * x_star + b_S
                sK_star = sum(int(s_star[i]) * (a**i) for i in range(n))
                yK_star = F_univar(sK_star)
                yvec_star = vector(F2, K(yK_star)._vector_())
                b_T = A_T * yvec_star

                # Final public polynomials: z = z0 + b_T
                z_fin = [p + R(int(b_T[i])) for i, p in enumerate(z0_R)]

                # Sanity check
                subst_eval = {x_R[i]: F2(int(x_star[i])) for i in range(n)}
                if not all(p.subs(subst_eval) == 0 for p in z_fin):
                    # Extremely unlikely; if it happens, just continue to next map
                    log_write(L, f"[map {amap}] WARNING: P(x*) != 0 after setting b_T. Continuing search.")
                    continue

                # SUCCESS — write .in and logs
                ensure_dir_for(out_infile)
                with open(out_infile, "w") as f:
                    f.write(", ".join(str(v) for v in x_R) + "\n")
                    f.write("2\n")
                    for p in z_fin: f.write(str(p) + "\n")
                    for v in x_R:   f.write(f"{v}^2 + {v}\n")

                log_write(L, f"[map {amap}] success with rank {rankJ}")
                log_write(L, f"OUTPUT: {out_infile}")
                log_write(L, f"Field modulus: {modulus}")
                log_write(L, f"Public degs: {degs}")
                log_write(L, f"Secret: {list(x_star)}")

                if secret_logfile is None:
                    secret_logfile = os.path.join("logs", f"HFE_n{n}_D{D}_secret.txt")
                write_secret_log(secret_logfile, n, K, F_univar, A_S, b_S, A_T, b_T, x_star)

                return {
                    "ok_zero": True, "rankJ": rankJ, "rank_min": rank_min,
                    "A_S": A_S, "b_S": b_S, "A_T": A_T, "b_T": b_T,
                    "F": F_univar, "modulus": modulus, "secret": x_star,
                    "public_polys": z_fin, "R": R,
                    "have_quad_in_F": have_quad, "have_lin_in_F": have_lin,
                    "public_degrees": degs, "map_attempts": amap,
                    "secret_logfile": secret_logfile
                }

        # If we reached here, no success. Try fallback?
        if allow_fallback and best["rank"] >= 0 and best["bundle"] is not None:
            A_S, b_S, A_T, z0_R, x_star, degs = best["bundle"]
            if x_star is None:
                raise RuntimeError("Fallback requested but no candidate point was recorded.")
            s_star  = A_S * x_star + b_S
            sK_star = sum(int(s_star[i]) * (a**i) for i in range(n))
            yK_star = F_univar(sK_star)
            yvec_star = vector(F2, K(yK_star)._vector_())
            b_T = A_T * yvec_star
            z_fin = [p + R(int(b_T[i])) for i, p in enumerate(z0_R)]

            subst_eval = {R.gens()[i]: F2(int(x_star[i])) for i in range(n)}
            ok_zero = all(p.subs(subst_eval) == 0 for p in z_fin)

            ensure_dir_for(out_infile)
            with open(out_infile, "w") as f:
                f.write(", ".join(str(v) for v in R.gens()) + "\n")
                f.write("2\n")
                for p in z_fin: f.write(str(p) + "\n")
                for v in R.gens(): f.write(f"{v}^2 + {v}\n")

            log_write(L, f"FALLBACK OUTPUT: {out_infile} (rank={best['rank']} < {rank_min})")
            log_write(L, f"Field modulus: {modulus}")
            log_write(L, f"Public degs: {degs}")
            log_write(L, f"Secret: {list(x_star)}")

            if secret_logfile is None:
                secret_logfile = os.path.join("logs", f"HFE_n{n}_D{D}_secret.txt")
            write_secret_log(secret_logfile, n, K, F_univar, A_S, b_S, A_T, b_T, x_star)

            return {
                "ok_zero": ok_zero, "rankJ": best["rank"], "rank_min": rank_min,
                "A_S": A_S, "b_S": b_S, "A_T": A_T, "b_T": b_T,
                "F": F_univar, "modulus": modulus, "secret": x_star,
                "public_polys": z_fin, "R": R,
                "have_quad_in_F": have_quad, "have_lin_in_F": have_lin,
                "public_degrees": degs, "map_attempts": max_maps,
                "secret_logfile": secret_logfile
            }

        # Hard failure (no fallback or no candidate at all)
        raise RuntimeError(f"Failed: best Jacobian rank seen: {best['rank']}.")

    finally:
        if L is not None:
            L.close()

# ============================================================================
# MAIN (CLI)
# ============================================================================

def main():
    if len(sys.argv) < 3:
        print("Usage: sage scripts/construct_HFE_system.sage <n> <D> [outfile.in] [seed] "
              "[prob_quad] [prob_lin] [rank_min] [allow_fallback] [max_maps] [max_secret_tries] [workers]")
        sys.exit(1)

    n = int(sys.argv[1]); D = int(sys.argv[2])

    # Output path (treat '', '-', 'default' as auto)
    auto_out = os.path.join("data", "hfe_instances", f"HFE_n{n}_D{D}.in")
    if len(sys.argv) >= 4:
        a3 = sys.argv[3].strip()
        if a3 in ("", "-", "default"):
            out_infile = auto_out; arg_offset = 4
        elif a3.isdigit():
            out_infile = auto_out; arg_offset = 3  # a3 is seed
        else:
            out_infile = a3; arg_offset = 4
    else:
        out_infile = auto_out; arg_offset = 3

    seed             = int(sys.argv[arg_offset])     if len(sys.argv) >= arg_offset+1 else None
    prob_quad        = float(sys.argv[arg_offset+1]) if len(sys.argv) >= arg_offset+2 else 0.70
    prob_lin         = float(sys.argv[arg_offset+2]) if len(sys.argv) >= arg_offset+3 else 0.60
    rank_min_arg     = sys.argv[arg_offset+3]        if len(sys.argv) >= arg_offset+4 else None
    rank_min         = int(rank_min_arg) if rank_min_arg is not None else None
    allow_fallback   = (sys.argv[arg_offset+4].lower() in ("1","true","yes")) if len(sys.argv) >= arg_offset+5 else False
    max_maps         = int(sys.argv[arg_offset+5])   if len(sys.argv) >= arg_offset+6 else 20
    max_secret_tries = int(sys.argv[arg_offset+6])   if len(sys.argv) >= arg_offset+7 else 2048
    workers_arg      = sys.argv[arg_offset+7]        if len(sys.argv) >= arg_offset+8 else None

    workers = pick_workers(workers_arg)

    logname    = os.path.join("logs", f"HFE_n{n}_D{D}_genlog.txt")
    secretlog  = os.path.join("logs", f"HFE_n{n}_D{D}_secret.txt")

    print(f"[{now_str()}] Generating classical HFE instance (n={n}, D={D}) ...")
    print(f"[{now_str()}] Output .in: {out_infile}")
    print(f"[{now_str()}] Workers: {workers}")

    info = build_and_export_instance(
        n, D, out_infile, seed=seed,
        prob_quad=prob_quad, prob_lin=prob_lin,
        max_maps=max_maps, max_secret_tries=max_secret_tries,
        rank_min=rank_min, allow_fallback=allow_fallback,
        logfile=logname, secret_logfile=secretlog, verbose=False,
        workers=workers, batch_size=256
    )
    print(f"[{now_str()}] Log written to:        {logname}")
    print(f"[{now_str()}] Secret log written to: {info.get('secret_logfile', secretlog)}")
    print(f"[{now_str()}] P(x*)=0: {info['ok_zero']}, rank(J_P(x*))={info['rankJ']}")

if __name__ == "__main__":
    main()
