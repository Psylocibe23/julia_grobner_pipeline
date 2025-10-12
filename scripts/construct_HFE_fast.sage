###############################################################################
# construct_HFE_fast.sage — Fast, bounded, fail-safe HFE instance generator over GF(2)
#
# Behavior:
#   • Tries to generate an instance with target degree D and Jacobian rank n-1.
#   • Hard internal timer and bounded attempts; always emits best-so-far.
#   • Outputs the .in file, a public genlog, and a secret log (F,S,T,secret).
###############################################################################

import sys, os, time, random as pyrand
from sage.all import *
import operator

# --------------------------- helpers ---------------------------
def ensure_dir_for(path):
    d = os.path.dirname(path)
    if d and not os.path.exists(d):
        os.makedirs(d)

def now_str():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def is_unset(s):
    return (s is None) or (str(s).strip().lower() in {"", "-", "none", "default"})

def split_to_F2_vec(elem, n):
    """
    Coordinates of elem ∈ GF(2^n) on the construction power basis as GF(2)^n.
    """
    F2 = elem.parent().base_ring()
    try:
        v = elem._vector_()
        v_list = list(v)
        if len(v_list) < n:
            v_list += [F2(0)]*(n-len(v_list))
        else:
            v_list = v_list[:n]
        return vector(F2, v_list)
    except Exception:
        pass
    pol = elem.polynomial()  # degree < n over GF(2)
    coeffs = [F2(c) for c in pol.list()]
    if len(coeffs) < n:
        coeffs += [F2(0)]*(n-len(coeffs))
    else:
        coeffs = coeffs[:n]
    return vector(F2, coeffs)

# --------------------------- logs ---------------------------
def open_genlog(n, D):
    log_path = os.path.join("logs", f"HFE_n{n}_D{D}_genlog.txt")
    ensure_dir_for(log_path)
    L = open(log_path, "w")
    return L, log_path

def write_secret_log(secret_logfile, n, K, F_univar, A_S, b_S, A_T, b_T, secret_vec):
    ensure_dir_for(secret_logfile)
    with open(secret_logfile, "w") as S:
        S.write(f"Secret log — {now_str()}\n")
        S.write(f"Field: GF(2^{n})\n")
        S.write(f"Modulus polynomial: {K.modulus()}\n")
        S.write(f"F(X) = {F_univar}\n")
        S.write("A_S =\n")
        for i in range(A_S.nrows()):
            S.write("[" + " ".join(str(int(A_S[i, j])) for j in range(A_S.ncols())) + "]\n")
        S.write("b_S = (" + ", ".join(str(int(b)) for b in b_S) + ")\n")
        S.write("A_T =\n")
        for i in range(A_T.nrows()):
            S.write("[" + " ".join(str(int(A_T[i, j])) for j in range(A_T.ncols())) + "]\n")
        S.write("b_T = (" + ", ".join(str(int(b)) for b in b_T) + ")\n")
        S.write("Secret = [" + ", ".join(str(int(b)) for b in secret_vec) + "]\n")

# ----- Boolean "degree ≤ 2" Jacobian via bitmasks -----
def _build_jacobian_bitmasks(polys, R):
    n = len(R.gens())
    out = []
    for p in polys:
        lin_mask = 0
        quad_rows = [0]*n
        for mon, coeff in p.dict().items():
            if coeff == 0:
                continue
            idxs = [t for t, e in enumerate(mon) if e]  # Boolean reduced → e∈{0,1}
            if not idxs:
                continue
            if len(idxs) == 1:
                j = idxs[0]
                lin_mask |= (1 << j)
            elif len(idxs) == 2:
                u, v = idxs
                if u != v:
                    quad_rows[u] |= (1 << v)
                    quad_rows[v] |= (1 << u)
        out.append((lin_mask, quad_rows))
    return out

def popcount_int(v):
    try:
        return int(v).bit_count()
    except AttributeError:
        v = int(v); c = 0
        while v:
            v &= v-1; c += 1
        return c

def _gf2_rank_from_bitrows(bitrows, n):
    rows = bitrows[:]
    rank = 0
    for col in range(n):
        mask = 1 << col
        piv = None
        for r in range(rank, len(rows)):
            if rows[r] & mask:
                piv = r; break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        pivrow = rows[rank]
        # XOR-eliminate this column from all other rows 
        for r in range(len(rows)):
            if r != rank and (rows[r] & mask):
                rows[r] = operator.xor(int(rows[r]), int(pivrow))
        rank += 1
        if rank == n:
            break
    return rank

def jacobian_rank_bitset(precomp, x_tuple):
    n = len(x_tuple)
    xmask = 0
    for j, b in enumerate(x_tuple):
        if b: xmask |= (1 << j)
    bitrows = []
    for (lin_mask, qrows) in precomp:
        rowmask = 0
        lin_mask_int = int(lin_mask)
        for j in range(n):
            val = (lin_mask_int >> j) & 1
            val = operator.xor(int(val), (popcount_int(int(qrows[j]) & xmask) & 1))
            if val:
                rowmask |= (1 << j)
        bitrows.append(rowmask)
    return _gf2_rank_from_bitrows(bitrows, n)

# ----- Booleanize: collapse x^k to x (k>0) in each coordinate -----
def boolean_reduce(poly, R):
    terms = {}
    for mon, coeff in poly.dict().items():
        if coeff == 0: continue
        red = tuple(1 if e > 0 else 0 for e in mon)
        terms[red] = (terms.get(red, R.base_ring()(0)) + coeff)
    return R(terms)

def booleanize_vec(vec_R, R):
    return [boolean_reduce(p, R) for p in vec_R]

# ----- Univariate HFE builder over K=GF(2^n) (≤ D) -----
def random_hfe_univariate(K, n, D, prob_quad=0.7, prob_lin=0.6, must_quad=True, must_lin=True):
    R.<X> = PolynomialRing(K)
    F = R(0)
    # linearized exponents 2^i
    lin_is = []
    e = 1; i = 0
    while i < n and e <= D:
        lin_is.append(i)
        i += 1; e <<= 1
    # quadratic exponents 2^i+2^j
    quad_pairs = []
    for i in range(n):
        e_i = 1 << i
        if e_i > D: break
        for j in range(i+1, n):
            e_ij = e_i + (1 << j)
            if e_ij > D: break
            quad_pairs.append((i, j))
    have_q = False; have_l = False
    used_max_exp = 0

    # quadratic terms
    for (i, j) in quad_pairs:
        if pyrand.random() < prob_quad:
            c = K.random_element()
            if c != 0:
                exp = (1 << i) + (1 << j)
                if exp > used_max_exp: used_max_exp = exp
                F += c * X^exp
                have_q = True
    # linearized terms
    for i in lin_is:
        if pyrand.random() < prob_lin:
            c = K.random_element()
            if c != 0:
                exp = (1 << i)
                if exp > used_max_exp: used_max_exp = exp
                F += c * X^exp
                have_l = True
    # constant term
    c0 = K.random_element()
    if c0 != 0:
        F += c0

    # ensure presence if requested and admissible
    if must_quad and (not have_q) and quad_pairs:
        i, j = quad_pairs[0]
        exp = (1 << i) + (1 << j)
        F += K(1) * X^exp
        used_max_exp = max(used_max_exp, exp)
        have_q = True
    if must_lin and (not have_l) and lin_is:
        i = lin_is[0]
        exp = (1 << i)
        F += K(1) * X^exp
        used_max_exp = max(used_max_exp, exp)
        have_l = True

    return F, used_max_exp, have_q, have_l

# ----- Fast public map evaluation: K → F2^n split via power basis -----
def fast_public_from_F_and_S(K, a, R, A_S, b_S, F_univar, n, D):
    # s(x) = c0 + Σ_v c_v x_v  in K, with c_• determined by A_S,b_S and basis {a^k}
    c = [K(0)]*(n+1)    # c[0]=const; c[v+1] for x_v
    for k in range(n):
        ak = a^k
        if int(b_S[k]) != 0:
            c[0] += ak
        for v in range(n):
            if int(A_S[k, v]) != 0:
                c[v+1] += ak

    # largest 2^t ≤ D
    max_t = 0; e = 1
    while (e << 1) <= D:
        e <<= 1; max_t += 1
    # precompute c_i^{2^t}
    c_pow = [[K(0)]*(n+1) for _ in range(max_t+1)]
    for i in range(n+1):
        val = c[i]
        c_pow[0][i] = val
        for t in range(1, max_t+1):
            val = val^2
            c_pow[t][i] = val

    # accumulate y ∈ (F2[x])^n via coefficient splitting
    y = [R(0) for _ in range(n)]
    zero_mon = tuple(0 for _ in range(n))

    def add_Kcoef_mon(Kcoef, mon_tuple):
        if Kcoef == 0:
            return
        vec = split_to_F2_vec(Kcoef, n)
        mono = R({mon_tuple: 1})
        for t in range(n):
            if int(vec[t]) != 0:
                y[t] += mono

    for exp, Kc in F_univar.dict().items():
        eX = int(exp)
        if eX == 0:
            add_Kcoef_mon(Kc, zero_mon); continue
        # power of two: 2^i
        if (eX & (eX - 1)) == 0:
            i = eX.bit_length() - 1
            c0i = c_pow[i][0]
            if c0i != 0:
                add_Kcoef_mon(Kc * c0i, zero_mon)
            for v in range(n):
                Cv = c_pow[i][v+1]
                if Cv != 0:
                    mon = tuple(1 if j == v else 0 for j in range(n))
                    add_Kcoef_mon(Kc * Cv, mon)
            continue
        # binary popcount == 2  → 2^i+2^j
        bits = []
        tmp = eX; pos = 0
        while tmp:
            if tmp & 1: bits.append(pos)
            tmp >>= 1; pos += 1
        if len(bits) != 2:
            continue
        i, j = bits[0], bits[1]
        c0i, c0j = c_pow[i][0], c_pow[j][0]
        if c0i != 0 and c0j != 0:
            add_Kcoef_mon(Kc * c0i * c0j, zero_mon)
        for v in range(n):
            Ci_v = c_pow[i][v+1]; Cj_v = c_pow[j][v+1]
            if c0i != 0 and Cj_v != 0:
                mon = tuple(1 if t == v else 0 for t in range(n))
                add_Kcoef_mon(Kc * c0i * Cj_v, mon)
            if c0j != 0 and Ci_v != 0:
                mon = tuple(1 if t == v else 0 for t in range(n))
                add_Kcoef_mon(Kc * c0j * Ci_v, mon)
        for v in range(n):
            Ci_v = c_pow[i][v+1]
            if Ci_v == 0: continue
            for u in range(n):
                Cj_u = c_pow[j][u+1]
                if Cj_u == 0: continue
                if v == u:
                    mon = tuple(1 if t == v else 0 for t in range(n))
                else:
                    mon = tuple(1 if (t == v or t == u) else 0 for t in range(n))
                add_Kcoef_mon(Kc * Ci_v * Cj_u, mon)
    return y

def rnd_invertible(n, F2):
    while True:
        M = random_matrix(F2, n, n)
        if M.is_invertible():
            return M

# --------------------------- core builder ---------------------------
def build_and_export_fast(n, D, out_infile,
                          seed=None,
                          time_budget_sec=5400,
                          max_maps=40,
                          secret_tries_per_map=2048,
                          rank_min=None,
                          force=False,
                          prob_quad=0.70, prob_lin=0.60):
    """
    Returns dict with summary; ALWAYS writes an instance (success or fallback).
    Also writes genlog and secret log in the classic format.
    """
    # logs
    L, genlog_path = open_genlog(n, D)
    def lw(s): L.write(s + "\n"); L.flush()

    # header identical to the parallel version
    lw(f"HFE generation start — {now_str()}")
    lw("Params: n=%d, D=%d, prob_quad=%.15f, prob_lin=%.15f, max_maps=%d, max_secret_tries=%d, rank_min=%s, allow_fallback=%s"
       % (n, D, prob_quad, prob_lin, max_maps, secret_tries_per_map,
          "auto" if rank_min is None else str(rank_min), "True"))

    # file existence guard
    if (not force) and os.path.exists(out_infile):
        lw(f"Found existing instance: {out_infile} — skipping.")
        L.close()
        return {"skipped": True, "outfile": out_infile, "genlog": genlog_path}

    if seed is not None:
        set_random_seed(int(seed)); pyrand.seed(int(seed))

    start = time.time()
    if rank_min is None: rank_min = max(0, n - 1)

    F2 = GF(2)
    K.<a> = GF(2^n)
    names = tuple(f"x{i}" for i in range(n))
    R = PolynomialRing(F2, n, names=names)

    # Univariate F (log quadratic/linearized presence)
    F_univar, used_exp, have_q, have_l = random_hfe_univariate(
        K, n, D, prob_quad=prob_quad, prob_lin=prob_lin, must_quad=True, must_lin=True
    )
    lw(f"F has true quadratic? {bool(have_q)}; linearized? {bool(have_l)}")

    # Workers line (fixed 1; batch size shown as tries per map)
    lw(f"Workers: 1   (batch_size per task: {secret_tries_per_map})")

    D_used = max(used_exp, 0)

    best = {"rank": -1, "bundle": None}  # (A_S, b_S, A_T, z0_R, x_best_tuple_or_None, degs)
    attempts = 0
    success = None

    while attempts < max_maps and (time.time() - start) < time_budget_sec:
        attempts += 1

        A_S = rnd_invertible(n, F2)
        b_S = vector(F2, [F2.random_element() for _ in range(n)])

        lw(f"[map {attempts}] building public polynomials via Frobenius fast path ...")
        t0 = time.time()
        y_vec_R = fast_public_from_F_and_S(K, a, R, A_S, b_S, F_univar, n, D)
        t_build = time.time() - t0

        # choose A_T with nontrivial Boolean degree
        chosen_A_T = None
        z0_R = None
        degs = None
        for _ in range(8):
            A_T_try = rnd_invertible(n, F2)
            cand = []
            for i in range(n):
                acc = R(0)
                for j in range(n):
                    if int(A_T_try[i, j]) != 0:
                        acc += y_vec_R[j]
                cand.append(acc)
            cand_bool = [boolean_reduce(p, R) for p in cand]
            degs_try = [p.total_degree() for p in cand_bool]
            if max(degs_try) >= 2:
                chosen_A_T = A_T_try
                z0_R = cand_bool
                degs = degs_try
                break

        lw(f"[map {attempts}] public system built in {t_build:.2f}s")
        if z0_R is None:
            lw(f"[map {attempts}] all A_T trials produced maxdeg<2; resampling.")
            continue

        lw(f"[map {attempts}] public max degree after Boolean: {max(degs)}")

        # Jacobian precompute & random probing
        bitJ = _build_jacobian_bitmasks(z0_R, R)
        best_rank = -1
        best_x_tuple = None

        for t in range(secret_tries_per_map):
            x_tuple = tuple(pyrand.getrandbits(1) for _ in range(n))
            rankJ = jacobian_rank_bitset(bitJ, x_tuple)
            if rankJ > best_rank:
                best_rank, best_x_tuple = rankJ, x_tuple
                if best_rank >= rank_min:
                    success = (A_S, b_S, chosen_A_T, z0_R, vector(F2, best_x_tuple), degs, best_rank)
                    break
            if (t & 255) == 0 and (time.time() - start) >= time_budget_sec:
                break

        if best_rank > best["rank"]:
            best["rank"] = best_rank
            best["bundle"] = (A_S, b_S, chosen_A_T, z0_R,
                              tuple(best_x_tuple) if best_x_tuple is not None else None,
                              degs)

        if success is not None:
            break

    # choose output (success or fallback), compute b_T, write files + logs + secret log
    secretlog_path = os.path.join("logs", f"HFE_n{n}_D{D}_secret.txt")
    modulus = K.modulus()

    if success is None:
        if best["bundle"] is None:
            # very unlikely: force one skeleton
            A_S = rnd_invertible(n, F2)
            b_S = vector(F2, [F2.random_element() for _ in range(n)])
            y_vec_R = fast_public_from_F_and_S(K, a, R, A_S, b_S, F_univar, n, D)
            A_T = rnd_invertible(n, F2)
            z0_raw = []
            for i in range(n):
                acc = R(0)
                for j in range(n):
                    if int(A_T[i, j]) != 0:
                        acc += y_vec_R[j]
                z0_raw.append(acc)
            z0_R = [boolean_reduce(p, R) for p in z0_raw]
            degs = [p.total_degree() for p in z0_R]
            x_star_tuple = None
        else:
            A_S, b_S, A_T, z0_R, x_star_tuple, degs = best["bundle"]

        x_star_vec = vector(F2, [0]*n) if x_star_tuple is None else vector(F2, x_star_tuple)
        s_star = A_S * x_star_vec + b_S
        sK_star = sum(int(s_star[i]) * (a^i) for i in range(n))
        yK_star = F_univar(sK_star)
        yvec_star = split_to_F2_vec(yK_star, n)
        b_T = A_T * yvec_star
        z_fin = [z0_R[i] + R(int(b_T[i])) for i in range(n)]

        ensure_dir_for(out_infile)
        with open(out_infile, "w") as f:
            f.write(", ".join(str(v) for v in R.gens()) + "\n")
            f.write("2\n")
            for p in z_fin: f.write(str(p) + "\n")
            for v in R.gens(): f.write(f"{v}^2 + {v}\n")

        lw(f"OUTPUT: {out_infile}")
        lw(f"Field modulus: {modulus}")
        lw(f"Public degs: {degs}")
        lw(f"Secret: [{', '.join(str(int(b)) for b in x_star_vec)}]")

        write_secret_log(secretlog_path, n, K, F_univar, A_S, b_S, A_T, b_T, x_star_vec)
        L.close()

        return {
            "ok_zero": True, "rankJ": best["rank"], "rank_min": rank_min,
            "D_target": D, "D_used": D_used, "public_maxdeg": max(degs),
            "outfile": out_infile, "attempts": attempts,
            "genlog": genlog_path, "secretlog": secretlog_path
        }

    else:
        A_S, b_S, A_T, z0_R, x_star_vec, degs, rankJ = success
        s_star = A_S * x_star_vec + b_S
        sK_star = sum(int(s_star[i]) * (a^i) for i in range(n))
        yK_star = F_univar(sK_star)
        yvec_star = split_to_F2_vec(yK_star, n)
        b_T = A_T * yvec_star
        z_fin = [z0_R[i] + R(int(b_T[i])) for i in range(n)]

        ensure_dir_for(out_infile)
        with open(out_infile, "w") as f:
            f.write(", ".join(str(v) for v in R.gens()) + "\n")
            f.write("2\n")
            for p in z_fin: f.write(str(p) + "\n")
            for v in R.gens(): f.write(f"{v}^2 + {v}\n")

        lw(f"[map {attempts}] success with rank {rankJ}")
        lw(f"OUTPUT: {out_infile}")
        lw(f"Field modulus: {modulus}")
        lw(f"Public degs: {degs}")
        lw(f"Secret: [{', '.join(str(int(b)) for b in x_star_vec)}]")

        write_secret_log(secretlog_path, n, K, F_univar, A_S, b_S, A_T, b_T, x_star_vec)
        L.close()

        return {
            "ok_zero": True, "rankJ": rankJ, "rank_min": rank_min,
            "D_target": D, "D_used": D_used, "public_maxdeg": max(degs),
            "outfile": out_infile, "attempts": attempts,
            "genlog": genlog_path, "secretlog": secretlog_path
        }

# --------------------------- CLI ---------------------------
def main():
    if len(sys.argv) < 3:
        print("Usage: sage scripts/construct_HFE_fast.sage n D [outfile.in|-] [seed|-] [time_budget_sec] [max_maps] [secret_tries_per_map] [rank_min|-] [force(true/false)]")
        sys.exit(1)

    n = int(sys.argv[1]); D = int(sys.argv[2])

    # outfile
    if len(sys.argv) >= 4 and not is_unset(sys.argv[3]):
        outfile = sys.argv[3]
        if outfile.strip() == "-":
            outfile = os.path.join("data", "hfe_instances", f"HFE_n{n}_D{D}.in")
    else:
        outfile = os.path.join("data", "hfe_instances", f"HFE_n{n}_D{D}.in")

    # seed
    seed = None
    if len(sys.argv) >= 5 and not is_unset(sys.argv[4]):
        seed = int(sys.argv[4])

    time_budget_sec = 5400
    if len(sys.argv) >= 6 and not is_unset(sys.argv[5]):
        time_budget_sec = int(sys.argv[5])

    max_maps = 40
    if len(sys.argv) >= 7 and not is_unset(sys.argv[6]):
        max_maps = int(sys.argv[6])

    secret_tries_per_map = 2048
    if len(sys.argv) >= 8 and not is_unset(sys.argv[7]):
        secret_tries_per_map = int(sys.argv[7])

    rank_min = None
    if len(sys.argv) >= 9 and not is_unset(sys.argv[8]):
        rank_min = int(sys.argv[8])

    force = False
    if len(sys.argv) >= 10 and not is_unset(sys.argv[9]):
        force = str(sys.argv[9]).strip().lower() in {"1","true","t","yes","y"}

    ensure_dir_for(outfile)  # ensure folders exist for the .in file

    print(f"[{now_str()}] Generating HFE(n={n}, D={D}) → {outfile}")
    info = build_and_export_fast(
        n, D, outfile, seed=seed, time_budget_sec=time_budget_sec,
        max_maps=max_maps, secret_tries_per_map=secret_tries_per_map,
        rank_min=rank_min, force=force
    )
    if info.get("skipped"):
        print(f"[{now_str()}] Skipped (exists).")
    else:
        print(f"[{now_str()}] Done. rank={info['rankJ']} (target {info['rank_min']}), D_used={info['D_used']}, public_maxdeg={info['public_maxdeg']}, attempts={info['attempts']}")
        print(f"[{now_str()}] Wrote: {info['outfile']}")
        print(f"[{now_str()}] Gen log: {info['genlog']}")
        print(f"[{now_str()}] Secret log: {info['secretlog']}")

if __name__ == "__main__":
    main()
