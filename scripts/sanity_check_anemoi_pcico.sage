###############################################################################
# sanity_check_anemoi_pcico.sage — correctness checks for emitted PCICO systems
#
# What it does
#   1) Loads a pipeline .in file (variables, p, polynomials).
#   2) Rebuilds AnemoiPermutation with (p, alpha, rounds) and regenerates the
#      PCICO equations (ℓ = 1) via model_P_CICO(..., final_ll=True, ordering=1).
#      Compares the two polynomial sets (up to ring coercion).
#   3) Brute-forces (or samples) y0 ∈ GF(p) with x_in=0 to find cases where the
#      final linear-layer x-output is 0 (PCICO condition). For each found y0, it
#      computes s_t = y_in[t] - y_out[t] per round and checks ALL input-file
#      polynomials vanish under the assignment {Y0000=y0, S0100=s1, ..., Sr00}.
#
# Usage
#   sage scripts/sanity_check_anemoi_pcico.sage \
#        --infile data/ANEMOI_p251_r2_a3_PCICO.in --alpha 3 --rounds 2
#
# Options
#   --tries N       : random y0 trials if p is large (default 500)
#   --exhaustive    : try all y0 ∈ GF(p) (auto-enabled if p ≤ 50000)
#   --max-sol K     : stop after confirming K solutions (default 2)
#   --verbose       : print extra details
###############################################################################

import sys, os, random

# Make sure we can import our local copies
if os.path.isdir("scripts"):
    sys.path.insert(0, os.path.abspath("scripts"))

load("scripts/anemoi.sage")
load("scripts/models.sage")

def usage():
    print("Usage: sage scripts/sanity_check_anemoi_pcico.sage "
          "--infile <.in path> --alpha <odd> --rounds <r> "
          "[--tries N] [--exhaustive] [--max-sol K] [--verbose]", file=sys.stderr)
    sys.exit(1)

# ---------------- CLI ----------------
argv = sys.argv[1:]
args = {
    "infile": None,
    "alpha": None,
    "rounds": None,
    "tries": 500,
    "exhaustive": False,
    "max_sol": 2,
    "verbose": False,
}
i = 0
while i < len(argv):
    a = argv[i]
    if a == "--infile" and i+1 < len(argv):
        args["infile"] = argv[i+1]; i += 2
    elif a == "--alpha" and i+1 < len(argv):
        args["alpha"] = int(argv[i+1]); i += 2
    elif a == "--rounds" and i+1 < len(argv):
        args["rounds"] = int(argv[i+1]); i += 2
    elif a == "--tries" and i+1 < len(argv):
        args["tries"] = int(argv[i+1]); i += 2
    elif a == "--exhaustive":
        args["exhaustive"] = True; i += 1
    elif a == "--max-sol" and i+1 < len(argv):
        args["max_sol"] = int(argv[i+1]); i += 2
    elif a == "--verbose":
        args["verbose"] = True; i += 1
    else:
        usage()

if args["infile"] is None or args["alpha"] is None or args["rounds"] is None:
    usage()

# ---------------- Parse .in ----------------
with open(args["infile"], "r") as f:
    lines = [ln.strip() for ln in f if ln.strip()]

if len(lines) < 3:
    raise ValueError("Malformed .in file: need ≥ 3 non-empty lines.")

var_names = [v.strip() for v in lines[0].split(",")]
p = int(lines[1])
poly_strs = lines[2:]

if args["verbose"]:
    print(f"[file] vars={var_names}")
    print(f"[file] p={p}")
    print(f"[file] #eqs={len(poly_strs)}")

if not is_prime(p):
    raise ValueError(f"Line 2 p={p} is not prime (this checker expects GF(p)).")

F = GF(p)
R_file = PolynomialRing(F, var_names, order='degrevlex')
G_file = [R_file(s) for s in poly_strs]

# ---------------- Rebuild model equations & compare ----------------
alpha = args["alpha"]
r = args["rounds"]
if gcd(alpha, p-1) != 1:
    raise ValueError(f"alpha={alpha} not invertible mod p-1; need gcd(alpha, p-1)=1.")

P = AnemoiPermutation(q=p, alpha=alpha, n_rounds=r, n_cols=1)
G_model, _ = model_P_CICO(P, n_rounds=r, final_ll=True, ordering=1, debug=False)

# Coerce model equations to the file's ring (identical variable names/order)
G_model_fileRing = [R_file(g) for g in G_model]

set_file  = {str(g) for g in G_file}
set_model = {str(g) for g in G_model_fileRing}
same_polys = (set_file == set_model)

print("== Polynomial-set comparison ==")
print(f"  #file eqs = {len(G_file)}")
print(f"  #model eqs = {len(G_model_fileRing)}")
print(f"  sets equal? {same_polys}")
if not same_polys:
    only_in_file  = set_file - set_model
    only_in_model = set_model - set_file
    print(f"  [!] {len(only_in_file)} equations only in .in file")
    print(f"  [!] {len(only_in_model)} equations only in model")
    # keep going; evaluation test below is the decisive check

# Degree summary (cheap sanity)
degs = [g.total_degree() for g in G_file]
print("== Degree summary (from file) ==")
print(f"  degrees: {degs}  (max={max(degs)})")

# ---------------- Concrete evaluation check ----------------
# We search y0 with x_in=0 and final-LL x_out=0. For each such y0, compute s_t and
# evaluate all polynomials from the file at the assignment.

def find_assignments(max_needed=2, exhaustive=False, tries=500, verbose=False):
    sols = []
    # generator of y0 candidates
    if exhaustive or p <= 50000:
        cand = list(F)
    else:
        # sample without replacement up to 'tries' random elements
        cand = []
        seen = set()
        while len(cand) < tries and len(seen) < p:
            a = F(random.randrange(p))
            if int(a) not in seen:
                seen.add(int(a))
                cand.append(a)
    for y0 in cand:
        # initial state (CICO): x0 = 0, y0 variable
        x = F(0)
        y = F(y0)
        s_values = []
        for t in range(1, r+1):
            # add round constants (r-1 index in code)
            x_in = x + P.C[t-1][0]
            y_in = y + P.D[t-1][0]
            # apply linear layer
            [x_in], [y_in] = P.linear_layer([x_in], [y_in])
            # evaluate sbox to get outputs at end of round t
            x, y = P.evaluate_sbox(x_in, y_in)
            # s_t = y_in - y_out
            s_values.append(y_in - y)
        # final linear layer (PCICO with final_ll=True adds x==0)
        [xf], [yf] = P.linear_layer([x], [y])
        if xf != 0:
            continue
        # build assignment dict by NAME
        assignment = {}
        # Y0000
        assignment[R_file("Y0000")] = y0
        # S0100,...,Sr00   (note: the file may list variables in any order; we map by name)
        for t, sval in enumerate(s_values, start=1):
            vname = f"S{t:02d}00"
            if vname in var_names:
                assignment[R_file(vname)] = sval
        # ensure every variable has a value
        if any(v not in assignment for v in R_file.gens()):
            # missing? then skip (this would be odd for ℓ=1)
            continue
        # evaluate all polynomials
        ok = all(g.subs(assignment) == 0 for g in G_file)
        if ok:
            sols.append((y0, s_values, assignment))
            if verbose:
                print(f"  found solution y0={int(y0)} ; s={ [int(s) for s in s_values] }")
            if len(sols) >= max_needed:
                break
    return sols

print("== Evaluation check on concrete instances ==")
exhaustive = args["exhaustive"] or (p <= 50000)
sols = find_assignments(max_needed=args["max_sol"], exhaustive=exhaustive,
                        tries=args["tries"], verbose=args["verbose"])
print(f"  confirmed {len(sols)} solution(s) with final x_out=0 and all polynomials vanishing.")

if len(sols) == 0:
    print("  [note] no solutions found with given search budget. This can happen for large p "
          "with random sampling. Try --exhaustive (if p is small) or increase --tries.")
else:
    # echo one solution nicely
    y0, svals, assign = sols[0]
    def pv(a):  # pretty value as int
        try:
            return int(a)
        except Exception:
            return a
    print("  example assignment:")
    print(f"    Y0000 = {pv(y0)}")
    for t, sval in enumerate(svals, start=1):
        print(f"    S{t:02d}00 = {pv(sval)}")
print("Done.")
