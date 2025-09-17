###############################################################################
# emit_anemoi_pcico.sage — Write PCICO (l=1) Anemoi systems over GF(p) to .in
#
# INPUT  (CLI)
#   --p prime p, or hex literal, or a name from constants.py
#   --alpha <odd> S-box alpha (must satisfy gcd(alpha, p-1) = 1)
#   --rounds <r> number of rounds r (l=1 only)
#   [--outfile path] output file path (default: data/ANEMOI_p<p>_r<r>_a<alpha>_PCICO.in)
#   [--ordering {1,2,3}] variable ordering in model_P_CICO (default 1 = recommended)
#   [--no-final-ll] omit final linear layer constraint (default: include it)
#e.g.: sage scripts/emit_anemoi_pcico.sage --p 251 --alpha 3 --rounds 2
#
# OUTPUT (.in)
#   Line 1: comma-separated variable names (in the exact DRL tie-break order we use)
#   Line 2: prime p
#   Line 3+: one polynomial per line (Sage textual form; coefficients in [0..p-1])
#
# NOTES
#   • Model: P_CICO, l=1 (as in the paper’s best-performing prime-field model)
#   • Field: odd prime GF(p) only 
###############################################################################

import sys, os


if os.path.isdir("scripts"):
    sys.path.insert(0, os.path.abspath("scripts"))

# Load modules
load("scripts/anemoi.sage")
load("scripts/models.sage")

# Import constants.py as a Python module to resolve names given to --p
try:
    import importlib
    constants_mod = importlib.import_module("constants")  # from scripts/constants.py
except Exception:
    constants_mod = None

# ================= CLI parsing =================
def _usage_and_exit():
    print("Usage: sage scripts/emit_anemoi_pcico.sage --p <prime|0xHEX|NAME> --alpha <odd> --rounds <r> "
          "[--outfile path] [--ordering 1|2|3] [--no-final-ll]", file=sys.stderr)
    sys.exit(1)

argv = sys.argv[1:]
args = {"p": None, "alpha": None, "rounds": None,
        "outfile": None, "ordering": 1, "final_ll": True}

i = 0
while i < len(argv):
    a = argv[i]
    if a == "--p" and i+1 < len(argv):
        args["p"] = argv[i+1]; i += 2
    elif a == "--alpha" and i+1 < len(argv):
        args["alpha"] = int(argv[i+1]); i += 2
    elif a == "--rounds" and i+1 < len(argv):
        args["rounds"] = int(argv[i+1]); i += 2
    elif a == "--outfile" and i+1 < len(argv):
        args["outfile"] = argv[i+1]; i += 2
    elif a == "--ordering" and i+1 < len(argv):
        args["ordering"] = int(argv[i+1]); i += 2
    elif a == "--no-final-ll":
        args["final_ll"] = False; i += 1
    else:
        _usage_and_exit()

if args["p"] is None or args["alpha"] is None or args["rounds"] is None:
    _usage_and_exit()

# ================= Convert/validate p =================
def resolve_prime(p_spec):
    """
    p_spec may be decimal int, 0x... hex, or a symbol in constants.py.
    Returns p as Python int (Sage Integer is fine too).
    """
    # hex literal
    if isinstance(p_spec, str) and p_spec.lower().startswith("0x"):
        try:
            return int(p_spec, 16)
        except Exception:
            pass
    # decimal int
    try:
        return int(p_spec)
    except Exception:
        pass
    # symbol in constants.py
    if constants_mod is not None and isinstance(p_spec, str):
        if hasattr(constants_mod, p_spec):
            return int(getattr(constants_mod, p_spec))
    raise ValueError(f"Could not resolve prime from: {p_spec}")

p = resolve_prime(args["p"])
alpha = int(args["alpha"])
r = int(args["rounds"])

if not is_prime(p):
    raise ValueError(f"Provided p = {p} is not prime (GF(p) required by the pipeline).")

if gcd(alpha, p-1) != 1:
    raise ValueError(f"alpha = {alpha} is not invertible modulo p-1; need gcd(alpha, p-1) = 1.")

if args["ordering"] not in (1,2,3):
    raise ValueError("ordering must be one of {1,2,3} (recommended: 1)")

# ================= Build the PCICO system =================
# l = 1 (n_cols = 1)
P = AnemoiPermutation(q=p, alpha=alpha, n_rounds=r, n_cols=1)

F_CICO, _dict_vars = model_P_CICO(
    P, n_rounds=r, final_ll=args["final_ll"], ordering=args["ordering"], debug=False
)

if not F_CICO:
    raise RuntimeError("Model returned an empty equation set; this should not happen.")

# Variable names in the EXACT order of the DRL ring produced by the model
R = F_CICO[0].parent()
var_names = list(R.variable_names())

# Sanity: expecting r + (1 if final_ll else 0) equations for l=1
expected_eqs = r + (1 if args["final_ll"] else 0)
if len(F_CICO) != expected_eqs:
    print(f"[warn] Got {len(F_CICO)} equations; expected {expected_eqs}. Proceeding.", file=sys.stderr)

# ================= Write the .in file =================
if args["outfile"] is None:
    stem = f"ANEMOI_p{p}_r{r}_a{alpha}_PCICO" + ("" if args["final_ll"] else "_noFLL")
    outdir = "data"
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, stem + ".in")
else:
    outfile = args["outfile"]
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)

with open(outfile, "w") as f:
    # Line 1: variables (comma-separated, NO spaces)
    f.write(",".join(var_names) + "\n")
    # Line 2: field prime
    f.write(str(p) + "\n")
    # Lines 3+: equations
    for g in F_CICO:
        # Ensure coefficients are printed as integers in [0..p-1]
        f.write(str(g) + "\n")

print(f"Wrote PCICO system to: {outfile}")
print(f"Variables (|V|={len(var_names)}): {', '.join(var_names)}")
print(f"Field: GF({p})  (alpha={alpha}, rounds={r}, ordering={args['ordering']}, final_ll={args['final_ll']})")
print(f"#Equations: {len(F_CICO)}")
