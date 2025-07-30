import random
import os 
from sage.crypto.boolean_function import BooleanFunction

################################################################################
# HFE SYSTEM GENERATOR FOR CRYPTANALYSIS
################################################################################
#
# This script constructs a random instance of a general HFE (Hidden Field Equations)
# system, ready for algebraic cryptanalysis.
# HFE is a multivariate public-key cryptosystem proposed by Patarin, built as follows:
# - It uses a univariate polynomial of special form F(X) over an extension field GF(q^n)
# - The public key is a multivariate system (over GF(q)) obtained by composing F(X)
#   with secret affine maps S (input), T (output), and embedding variables as field elements.
# - The security reduction relies on the difficulty of solving the resultant multivariate
#   system in the public variables (typically quadratic, over GF(2)).
#
# This script implements:
#   1. Field construction (including primitive/irreducible polynomials)
#   2. Random HFE polynomial generation
#   3. Random secret affine maps S, T
#   4. Public key map construction and evaluation
#   5. Exporting the multivariate system (as Boolean polynomials) for attack pipelines
################################################################################

# --- 1. Field Construction ----------------------------------------------------
def construct_extension_field(q, n, prim_poly=None):
    """
    Construct the extension field F_{q^n} used for HFE.
    - q: base field characteristic (usually 2)
    - n: extension degree (number of variables in public system)
    - prim_poly: if provided, an irreducible polynomial over F_q of degree n
    Returns: (K, a, prim_poly)
        K: the extension field GF(q^n)
        a: the primitive element (generator)
        prim_poly: the modulus/irreducible polynomial (for reproducibility)
    """
    Fq = GF(q)
    if prim_poly is not None:
        # prim_poly is an irreducible polynomial over F_q of degree n
        K.<a> = GF(Fq**n, modulus=prim_poly)
        return K, a, prim_poly
    else:
        # Use Sage's built-in field constructor (Conway polynomial by default).
        K = GF(q**n, 'a')
        a = K.gen()
        # Get modulus
        PR = PolynomialRing(Fq, 'z')
        mod_poly = K.modulus() if hasattr(K, "modulus") else None
        return K, a, mod_poly

# --- 2. Random HFE Polynomial ------------------------------------------------
def random_hfe_poly(K, n, d=None, q=2):
    """
    Generate a random HFE polynomial F(x) over K = GF(q^n), degree at most d (if specified).
    If d is None, allows maximal degree.
    """
    R.<x> = K[]
    F = R(0)
    exponents_seen = set()

    # Quadratic terms: exponents of the form q^i + q^j
    for i in range(n):
        for j in range(i, n):
            exp = q**i + q**j
            if (d is not None) and (exp > d):
                continue
            coeff = K.random_element()
            if coeff != 0:
                F += coeff * x^exp
                exponents_seen.add(exp)
    # Linear terms: exponents of the form q^k
    for k in range(n):
        exp = q**k
        if (d is not None) and (exp > d):
            continue
        coeff = K.random_element()
        if coeff != 0:
            F += coeff * x^exp
            exponents_seen.add(exp)
    # Constant term (always present)
    F += K.random_element()
    return F


# --- 3. Random (Secret) Affine Maps S, T -------------------------------------
def random_affine_map(n, Fp):
    """
    Generate a random invertible affine map over Fp^n:
      - Matrix A ∈ GL(n, Fp)
      - Vector b ∈ Fp^n
    Used for both S (input transformation) and T (output).
    Secret in HFE; critical for hiding the algebraic structure from attackers.
    """
    A = random_matrix(Fp, n, n)
    while not A.is_invertible():
        A = random_matrix(Fp, n, n)
    b = vector(Fp, [Fp.random_element() for _ in range(n)])
    return (A, b)

# --- 4. Coordinate and Field Conversion Utilities -----------------------------
def vec_to_field(xvec):
    """
    Given a vector over Fp of length n, interpret it as an element in GF(p^n)
    using the canonical basis (1, a, a^2, ..., a^{n-1}) for the extension field.
    This is the standard embedding for HFE construction and attacks.
    """
    return sum([int(xvec[i]) * a**i for i in range(n)])

def field_to_vec(val):
    """
    Convert a GF(p^n) element to its vector of coordinates in Fp^n
    (with respect to the canonical basis). Required for switching between
    field element representation and public variable representation.
    """
    V = Fpn.vector_space()[1]
    return V(val)

# --- 5. Secret S, T, HFE Poly Setup (assigned globally in __main__) ----------
def S_map(xvec):
    """
    Secret input affine map: S(x) = A_S * x + b_S
    (x is a vector in Fp^n)
    """
    return A_S * vector(Fp, [Fp(x) for x in xvec]) + b_S

def T_map(yvec):
    """
    Secret output affine map: T(y) = A_T * y + b_T
    """
    return A_T * vector(Fp, [Fp(x) for x in yvec]) + b_T


def public_map(xvec):
    """
    The full public HFE map: z = T(F(S(x)))
    Steps:
        1. Apply input affine map S to xvec ∈ Fp^n
        2. Embed S(x) as a field element in GF(p^n)
        3. Evaluate F at this field element (univariate HFE poly over extension field)
        4. Expand result back to Fp^n coordinates
        5. Apply output affine map T
    Returns: z ∈ Fp^n (vector)
    """
    # 1. S: secret affine transformation of input
    xs = S_map(xvec)

    # 2. Embed S(x) into the extension field as a single element
    xs_field = vec_to_field(xs)

    # 3. Evaluate the HFE polynomial F
    y_field = F(xs_field)

    # 4. Expand the output field element to a vector over Fp
    y_vec = field_to_vec(y_field)

    # 5. T: secret output affine transformation
    z = T_map(y_vec)
    return z

# --- 6. Export to Public Multivariate System File ----------------------------
def export_public_map_to_infile(n, public_map, outfilename):
    """
    For cryptanalysis: enumerate the truth table of the HFE public map,
    interpolate each output coordinate as a Boolean polynomial (Algebraic Normal Form),
    and save the system in standard ".in" format for pipelines.
      Line 1: variable names (x_0,...,x_{n-1})
      Line 2: field (2 for GF(2))
      Next n lines: ANF polynomials for each output variable z_i in terms of x_j

    The ANF is essential: it expresses each z_i as an explicit (possibly quadratic)
    multivariate polynomial in the public variables. This is the input for
    algebraic attacks and Gröbner basis computations.
    """
    Fp = GF(2)
    R = BooleanPolynomialRing(n, 'x')
    xvars = R.gens()
    var_names = [str(v) for v in xvars]

    polys = []
    # For each output coordinate (z_i)
    for idx in range(n):
        # Build the truth table for the ith coordinate as a function of x_0,...,x_{n-1}
        table = []
        for val in cartesian_product_iterator([Fp]*n):
            xvec = vector(Fp, val)
            zout = public_map(xvec)
            table.append(zout[idx])
        # Use interpolation (ANF) to get the polynomial for z_i
        p = BooleanFunction(table).algebraic_normal_form()
        polys.append(str(p))

    # Write .in file
    with open(outfilename, "w") as f:
        f.write(", ".join(var_names) + "\n")
        f.write("2\n")
        for poly in polys:
            f.write(poly + "\n")
    print(f"Exported system to {outfilename}")



################################################################################
# MAIN LOGIC
################################################################################
if __name__=="__main__":
    # --- USER-ADJUSTABLE PARAMETERS ------------------------------------------
    q = 2  # Usually HFE is used over GF(2), but can generalize
    n = 80  # Number of variables (and extension degree)
    d = 96
    # For cryptanalysis research, n = 3..8 for experiments; real-world HFE uses larger n


    # --- 1. Field construction -----------------------------------------------
    K, a, modulus = construct_extension_field(q, n)
    print(f"Field: GF({q}^{n}) with primitive element 'a'.")
    print(f"Modulus polynomial: {modulus}")
    print(f"a^({q**n}) = {a**(q**n)} (should equal a for the correct minimal polynomial)")

    # --- 2. Random HFE polynomial generation --------------------------------
    F = random_hfe_poly(K, n, d)
    print(f"Random HFE poly F(X): {F}")

    # --- 3. Random affine maps S, T -----------------------------------------
    Fp = GF(q)
    A_S, b_S = random_affine_map(n, Fp)
    A_T, b_T = random_affine_map(n, Fp)
    print("Affine input map S(x) = A_S * x + b_S")
    print(A_S)
    print(b_S)
    print("Affine output map T(y) = A_T * y + b_T")
    print(A_T)
    print(b_T)

    # --- 3b. Save secret HFE instance information for debugging ---
    os.makedirs("logs", exist_ok=True)
    log_filename = f"logs/HFE_n{n}_system_information.txt"
    with open(log_filename, "w") as logf:
        logf.write(f"=== HFE Instance Information (n={n}) ===\n\n")
        logf.write("Field:\n")
        logf.write(f"  GF({q}^{n}), primitive element 'a'\n")
        logf.write(f"  Modulus polynomial: {modulus}\n\n")
        logf.write("Secret HFE Polynomial F(X):\n")
        logf.write(f"  F(X) = {F}\n\n")
        logf.write("Secret Affine Map S(x) = A_S * x + b_S:\n")
        logf.write(f"A_S =\n{A_S}\n")
        logf.write(f"b_S = {b_S}\n\n")
        logf.write("Secret Affine Map T(y) = A_T * y + b_T:\n")
        logf.write(f"A_T =\n{A_T}\n")
        logf.write(f"b_T = {b_T}\n")
    print(f"Saved HFE system information to {log_filename}")

    # --- 4. Define global variables for use in functions above ---------------
    Fpn = K

    # --- 5. Test public map on a random input --------------------------------
    xvec = [Fp.random_element() for _ in range(n)]
    print(f"\nSample input x = {xvec}")
    z = public_map(xvec)
    print(f"Public map output (z): {list(z)}")

    # --- 6. Enumerate all possible inputs (if small n) ----------------------
    if n <= 4:
        print("\nAll possible inputs and outputs:")
        from itertools import product
        for xv in product([0,1], repeat=n):
            out = public_map(list(xv))
            print(f"x={list(xv)} --> z={list(out)}")

    # --- 7. Export public system to .in file for algebraic attacks -----------
    outname = f"data/hfe_instances/HFE_n{n}_system.in"
    export_public_map_to_infile(n, public_map, outname)
