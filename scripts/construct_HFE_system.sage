################################################################################
# SAGE: HFE INSTANCE GENERATOR (QUADRATIC, LOG SECRET, OUTPUT PUBLIC SYSTEM)
################################################################################
# Generates an HFE system:
#   - n variables (public key size)
#   - degree d (max degree of the secret univariate poly)
#   - only quadratic, linear, and constant terms (as in classic HFE)
#   - secret affine maps S (input), T (output), secret F
#   - Saves secret parameters to log file
#   - Saves public system as n polynomials in n variables (for algebraic attack)
################################################################################

import os

def hamming_weight(k):
    """Return the Hamming weight (number of 1's) in binary of integer k."""
    return bin(k).count("1")

def construct_extension_field(q, n, prim_poly=None):
    """
    Construct GF(q^n) with generator 'a', optionally with a given irreducible polynomial.
    Returns: (K, a, modulus_poly)
    """
    Fq = GF(q)
    if prim_poly is not None:
        K.<a> = GF(Fq**n, modulus=prim_poly)
        return K, a, prim_poly
    else:
        K = GF(q**n, 'a')
        a = K.gen()
        mod_poly = K.modulus() if hasattr(K, "modulus") else None
        return K, a, mod_poly

def random_hfe_poly_quadratic(K, n, d, q=2):
    """
    Return a random univariate HFE polynomial F(x) of degree at most d in K[x]
    with only terms of the form x^{q^i + q^j}, x^{q^k}, and 1,
    i.e., only quadratic, linear, and constant terms.
    """
    R.<x> = K[]
    F = R(0)
    # Quadratic terms: exponents q^i + q^j, for all 0 <= i <= j < n, <= d
    for i in range(n):
        for j in range(i, n):
            exp = q**i + q**j
            if exp > d:
                continue
            coeff = K.random_element()
            if coeff != 0:
                F += coeff * x**exp
    # Linear terms: x^{q^k}
    for k in range(n):
        exp = q**k
        if exp > d:
            continue
        coeff = K.random_element()
        if coeff != 0:
            F += coeff * x**exp
    # Constant term
    coeff = K.random_element()
    if coeff != 0:
        F += coeff
    return F

def random_affine_map(n, Fp):
    """
    Generate a random invertible affine map (A, b) over Fp^n.
    Returns: (A, b) where A is invertible nxn matrix, b is vector.
    """
    while True:
        A = random_matrix(Fp, n, n)
        if A.is_invertible():
            break
    b = vector(Fp, [Fp.random_element() for _ in range(n)])
    return (A, b)

def extract_coordinate_polys(fy, a, n, xvars):
    """
    Given fy: a polynomial over K[x_0,...,x_{n-1}],
    express it as a vector of n polynomials over Fq[x_0,...,x_{n-1}],
    corresponding to the canonical basis 1, a, ..., a^{n-1}.
    """
    K = a.parent()
    Fq = K.base_ring()
    R = xvars[0].parent()
    result = [R(0) for _ in range(n)]  # polynomials for each coordinate
    for monom, coeff in fy.dict().items():
        # coeff is in K, expand in basis
        coeff_vec = K(coeff)._vector_()  # <-- THIS IS THE FIX!
        for idx in range(n):
            if coeff_vec[idx] != 0:
                mono_poly = R({monom: 1})
                result[idx] += Fq(coeff_vec[idx]) * mono_poly
    return result





def export_symbolic_public_map(n, q, F, A_S, b_S, A_T, b_T, K, a, outfilename):
    """
    Export the public HFE map as n explicit polynomials in n variables,
    as expected by algebraic cryptanalysis tools.
    """
    Fp = GF(q)
    R = PolynomialRing(Fp, n, 'x')
    xvars = R.gens()
    var_names = [str(v) for v in xvars]
    basis = [a^i for i in range(n)]

    print(f"Exporting public key: {n} polynomials in {n} variables...")

    public_polys = []
    for idx in range(n):
        xvec = list(xvars)
        s = A_S * vector(R, xvec) + vector(R, list(b_S))
        s_field = sum([s[i] * a**i for i in range(n)])
        fy = F(s_field)
        # << THIS IS THE KEY FIX >>
        fy_vec = extract_coordinate_polys(fy, a, n, xvars)
        fy_vec = vector(R, fy_vec)
        z = A_T * fy_vec + vector(R, list(b_T))
        public_polys.append(z[idx])

    os.makedirs(os.path.dirname(outfilename), exist_ok=True)
    with open(outfilename, "w") as f:
        f.write(", ".join(var_names) + "\n")
        f.write(str(q) + "\n")
        for poly in public_polys:
            f.write(str(poly) + "\n")
    print(f"Exported system to {outfilename}")







################################################################################
# MAIN LOGIC
################################################################################
if __name__ == "__main__":
    # --- User-set parameters (adjust as needed) ---
    q = 2             # Field characteristic (2 for binary HFE)
    n = 4             # Number of variables 
    d = 4             # Maximum degree of secret univariate polynomial F
    # ------------------------------------------------

    # --- 1. Construct extension field and generator ---
    K, a, modulus = construct_extension_field(q, n)
    print(f"Field: GF({q}^{n}), primitive element 'a'")
    print(f"Modulus polynomial: {modulus}")

    # --- 2. Generate random HFE polynomial F (quadratic, degree ≤ d) ---
    F = random_hfe_poly_quadratic(K, n, d, q)
    print(f"Random HFE poly F(x) (degree ≤ {d}):")
    print(F)

    # --- 3. Generate random affine maps S, T (secret) ---
    Fp = GF(q)
    A_S, b_S = random_affine_map(n, Fp)
    A_T, b_T = random_affine_map(n, Fp)
    print("Affine input map S(x) = A_S * x + b_S:")
    print(A_S)
    print(b_S)
    print("Affine output map T(y) = A_T * y + b_T:")
    print(A_T)
    print(b_T)

    # --- 4. Save all secret instance information (for later verification/debug) ---
    os.makedirs("logs", exist_ok=True)
    log_filename = f"logs/HFE_n{n}_d{d}_instance_info.txt"
    with open(log_filename, "w") as logf:
        logf.write(f"=== HFE Instance Information (n={n}, d={d}) ===\n\n")
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
    print(f"Saved secret info to {log_filename}")

    # --- 5. Export public system for algebraic cryptanalysis (pipeline input) ---
    outname = f"data/hfe_instances/HFE_n{n}_d{d}_system.in"
    export_symbolic_public_map(n, q, F, A_S, b_S, A_T, b_T, K, a, outname)

    print("HFE system generation complete.")
