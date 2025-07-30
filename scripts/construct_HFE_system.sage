import os
from random import randint

################################################################################
# HFE SYSTEM GENERATOR FOR LARGE n (SYMBOLIC, NO ENUMERATION)
################################################################################
#
# This script constructs a random instance of a general HFE (Hidden Field Equations)
# system, ready for algebraic cryptanalysis.
# - For large n (e.g., n=80), it **symbolically expands** the public key polynomials,
#   writing them explicitly as multivariate polynomials in the public variables.
#
# For small n (e.g., n <= 5), you can still use the old enumeration/interpolation.
################################################################################

def construct_extension_field(q, n, prim_poly=None):
    """
    Construct the extension field F_{q^n} used for HFE.
    - q: base field characteristic (usually 2)
    - n: extension degree
    - prim_poly: (optional) irreducible polynomial for reproducibility
    Returns: (K, a, prim_poly)
        K: extension field GF(q^n)
        a: primitive element (generator)
        prim_poly: modulus/irreducible polynomial
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

def random_hfe_poly(K, n, d=None, q=2):
    """
    Generate a random HFE polynomial F(x) over K = GF(q^n), degree at most d.
    Only nonzero terms up to degree d are included.
    """
    R.<x> = K[]
    F = R(0)
    # Quadratic terms: exponents q^i + q^j
    for i in range(n):
        for j in range(i, n):
            exp = q**i + q**j
            if (d is not None) and (exp > d):
                continue
            coeff = K.random_element()
            if coeff != 0:
                F += coeff * x**exp
    # Linear terms: exponents q^k
    for k in range(n):
        exp = q**k
        if (d is not None) and (exp > d):
            continue
        coeff = K.random_element()
        if coeff != 0:
            F += coeff * x**exp
    # Constant term
    F += K.random_element()
    return F

def random_affine_map(n, Fp):
    """
    Generate a random invertible affine map over Fp^n:
      - Matrix A ∈ GL(n, Fp)
      - Vector b ∈ Fp^n
    Used for S (input transformation) and T (output).
    """
    A = random_matrix(Fp, n, n)
    while not A.is_invertible():
        A = random_matrix(Fp, n, n)
    b = vector(Fp, [Fp.random_element() for _ in range(n)])
    return (A, b)

def vec_to_field(xvec, a):
    """
    Embed a vector over Fp^n as an element of GF(p^n) using the canonical basis.
    """
    return sum([int(xvec[i]) * a**i for i in range(len(xvec))])

def field_to_vec(val, Fpn, n):
    """
    Convert a GF(p^n) element to its vector of coordinates in Fp^n.
    """
    return vector(GF(2), Fpn.polynomial()(val).list() + [0]*(n - len(Fpn.polynomial()(val).list())))

def export_symbolic_public_map(n, q, F, A_S, b_S, A_T, b_T, K, a, outfilename):
    """
    Symbolically expand and export the public HFE map as n explicit multivariate polynomials
    in the public variables, for large n (e.g., n=80).
    """
    # Set up the base field and polynomial ring for n public variables
    Fp = GF(q)
    R = PolynomialRing(Fp, n, 'x')
    xvars = R.gens()
    var_names = [str(v) for v in xvars]

    print(f"Exporting symbolic HFE public key as {n} polynomials in {n} variables...")

    public_polys = []
    for idx in range(n):
        # 1. Input as symbolic vector of variables
        xvec_sym = list(xvars)
        # 2. S: affine input transformation
        sx = A_S * vector(Fp, xvec_sym) + b_S
        # 3. Embed into extension field
        sx_field = sum([sx[i] * a**i for i in range(n)])
        # 4. Evaluate HFE polynomial (symbolically)
        fy = F(sx_field)
        # 5. Convert back to vector over Fp^n (coordinate-wise, with respect to basis)
        # (Use K.polynomial() to get vector space representation)
        fy_coeffs = K.polynomial()(fy).list()
        fy_vec = vector(Fp, fy_coeffs + [0]*(n - len(fy_coeffs)))
        # 6. Output affine map
        z = A_T * fy_vec + b_T
        # 7. Each coordinate z[idx] is a symbolic polynomial in x_0,...,x_{n-1}
        public_polys.append(z[idx].polynomial())

    # Write to .in file
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
    # --- USER-ADJUSTABLE PARAMETERS ------------------------------------------
    q = 2      # Field characteristic
    n = 80     # Number of variables (extension degree)
    d = 96     # Degree bound for HFE polynomial

    # --- 1. Field construction -----------------------------------------------
    K, a, modulus = construct_extension_field(q, n)
    print(f"Field: GF({q}^{n}) with primitive element 'a'.")
    print(f"Modulus polynomial: {modulus}")
    print(f"a^({q**n}) = {a**(q**n)} (should equal a for the correct minimal polynomial)")

    # --- 2. Random HFE polynomial generation ---------------------------------
    F = random_hfe_poly(K, n, d)
    print(f"Random HFE poly F(X): {F}")

    # --- 3. Random affine maps S, T ------------------------------------------
    Fp = GF(q)
    A_S, b_S = random_affine_map(n, Fp)
    A_T, b_T = random_affine_map(n, Fp)
    print("Affine input map S(x) = A_S * x + b_S")
    print(A_S)
    print(b_S)
    print("Affine output map T(y) = A_T * y + b_T")
    print(A_T)
    print(b_T)

    # --- 3b. Save secret HFE instance information for debugging --------------
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

    # --- 4. Export public system to .in file for algebraic attacks -----------
    outname = f"data/hfe_instances/HFE_n{n}_system.in"
    export_symbolic_public_map(n, q, F, A_S, b_S, A_T, b_T, K, a, outname)

    print("HFE system generation complete.")
