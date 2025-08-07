################################################################################
# SAGE: PARALLELIZED HFE INSTANCE GENERATOR (QUADRATIC BOOLEAN, AT LEAST ONE SOLUTION)
################################################################################
# Generates an HFE system with guaranteed secret solution for cryptanalysis:
#   - n variables over GF(2)
#   - degree d HFE univariate poly over GF(2^n)
#   - applies secret invertible affine maps S, T
#   - outputs public system (quadratic Boolean, x_i^2 + x_i for all i)
#   - pipeline-compatible .in and log format
################################################################################

import os
import sys
from joblib import Parallel, delayed
import multiprocessing
from sage.all import *

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
        K.<a> = GF(q**n, modulus=prim_poly)
    else:
        K = GF(q**n, 'a')
        a = K.gen()
    mod_poly = K.modulus() if hasattr(K, "modulus") else None
    return K, K.gen(), mod_poly

def random_hfe_poly_quadratic(K, n, d, q=2):
    """
    Generate a random univariate HFE polynomial F(x) of degree at most d (≤ 2^n),
    with only terms of the form x^{q^i + q^j}, x^{q^k}, and constant.
    """
    R.<x> = K[]
    F = R(0)
    # Quadratic terms x^{q^i + q^j} (≤ d)
    for i in range(n):
        for j in range(i, n):
            exp = q**i + q**j
            if exp <= d:
                coeff = K.random_element()
                if coeff != 0:
                    F += coeff * x**exp
    # Linear terms x^{q^k} (≤ d)
    for k in range(n):
        exp = q**k
        if exp <= d:
            coeff = K.random_element()
            if coeff != 0:
                F += coeff * x**exp
    # Constant term (optional)
    const = K.random_element()
    if const != 0:
        F += const
    return F

def random_affine_map(n, Fp):
    """Generate a random invertible affine map (A, b) over Fp^n."""
    while True:
        A = random_matrix(Fp, n, n)
        if A.is_invertible():
            break
    b = vector(Fp, [Fp.random_element() for _ in range(n)])
    return A, b

def extract_coordinate_polys(fy, a, n, xvars):
    """
    Convert an element of GF(2^n)[xvars] to an n-tuple of polynomials over GF(2).
    """
    K = a.parent()
    Fq = K.base_ring()
    R = xvars[0].parent()
    result = [R(0) for _ in range(n)]
    for monom, coeff in fy.dict().items():
        coeff_vec = K(coeff)._vector_()
        for idx in range(n):
            if coeff_vec[idx] != 0:
                mono_poly = R({monom: 1})
                result[idx] += Fq(coeff_vec[idx]) * mono_poly
    return result

def compute_public_poly(idx, F, A_S, b_S, A_T, b_T, K, a, xvars):
    """
    Compute the idx-th coordinate of the public HFE map as a multivariate polynomial over GF(2).
    """
    Fp = GF(2)
    R = xvars[0].parent()
    n = len(xvars)
    s = A_S * vector(R, xvars) + vector(R, list(b_S))
    s_field = sum([s[i] * a**i for i in range(n)])
    fy = F(s_field)
    fy_vec = extract_coordinate_polys(fy, a, n, xvars)
    fy_vec = vector(R, fy_vec)
    z = A_T * fy_vec + vector(R, list(b_T))
    return z[idx]

def is_solution(poly_list, secret_vec, xvars):
    """
    Given a list of Boolean polynomials and a secret vector, check if all vanish at secret.
    """
    subst = {str(x): int(v) for x, v in zip(xvars, secret_vec)}
    for poly in poly_list:
        if poly(**subst) != 0:
            return False
    return True

def reduce_to_quadratic_boolean(poly, xvars):
    """
    Reduce a multivariate polynomial modulo x_i^2 = x_i and keep only degree ≤2 terms.
    Returns: BooleanPolynomial (in B) of degree at most 2.
    """
    B = BooleanPolynomialRing(len(xvars), names=[str(v) for v in xvars])
    # Coerce to Boolean polynomial (automatic reduction mod field equations)
    poly_bool = B(str(poly))
    # Only keep terms of degree at most 2
    result = B.zero()
    for mon in poly_bool.monomials():
        if mon.degree() <= 2:
            coeff = poly_bool.monomial_coefficient(mon)
            if coeff != 0:
                result += coeff * mon
    return result

def export_public_system_quadratic(n, q, F, A_S, b_S, A_T, b_T, K, a, outfilename, secret, num_jobs=None):
    """
    Generate quadratic Boolean public polynomials, guarantee secret is a root,
    and export to file.
    """
    Fp = GF(q)
    R = PolynomialRing(Fp, n, 'x')
    xvars = R.gens()
    B = BooleanPolynomialRing(n, names=[str(v) for v in xvars])

    print(f"Computing {n} public polynomials in parallel using {num_jobs or multiprocessing.cpu_count()} cores...")

    # Compute each output coordinate polynomial (public key), in parallel
    public_polys = Parallel(n_jobs=num_jobs or multiprocessing.cpu_count())(
        delayed(compute_public_poly)(i, F, A_S, b_S, A_T, b_T, K, a, xvars) for i in range(n)
    )

    # Ensure known solution by evaluating P(secret)
    s_secret = A_S * secret + b_S
    s_field_secret = sum([s_secret[i] * a**i for i in range(n)])
    fy_secret = F(s_field_secret)
    fy_vec_secret = vector(Fp, K(fy_secret)._vector_())
    public_out = A_T * fy_vec_secret + b_T

    # Shifted, reduce to quadratic Boolean
    quad_public_polys = []
    for i in range(n):
        poly = public_polys[i] - public_out[i]
        poly_quad = reduce_to_quadratic_boolean(poly, xvars)
        quad_public_polys.append(poly_quad)

    # --- Write output file (pipeline format: no comments!) ---
    os.makedirs(os.path.dirname(outfilename), exist_ok=True)
    with open(outfilename, "w") as f:
        # Variable names (comma-separated, no comment)
        f.write(", ".join(str(v) for v in xvars) + "\n")
        # Field characteristic
        f.write(f"{q}\n")
        # Public quadratic Boolean polynomials
        for poly in quad_public_polys:
            f.write(str(poly) + "\n")
        # Boolean field equations
        for v in xvars:
            f.write(f"{str(v)}^2 + {str(v)}\n")

    print(f"Exported quadratic Boolean system to {outfilename}")
    return quad_public_polys, xvars  # for solution checking

################################################################################
# MAIN: Parse arguments, orchestrate, guarantee secret solution
################################################################################

if __name__ == "__main__":
    # Parse command-line arguments: n, d, [num_threads]
    if len(sys.argv) < 3:
        print("Usage: sage construct_HFE_system.sage <n> <d> [num_threads]")
        print("  <n>: number of variables (extension degree)")
        print("  <d>: max degree of secret HFE polynomial")
        print("  [num_threads]: number of CPU threads (default: 8)")
        sys.exit(1)

    q = 2
    n = int(sys.argv[1])
    d = int(sys.argv[2])
    num_threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8

    set_random_seed(42)

    # 1. Construct extension field
    K, a, modulus = construct_extension_field(q, n)
    print(f"Field: GF({q}^{n}), modulus: {modulus}")

    Fp = GF(q)
    attempts = 0
    found = False

    while not found:
        attempts += 1
        # 2. Generate random HFE polynomial of degree ≤ d
        F = random_hfe_poly_quadratic(K, n, d, q)
        print("Random HFE polynomial F(X):")
        print(F)

        # 3. Random secret affine maps
        A_S, b_S = random_affine_map(n, Fp)
        A_T, b_T = random_affine_map(n, Fp)
        print("Affine input map S(x) = A_S * x + b_S:")
        print(A_S)
        print(b_S)
        print("Affine output map T(y) = A_T * y + b_T:")
        print(A_T)
        print(b_T)

        # 4. Choose secret
        secret = vector(Fp, [Fp.random_element() for _ in range(n)])
        print(f"Secret vector: {secret}")

        # 5. Output paths
        outname = f"data/hfe_instances/HFE_n{n}_d{d}_system.in"
        logname = f"logs/HFE_n{n}_d{d}_instance_info.txt"
        os.makedirs("logs", exist_ok=True)

        # 6. Export system and get quad_public_polys
        quad_public_polys, xvars = export_public_system_quadratic(
            n, q, F, A_S, b_S, A_T, b_T, K, a, outname, secret, num_jobs=num_threads
        )

        # 7. Check if secret is a root of the quadratic Boolean system
        if is_solution(quad_public_polys, list(secret), xvars):
            print(f"Secret is a root! Found after {attempts} attempt(s).")
            found = True
        else:
            print(f"Secret not a root, regenerating... (attempt {attempts})")

    # 8. Log parameters (unchanged format)
    with open(logname, "w") as logf:
        logf.write(f"HFE Instance Info (n={n}, d={d})\n")
        logf.write(f"Field: GF({q}^{n}), modulus: {modulus}\n\n")
        logf.write(f"F(X) = {F}\n\n")
        logf.write("Affine input map (A_S, b_S):\n")
        logf.write(str(A_S) + "\n" + str(b_S) + "\n\n")
        logf.write("Affine output map (A_T, b_T):\n")
        logf.write(str(A_T) + "\n" + str(b_T) + "\n\n")
        logf.write(f"Secret = {list(secret)}\n")
    print(f"Logged system info to {logname}")
    print("HFE instance generation complete.")
