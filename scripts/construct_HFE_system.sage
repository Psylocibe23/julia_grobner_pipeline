################################################################################
# SAGE: PARALLELIZED HFE INSTANCE GENERATOR (AT LEAST ONE SOLUTION)
################################################################################
# Generates an HFE system:
#   - n variables over GF(2)
#   - degree d univariate polynomial over GF(2^n)
#   - applies secret affine maps S (input), T (output)
#   - guarantees at least one solution (the secret)
#   - exports system in .in format (for Gröbner attack)
################################################################################

import os
from joblib import Parallel, delayed
import multiprocessing
from sage.all import *

def hamming_weight(k):
    return bin(k).count("1")

def construct_extension_field(q, n, prim_poly=None):
    Fq = GF(q)
    if prim_poly is not None:
        K.<a> = GF(q**n, modulus=prim_poly)
    else:
        K = GF(q**n, 'a')
        a = K.gen()
    mod_poly = K.modulus() if hasattr(K, "modulus") else None
    return K, K.gen(), mod_poly

def random_hfe_poly_quadratic(K, n, d, q=2):
    R.<x> = K[]
    F = R(0)
    for i in range(n):
        for j in range(i, n):
            exp = q**i + q**j
            if exp <= d:
                coeff = K.random_element()
                if coeff != 0:
                    F += coeff * x**exp
    for k in range(n):
        exp = q**k
        if exp <= d:
            coeff = K.random_element()
            if coeff != 0:
                F += coeff * x**exp
    const = K.random_element()
    if const != 0:
        F += const
    return F

def random_affine_map(n, Fp):
    while True:
        A = random_matrix(Fp, n, n)
        if A.is_invertible():
            break
    b = vector(Fp, [Fp.random_element() for _ in range(n)])
    return A, b

def extract_coordinate_polys(fy, a, n, xvars):
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

def export_public_system_parallel(n, q, F, A_S, b_S, A_T, b_T, K, a, outfilename, secret, num_jobs=None):
    Fp = GF(q)
    R = PolynomialRing(Fp, n, 'x')
    xvars = R.gens()

    print(f"Computing {n} public polynomials in parallel using {num_jobs or multiprocessing.cpu_count()} cores...")

    public_polys = Parallel(n_jobs=num_jobs or multiprocessing.cpu_count())(
        delayed(compute_public_poly)(i, F, A_S, b_S, A_T, b_T, K, a, xvars) for i in range(n)
    )

    # Ensure known solution by evaluating P(secret)
    s_secret = A_S * secret + b_S
    s_field_secret = sum([s_secret[i] * a**i for i in range(n)])
    fy_secret = F(s_field_secret)
    fy_vec_secret = vector(Fp, K(fy_secret)._vector_())
    public_out = A_T * fy_vec_secret + b_T

    os.makedirs(os.path.dirname(outfilename), exist_ok=True)
    with open(outfilename, "w") as f:
        f.write("# Variables: " + ", ".join(str(v) for v in xvars) + "\n")
        f.write(f"# Field: GF({q})\n")
        for idx in range(n):
            f.write(str(public_polys[idx] - public_out[idx]) + "\n")

    print(f"Exported system to {outfilename}")

################################################################################
# MAIN: Edit parameters here
################################################################################

if __name__ == "__main__":
    # Set to HFE(80,96) for full system
    q = 2
    n = 3       
    d = 3        
    num_threads = 16  # Set to number of cores on your machine

    set_random_seed(42)

    # 1. Field GF(2^n)
    K, a, modulus = construct_extension_field(q, n)
    print(f"Field: GF({q}^{n}), modulus: {modulus}")

    # 2. Random HFE polynomial
    F = random_hfe_poly_quadratic(K, n, d, q)
    print("Random HFE polynomial F(X):")
    print(F)

    # 3. Secret affine maps
    Fp = GF(q)
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

    # 6. Export system
    export_public_system_parallel(n, q, F, A_S, b_S, A_T, b_T, K, a, outname, secret, num_jobs=num_threads)

    # 7. Log parameters
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
