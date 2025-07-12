import random
import os 
from sage.crypto.boolean_function import BooleanFunction


def construct_extension_field(q, n, prim_poly=None):
    """
    Construct the extension field F_{q^n}
    If prim_poly is None use a random irreducible polynomial of degree n
    Returns (K,a, prim_poly)=(F_{q^n}, primitive element, irreducible polynomial)
    """
    Fq = GF(q)
    if prim_poly is not None:
        # prim_poly is an irreducible polynomial over F_q of degree n
        K.<a> = GF(Fq**n, modulus=prim_poly)
        return K, a, prim_poly
    else:
        # Let Sage choose the Conway polynomial
        K = GF(q**n, 'a')
        a = K.gen()
        # Get modulus
        PR = PolynomialRing(Fq, 'z')
        mod_poly = K.modulus() if hasattr(K, "modulus") else None
        return K, a, mod_poly


def random_hfe_poly(K, n):
    """
    Generate a random HFE polynomial F(x) over K = GF(2^n) with degree up to 2^n.
    """
    R.<x> = K[]
    F = R(0)
    # Quadratic terms
    for i in range(n):
        for j in range(i, n):
            coeff = K.random_element()
            if coeff != 0:  
                exp = 2^i + 2^j
                F += coeff * x^exp
    # Linear terms
    for k in range(n):
        coeff = K.random_element()
        if coeff != 0:
            F += coeff * x^(2^k)
    # Constant term
    F += K.random_element()
    return F


def random_affine_map(n, Fp):
    A = random_matrix(Fp, n, n)
    while not A.is_invertible():
        A = random_matrix(Fp, n, n)
    b = vector(Fp, [Fp.random_element() for _ in range(n)])
    return (A, b)


def vec_to_field(xvec):
    return sum([int(xvec[i]) * a**i for i in range(n)])

def field_to_vec(val):
    V = Fpn.vector_space()[1]
    return V(val)

def S_map(xvec):
    return A_S * vector(Fp, [Fp(x) for x in xvec]) + b_S

def T_map(yvec):
    return A_T * vector(Fp, [Fp(x) for x in yvec]) + b_T


def public_map(xvec):
    # 1. S
    xs = S_map(xvec)
    # 2. Embed to F_{2^n}
    xs_field = vec_to_field(xs)
    # 3. Apply HFE poly
    y_field = F(xs_field)
    # 4. Back to F_2^n
    y_vec = field_to_vec(y_field)
    # 5. T
    z = T_map(y_vec)
    return z


def export_public_map_to_infile(n, public_map, outfilename):
    """
    For each output coordinate, constructs its algebraic normal form (ANF) as a polynomial in x_0,...,x_{n-1}.
    Saves to the .in file format:
        Line 1: x_0, ..., x_{n-1}
        Line 2: 2
        Next n lines: polynomials for each output coordinate
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




if __name__=="__main__":
    q = 2
    n = 4

    # 1. Construct extension field
    K, a, modulus = construct_extension_field(q, n)
    print(f"Field: GF({q}^{n}) with primitive element 'a'.")
    print(f"Modulus polynomial: {modulus}")
    print(f"a^({q**n}) = {a**(q**n)} (should equal a for the correct minimal polynomial)")

    # 2. Random HFE polynomial
    F = random_hfe_poly(K, n)
    print(f"Random HFE poly F(X): {F}")

    # 3. Random affine maps S, T
    Fp = GF(q)
    A_S, b_S = random_affine_map(n, Fp)
    A_T, b_T = random_affine_map(n, Fp)
    print("Affine input map S(x) = A_S * x + b_S")
    print(A_S)
    print(b_S)
    print("Affine output map T(y) = A_T * y + b_T")
    print(A_T)
    print(b_T)

    # 4. Global variables for functions
    Fpn = K

    # 5. Test public_map on random input
    xvec = [Fp.random_element() for _ in range(n)]
    print(f"\nSample input x = {xvec}")
    z = public_map(xvec)
    print(f"Public map output (z): {list(z)}")

    # 6. (Optional: test on all inputs if n is small)
    if n <= 4:
        print("\nAll possible inputs and outputs:")
        from itertools import product
        for xv in product([0,1], repeat=n):
            out = public_map(list(xv))
            print(f"x={list(xv)} --> z={list(out)}")

    outname = f"data/hfe_instances/HFE_n{n}_system.in"
    export_public_map_to_infile(n, public_map, outname)
