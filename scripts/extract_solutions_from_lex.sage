# fast_extract_zerosdim_solutions.sage

import sys

def parse_basis_file(filename):
    variables, p, polys = None, None, []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                if line.startswith("# Variables:"):
                    variables = [v.strip() for v in line.split(":",1)[1].split(",")]
                elif line.startswith("# Field: GF("):
                    field_str = line.split("GF(",1)[1].split(")")[0]
                    if "^" in field_str:
                        p, n = [int(x.strip()) for x in field_str.split("^")]
                        p = p**n
                    else:
                        p = int(field_str)
                continue
            polys.append(line)
    if variables is None or p is None:
        raise ValueError("Variables or field not found")
    return variables, p, polys

def reduce_exponents(poly, p):
    R = poly.parent()
    new_monomials = {}
    for mon, coeff in poly.dict().items():
        new_mon = tuple(e % p for e in mon)
        if new_mon in new_monomials:
            new_monomials[new_mon] += coeff
        else:
            new_monomials[new_mon] = coeff
    # Remove zeros
    new_monomials = {mon: coeff for mon, coeff in new_monomials.items() if coeff != 0}
    return R(new_monomials)

def main():
    if len(sys.argv) < 2:
        print("Usage: sage fast_extract_zerosdim_solutions.sage <basis_file>")
        sys.exit(1)
    filename = sys.argv[1]
    variables, p, polys_str = parse_basis_file(filename)
    F = GF(p)
    R = PolynomialRing(F, variables)
    polys = [reduce_exponents(R(s), p) for s in polys_str if s]
    I = R.ideal(polys)
    print(f"Extracting solutions for {filename} using state-of-the-art methods.")
    sols = I.variety()
    print(f"Found {len(sols)} solutions.")
    out_file = filename.replace('.txt', '_sols.txt')
    with open(out_file, 'w') as out:
        for sol in sols:
            out.write("{" + ", ".join(f"{k}: {v}" for k, v in sol.items()) + "}\n")
    print(f"Solutions saved to {out_file}")

if __name__ == "__main__":
    main()
