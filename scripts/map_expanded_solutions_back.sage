#!/usr/bin/env sage
import sys
import os

def parse_in_file(filename):
    with open(filename) as f:
        lines = [l.strip() for l in f if l.strip()]
    var_names = [v.strip() for v in lines[0].split(",")]
    field_line = lines[1]
    equations = lines[2:]
    if "^" in field_line:
        p, n = [int(x) for x in field_line.split("^")]
    else:
        p = int(field_line)
        n = 1
    return var_names, p, n

def read_solutions(filename):
    """Read {x_0: 1, x_1: 0, ...} style solutions from file."""
    sols = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            d = {}
            # Remove braces and spaces
            inner = line.strip("{} ")
            if inner:
                for part in inner.split(","):
                    k, v = part.split(":")
                    d[k.strip()] = int(v.strip())
            sols.append(d)
    return sols

def main():
    if len(sys.argv) != 4:
        print("Usage: sage scripts/map_expanded_solutions_back.sage <orig_infile> <expanded_infile> <expanded_solutions_file>")
        sys.exit(1)
    orig_infile, expanded_infile, solfile = sys.argv[1:4]
    orig_vars, p, n = parse_in_file(orig_infile)
    if n == 1:
        print("Original field is not an extension field; mapping not needed.")
        sys.exit(0)
    expanded_vars, _, _ = parse_in_file(expanded_infile)

    # Build field and polynomial ring
    Fpn = GF(p**n, name="a")
    a = Fpn.gen()
    Fp = GF(p)

    # Parse solutions
    sols = read_solutions(solfile)
    print(f"Read {len(sols)} solutions from expanded system.")

    # For each solution, reconstruct each original variable as x = x_0 + a*x_1 + ... + a^{n-1}*x_{n-1}
    mapped_solutions = []
    for sol in sols:
        mapped = {}
        for v in orig_vars:
            coeffs = [sol.get(f"{v}_{j}", 0) for j in range(n)]
            val = sum(Fpn(coeffs[j]) * a**j for j in range(n))
            mapped[v] = val
        mapped_solutions.append(mapped)

    # Write mapped solutions
    out_name = os.path.splitext(os.path.basename(solfile))[0] + "_MAPPED.txt"
    out_path = os.path.join("results", out_name)
    with open(out_path, "w") as f:
        f.write(f"# Mapped solutions for {solfile} (original vars over GF({p}^{n}))\n")
        for mapped in mapped_solutions:
            f.write("{" + ", ".join(f"{k}: {mapped[k]}" for k in orig_vars) + "}\n")
    print(f"Mapped solutions saved to {out_path}")

if __name__ == "__main__" or "sage" in __name__:
    main()
