#!/usr/bin/env sage
###############################################################################
# check_expansion_correctness.sage
#
# This script verifies the correctness of an *expansion* from a system over an
# extension field GF(p^n) (using n variables of the base field GF(p) to represent
# each extension variable) to a system entirely over GF(p). It checks, by enumeration,
# that for each assignment to the extension field variables, the solutions to the
# original and expanded systems agree.
###############################################################################

import sys
import itertools

###############################################################################
# 1. Parse a system file (variables, field, equations)
###############################################################################
def parse_input_file(filename):
    """
    Parse a polynomial system file:
        - First line: variable names (comma separated)
        - Second line: field (either p or p^n)
        - Rest: polynomial equations (as strings)
    Returns: (variable names, base p, extension degree n, equations as strings)
    """
    with open(filename) as f:
        lines = [l.strip() for l in f if l.strip()]
    var_names = [v.strip() for v in lines[0].split(",")]
    field_line = lines[1]
    equations = lines[2:]
    # Parse field
    if "^" in field_line:
        base, ext = field_line.split("^")
        p = int(base.strip())
        n = int(ext.strip())
    else:
        p = int(field_line)
        n = 1
    return var_names, p, n, equations

###############################################################################
# 2. Get the vector space representation for elements of GF(p^n) over GF(p)
###############################################################################
def get_vector_space(Fpn):
    """
    For a given extension field Fpn = GF(p^n), return the vector space over GF(p)
    corresponding to its elements. This allows us to represent any element of Fpn
    as an n-tuple over Fp (using a fixed basis).
    """
    vs = Fpn.vector_space()
    if isinstance(vs, tuple):
        # Sometimes returns (field, VS)
        V = vs[-1]
    else:
        V = vs
    return V

###############################################################################
# 3. Main expansion correctness check
###############################################################################
def verify_expansion(orig_file, exp_file, max_check=100000):
    """
    For each assignment to the original system's variables over GF(p^n),
    compute the corresponding base field coordinates, substitute into the
    expanded system, and check if the solution status matches.

    Stops after max_check assignments for large search spaces.
    """
    # Parse original (extension field) system
    orig_vars, p, n, orig_polys = parse_input_file(orig_file)
    exp_var_names, p_exp, n_exp, exp_polys = parse_input_file(exp_file)
    print(f"Original system: {len(orig_vars)} variables over GF({p}^{n}), {len(orig_polys)} equations.")
    print(f"Expanded system: {len(exp_var_names)} variables over GF({p_exp}), {len(exp_polys)} equations.")

    if n == 1:
        print("Original system is already in base field. Nothing to check.")
        return

    Fpn = GF(p**n, name="a")
    Fp  = GF(p)
    a = Fpn.gen()

    # Polynomial rings
    PR_K = PolynomialRing(Fpn, orig_vars)
    PR_p = PolynomialRing(Fp, exp_var_names)
    orig_poly_objs = [PR_K(s) for s in orig_polys]
    exp_poly_objs  = [PR_p(s) for s in exp_polys]

    # Prepare coordinate map
    V = get_vector_space(Fpn)
    total = 0
    match = 0

    print("Testing all assignments...\n")

    # Assignments to original variables
    all_assignments = itertools.product(Fpn, repeat=len(orig_vars))
    for orig_vals in all_assignments:
        subst_orig = dict(zip(orig_vars, orig_vals))
        # Expand to base field coordinates
        base_coords = []
        for val in orig_vals:
            coords = list(V(val))
            base_coords.extend(int(c) for c in coords)
        subst_exp = dict(zip(exp_var_names, base_coords))

        # Check both systems
        ok_orig = all(poly(**subst_orig) == 0 for poly in orig_poly_objs)
        ok_exp  = all(poly(**subst_exp) == 0 for poly in exp_poly_objs)
        if ok_orig != ok_exp:
            print(f"Discrepancy detected for assignment:\n  orig: {subst_orig}\n  exp: {subst_exp}\n")
            print(f"  Orig system: {'OK' if ok_orig else 'FAIL'}; Expanded: {'OK' if ok_exp else 'FAIL'}")
            sys.exit(1)
        if ok_orig and ok_exp:
            match += 1
        total += 1
        if total % 10000 == 0 and total > 0:
            print(f"Checked {total} assignments...")

        if total >= max_check:
            print(f"Reached max_check={max_check}. Stopping.")
            break

    print(f"\nChecked {total} assignments. Matches: {match}.")
    print("Expansion is CORRECT: original and expanded systems agree for all assignments checked.")

###############################################################################
# 4. Main entry point
###############################################################################
def main():
    if len(sys.argv) != 3:
        print("Usage: sage scripts/check_expansion_correctness.sage <original.in> <expanded.in>")
        sys.exit(1)
    orig_file = sys.argv[1]
    exp_file = sys.argv[2]
    verify_expansion(orig_file, exp_file)

if __name__ == "__main__" or "sage" in __name__:
    main()
