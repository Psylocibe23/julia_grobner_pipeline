def log2_ceil(D):
    D = ZZ(D)
    if D <= 0:
        raise ValueError("D must be ≥ 1")
    if D == 1:
        return 0
    # ceil(log2 D) = floor(log2(D-1)) + 1 = (D-1).nbits()
    return (D - 1).nbits()

def dlf_bound_hfe(D):
    """Theorem bound for HFE over F_{2^n}."""
    t = log2_ceil(D)
    return t + 2

# demo
for D in [1,2,3,6,8,16,20,32,33,48,64,72,80,96]:
    print(f"D={D:>3}  t=ceil(log2 D)={log2_ceil(D)}  bound d_lf ≤ {dlf_bound_hfe(D)}")

