################################################################################
# solve_F4_from_file.jl  —  F4 Gröbner basis stage (prime fields, DRL)
#
# PURPOSE
#   Read a polynomial system over GF(p) from a text file, construct the ideal,
#   compute a Gröbner basis with F4 (via AlgebraicSolving.jl), verify it, and
#   write the basis (with metadata) to disk for the next pipeline stage (FGLM).
#
# KEY DESIGN CHOICES
#   • Field: PRIME fields GF(p) ONLY in this stage (no extensions).
#   • Monomial order: DRL = degrevlex (explicitly recorded in the output).
#   • Correctness: Each input generator is reduced modulo Ideal(G)
#
# INPUT FILE FORMAT
#   Line 1: comma-separated variables, e.g. "x0, x1, x2"
#   Line 2: field "p" 
#   Line 3+: one polynomial per line, in the variables declared on line 1
#
# USAGE
#   julia solve_F4_from_file.jl <inputfile> [nthreads]
#
# OUTPUT
#   results/<name>_F4_<timestamp>.txt: Gröbner basis + metadata
#   logs/<name>_F4_<timestamp>.log: full computation log/trace
#
################################################################################

using AlgebraicSolving             
using Dates                        
using Printf                       
using Base.Filesystem: mkpath, basename, splitext

# ================ Utils ================
function isprime64(n::Integer)::Bool
    n ≤ 1 && return false
    n ≤ 3 && return true           
    iseven(n) && return false      
    for p in (3,5,7,11,13,17,19,23,29,31,37)
        if n == p; return true; end
        if n % p == 0; return false; end
    end

    # Write n - 1 = d * 2^s with d odd
    d = n - 1
    s = 0
    while iseven(d)
        d ÷= 2
        s += 1
    end

    # Fast modular exponentiation (binary exponentiation) with Int arithmetic
    function modexp(a::Integer, e::Integer, m::Integer)
        a %= m
        r = 1 % m
        ee = e
        while ee > 0
            if (ee & 1) == 1
                r = (r * a) % m
            end
            a = (a * a) % m
            ee ÷= 2
        end
        return r
    end

    function is_witness(a::Integer)
        x = modexp(a, d, n)
        (x == 1 || x == n - 1) && return false
        for _ in 1:(s-1)
            x = (x * x) % n
            x == n - 1 && return false
        end
        return true
    end

    # Deterministic base set for 64-bit integers
    for a in (2, 3, 5, 7, 11, 13, 17)
        a % n == 0 && continue  # if base coincides with n, skip
        if is_witness(a)
            return false        
        end
    end
    return true # prime
end


function parse_args()
    if length(ARGS) < 1
        println("Usage: julia solve_F4_from_file.jl <inputfile> [nthreads]")
        exit(1)
    end
    filename = ARGS[1]
    nthreads = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 1
    nthreads ≥ 1 || error("nthreads must be ≥ 1 (got $nthreads).")
    return filename, nthreads
end

filename, nthreads = parse_args()

mkpath("logs")
mkpath("results")


function parse_input_system(path::AbstractString)
    lines = readlines(path)
    length(lines) ≥ 3 || error("Input file must have ≥ 3 lines: variables, field, at least one polynomial.")

    # ---- Line 1: Variables
    var_line = strip(lines[1])
    var_line ≠ "" || error("Line 1 (variables) is empty.")
    var_names = [strip(v) for v in split(var_line, ",")]
    any(isempty, var_names) && error("Empty variable name found on line 1.")

    # ---- Line 2: Field (GF(p))
    field_spec = strip(lines[2])
    occursin("^", field_spec) && error("GF(p^n) not supported in this stage. Use a prime field GF(p).")
    # Parse p (as Int) and check primality 
    p = try
        parse(Int, field_spec)
    catch
        error("Could not parse field characteristic from line 2: '$field_spec'.")
    end
    p ≥ 2 || error("Field characteristic must be ≥ 2 (got $p).")
    isprime64(p) || error("Field characteristic p must be prime (got $p).")

    # ---- Lines 3+: Polynomials (skip blanks and comment lines)
    poly_strs = String[]
    for (j, line) in enumerate(lines[3:end])
        _lineno = j + 2                      
        s = strip(line)
        if isempty(s) || startswith(s, "#")
            continue                         
        end
        push!(poly_strs, s)
    end
    !isempty(poly_strs) || error("No polynomials found after line 2.")

    return var_names, p, poly_strs
end

var_names, p, poly_strs = parse_input_system(filename)

# ================ Build ring and ideal ================
# We work in the finite prime field GF(p)
K = GF(p)
R, vars = polynomial_ring(K, var_names; internal_ordering = :degrevlex)
for (i, v) in enumerate(var_names)
    @eval Main $(Symbol(v)) = vars[$i]
end

function parse_polynomials(poly_strs::Vector{String})
    P = Vector{typeof(vars[1])}() # vector of polynomials in R
    for (k, s) in enumerate(poly_strs)
        try
            # Evaluate the string using the bound symbols x0, x1, ...
            push!(P, eval(Meta.parse(s)))
        catch err
            @error "Error while parsing polynomial #$k (input line $(k+2)):" s
            rethrow(err)
        end
    end
    return P
end

polys = parse_polynomials(poly_strs)

# Form the ideal I = ⟨polys⟩ in R
I = Ideal(polys)


input_id  = splitext(basename(filename))[1]
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
log_file    = joinpath("logs",    "$(input_id)_F4_$timestamp.log")
result_file = joinpath("results", "$(input_id)_F4_$timestamp.txt")

# ================ F4 & Outputs ================
open(log_file, "w") do logio
    # Redirect both STDOUT and STDERR into the log for a compact trace.
    redirect_stdout(logio) do
        redirect_stderr(logio) do
            println("================================================================")
            println(" solve_F4_from_file.jl — F4 Gröbner basis computation (DRL)")
            println("================================================================")
            println("Input file: $filename")
            println("Variables: ", join(var_names, ", "))
            println("Field: GF($p)")
            println("Assumed order: degrevlex (DRL)")
            println("Threads requested: $nthreads")
            println("Number of input polynomials: ", length(polys))
            println("----------------------------------------------------------------")

            # F4 computation
            println("Starting F4 computation (AlgebraicSolving.groebner_basis)...")
            t0 = time()
            G = nothing
            try
                # info_level=2 provides verbosity inside the F4 engine
                @time G = groebner_basis(I; nr_thrds = nthreads, info_level = 2)
            catch err
                println("\nF4 computation FAILED with error:\n$err")
                rethrow(err)
            end
            t1 = time()
            @printf("F4 wall time: %.3f seconds\n", t1 - t0)
            println("Computed Gröbner basis size: ", length(G))

            # Correctness check
            # Verify that G is a Gröbner basis for I by reducing each original
            # generator f in polys modulo J = Ideal(G); all normal forms must be 0
            println("Verifying correctness (reduce each input generator modulo Ideal(G))...")
            J = Ideal(G)
            for (idx, f) in enumerate(polys)
                r = normal_form(f, J)  # normal form of f modulo G
                if r != 0
                    println("FAILED generator #$idx: nf ≠ 0\n  f   = $f\n  nf  = $r")
                    error("Correctness check failed: generator #$idx did not reduce to 0.")
                end
            end
            println("Correctness check PASSED.")

            # Write result file 
            println("----------------------------------------------------------------")
            println("Writing Gröbner basis to: $result_file")
            open(result_file, "w") do fio
                println(fio, "# Groebner basis (F4) computed for $filename")
                println(fio, "# Variables: ", join(var_names, ", "))
                println(fio, "# Field: GF($p)")
                println(fio, "# Order: degrevlex")
                println(fio, "# Basis: reduced=UNKNOWN")
                println(fio, "# Number of input equations: $(length(polys))")
                println(fio, "# --- Groebner basis ---")
                for g in G
                    println(fio, g)
                end
            end

            println("Done.")
            println("----------------------------------------------------------------")
            println("Log file:    $log_file")
            println("Result file: $result_file")
            println("================================================================")
        end
    end
end

println("Computation done. Log: $log_file   Result: $result_file")
