################################################################################
# solve_F4_from_file_parsing.jl — F4 Gröbner basis (DRL) with streaming parser
################################################################################

using AlgebraicSolving              # F4, Ideal, normal_form, etc.
import AbstractAlgebra               # MPolyBuildCtx
import Nemo                          # GF(p)
using Dates, Printf
using Base.Filesystem: mkpath, basename, splitext

# ---------------- CLI ----------------
function parse_args()
    if length(ARGS) < 1
        println("Usage: julia solve_F4_from_file_parsing.jl <inputfile> [nthreads]")
        exit(1)
    end
    filename = ARGS[1]
    nthreads = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 1
    nthreads ≥ 1 || error("nthreads must be ≥ 1 (got $nthreads).")
    return filename, nthreads
end

filename, nthreads = parse_args()
mkpath("logs"); mkpath("results")

# ---------------- Input (.in) reader ----------------
# Line 1: comma-separated variable names
# Line 2: prime p
# Line 3+: polynomials (expanded sums of monomials; no parentheses)
function isprime64(n::Integer)::Bool
    n ≤ 1 && return false
    n ≤ 3 && return true
    iseven(n) && return false
    for p in (3,5,7,11,13,17,19,23,29,31,37)
        if n == p; return true; end
        if n % p == 0; return false; end
    end
    d = n - 1; s = 0
    while iseven(d); d ÷= 2; s += 1; end
    function modexp(a,e,m)
        a%=m; r=1%m; ee=e
        while ee>0
            (ee&1)==1 && (r=(r*a)%m)
            a=(a*a)%m; ee>>=1
        end
        r
    end
    function is_witness(a)
        x = modexp(a,d,n)
        (x==1 || x==n-1) && return false
        for _ in 1:(s-1)
            x = (x*x) % n
            x == n-1 && return false
        end
        true
    end
    for a in (2,3,5,7,11,13,17)
        a % n == 0 && continue
        is_witness(a) && return false
    end
    true
end

function parse_input_system(path::AbstractString)
    lines = readlines(path)
    length(lines) ≥ 3 || error("Input file must have ≥ 3 lines.")
    # vars
    var_line = strip(lines[1]); var_line ≠ "" || error("Line 1 (variables) is empty.")
    var_names = [strip(v) for v in split(var_line, ",")]
    any(isempty, var_names) && error("Empty variable name in Line 1.")
    # field
    field_spec = strip(lines[2])
    occursin("^", field_spec) && error("GF(p^n) not supported in this stage.")
    p = try parse(Int, field_spec) catch; error("Bad field characteristic: '$field_spec'"); end
    (p ≥ 2 && isprime64(p)) || error("Field characteristic must be prime (got $p).")
    # polys
    poly_strs = String[]
    for (j, line) in enumerate(lines[3:end])
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        push!(poly_strs, s)
    end
    !isempty(poly_strs) || error("No polynomials after line 2.")
    return var_names, p, poly_strs
end

var_names, p, poly_strs = parse_input_system(filename)

# ---------------- Ring ----------------
K = Nemo.GF(p)
R, vars = AbstractAlgebra.polynomial_ring(K, var_names; internal_ordering = :degrevlex)
# Bind variables into Main so printer matches input names (not strictly needed)
for (i, v) in enumerate(var_names)
    @eval Main $(Symbol(v)) = vars[$i]
end

# ---------------- Streaming polynomial parser ----------------
# Grammar per term:   [±] [int coeff]* [*] (var[^exp] [*])*    until next ±
# Assumes the .in is fully expanded (no parentheses).
function parse_poly_line_build(s::String,
                               R::AbstractAlgebra.MPolyRing,
                               K, var_index::Dict{String,Int}, p::Int)
    s = replace(s, " " => "")
    nvars = length(var_index)
    exps  = zeros(Int, nvars)
    ctx   = AbstractAlgebra.MPolyBuildCtx(R)

    i  = firstindex(s); L = lastindex(s)
    while i <= L
        # sign
        sign = 1
        if s[i] == '+'
            i += 1
        elseif s[i] == '-'
            sign = -1; i += 1
        end
        # reset term
        c = 1 % p
        fill!(exps, 0)

        # factors until next +/-
        while i <= L && s[i] != '+' && s[i] != '-'
            if isdigit(s[i])
                j = i
                while j <= L && isdigit(s[j]); j += 1; end
                val = parse(Int, s[i:j-1]) % p
                c = (c * val) % p
                i = j
                if i <= L && s[i] == '*'; i += 1; end
            else
                # variable token
                j = i
                while j <= L && (isletter(s[j]) || isdigit(s[j]) || s[j] == '_'); j += 1; end
                vname = s[i:j-1]
                idx = get(var_index, vname, 0)
                idx == 0 && error("Unknown variable '$vname' in input polynomial.")
                e = 1
                if j <= L && s[j] == '^'
                    j += 1
                    k = j
                    while k <= L && isdigit(s[k]); k += 1; end
                    e = parse(Int, s[j:k-1])
                    i = k
                else
                    i = j
                end
                exps[idx] += e
                if i <= L && s[i] == '*'; i += 1; end
            end
        end

        coeff = mod(sign * c, p)
        if coeff != 0
            AbstractAlgebra.push_term!(ctx, K(coeff), copy(exps))  # copy exps, builder stores by ref
        end
        # do not consume the +/- here; outer loop handles it
    end

    return AbstractAlgebra.finish(ctx)
end

function parse_polynomials(poly_strs::Vector{String})
    P = Vector{typeof(vars[1])}()
    vmap = Dict(v => i for (i, v) in enumerate(var_names))
    t0 = time()
    for (k, s) in enumerate(poly_strs)
        try
            push!(P, parse_poly_line_build(s, R, K, vmap, p))
        catch err
            @error "Error while parsing polynomial #$k (input line $(k+2))" s
            rethrow(err)
        end
        # tiny heartbeat for very long inputs
        (k % 50 == 0) && @info "Parsed $k / $(length(poly_strs)) polynomials..."
    end
    @info @sprintf("Parsed %d polynomials in %.3f s", length(poly_strs), time() - t0)
    return P
end

polys = parse_polynomials(poly_strs)

# ---------------- Ideal & F4 ----------------
I = Ideal(polys)

input_id  = splitext(basename(filename))[1]
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
log_file    = joinpath("logs",    "$(input_id)_F4_$timestamp.log")
result_file = joinpath("results", "$(input_id)_F4_$timestamp.txt")

open(log_file, "w") do logio
    redirect_stdout(logio) do
        redirect_stderr(logio) do
            println("================================================================")
            println(" solve_F4_from_file_parsing.jl — F4 Gröbner basis (DRL, GF($p))")
            println("================================================================")
            println("Input file: $filename")
            println("Variables: ", join(var_names, ", "))
            println("Order: degrevlex (DRL)")
            println("Threads requested: $nthreads")
            println("Number of input polynomials: ", length(polys))
            println("----------------------------------------------------------------")
            println("Starting F4 (AlgebraicSolving.groebner_basis)...")
            t0 = time()
            G = nothing
            try
                @time G = groebner_basis(I; la_option=44, nr_thrds = nthreads, info_level = 2, initial_hts = 22, max_nr_pairs = 30000)
            catch err
                println("\nF4 FAILED:\n$err")
                rethrow(err)
            end
            t1 = time()
            @printf("F4 wall time: %.3f seconds\n", t1 - t0)
            println("Basis size: ", length(G))

            # (Optional) correctness: each input generator reduces to 0 mod Ideal(G)
            println("Verifying generators reduce to 0 modulo Ideal(G)...")
            J = Ideal(G)
            for (idx, f) in enumerate(polys)
                nf = normal_form(f, J)
                nf == 0 || error("Generator #$idx did not reduce to 0.\n  f  = $f\n  nf = $nf")
            end
            println("Verification passed.")

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
            println("Log file:    $log_file")
            println("Result file: $result_file")
            println("================================================================")
        end
    end
end

println("Computation done. Log: $log_file   Result: $result_file")
