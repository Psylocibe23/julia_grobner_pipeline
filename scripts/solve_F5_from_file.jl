using AlgebraicSolving
import Base.Filesystem: mkpath
using Dates

# --- Setup output folders ---
mkpath("logs")
mkpath("results")

# --- Parse command-line arguments ---
function parse_args()
    if length(ARGS) < 1
        println("Usage: julia solve_F4_from_file.jl <inputfile> [nthreads]")
        exit(1)
    end
    filename = ARGS[1]
    return filename
end

filename = parse_args()
lines = readlines(filename)
var_names = [strip(v) for v in split(strip(lines[1]), ",")]
p = parse(Int, strip(lines[2]))

# --- Field and polynomial ring ---
if p == 2
    K = GF(2)
else
    K = GF(p)
end

R, vars = polynomial_ring(K, var_names)

# --- Define variables in Main module for eval ---
for (i, v) in enumerate(var_names)
    @eval Main $(Symbol(v)) = vars[$i]
end

# --- Parse polynomials ---
polys = Vector{typeof(vars[1])}()
for line in lines[3:end]
    try
        poly = eval(Meta.parse(line))
        push!(polys, poly)
    catch err
        println("Error evaluating polynomial: $line")
        rethrow(err)
    end
end

input_id = splitext(basename(filename))[1]
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
log_file = "logs/$(input_id)_F5_$timestamp.log"
result_file = "results/$(input_id)_F5_$timestamp.txt"

open(log_file, "w") do logio
    redirect_stdout(logio) do
        redirect_stderr(logio) do
            println("Loaded system with $(length(var_names)) variables and $(length(polys)) equations over GF($p)")

            # F5 computation
            @time sig_basis = sig_groebner_basis(polys; info_level=1)
            G = [f[2] for f in sig_basis]
            println("\nComputed Groebner basis (F5):")
            for g in G
                println(g)
            end

            # Verify correctness
            I = Ideal(polys)
            for f in getfield(I, :gens)
                r = normal_form(f, I)
                @assert r == 0
            end
            println("\nThe computed Groebner basis is correct!")

            # Save Groebner basis to result file
            open(result_file, "w") do fio
                println(fio, "# Groebner basis (F5) computed for $(filename)")
                println(fio, "# Variables: ", join(var_names, ", "))
                println(fio, "# Field characteristic: $p")
                println(fio, "# Number of input equations: $(length(polys))")
                println(fio, "# --- Groebner basis ---")
                for g in G
                    println(fio, g)
                end
            end
        end
    end
end

println("Computation done. Log: $log_file   Result: $result_file")
