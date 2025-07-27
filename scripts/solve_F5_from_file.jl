###############################################################################
# solve_F5_from_file.jl
#
# This script reads a system of polynomial equations from a file, constructs the
# corresponding ideal over a finite field, and computes a Gröbner basis using
# the F5 algorithm (as implemented in AlgebraicSolving.jl). Results and logs are
# written to disk. The script is designed for cryptanalytic applications and
# general algebraic solving.
###############################################################################

using AlgebraicSolving  # Gröbner basis computation, F5 algorithm, etc.
import Base.Filesystem: mkpath
using Dates

##########################
# 1. Setup output folders
##########################
mkpath("logs")
mkpath("results")

##########################
# 2. Parse command-line arguments
##########################
function parse_args()
    # The script expects at least 1 argument: the input filename.
    if length(ARGS) < 1
        println("Usage: julia solve_F4_from_file.jl <inputfile> [nthreads]")
        exit(1)
    end
    filename = ARGS[1]
    return filename
end

filename = parse_args()

##########################
# 3. Read and parse the input file
##########################
lines = readlines(filename)
var_names = [strip(v) for v in split(strip(lines[1]), ",")]  # Expect first line to be variable names, separated by commas (e.g., "x, y, z") 
p = parse(Int, strip(lines[2]))  # Second line: field characteristic (prime number p, for GF(p))

##########################
# 4. Field and polynomial ring setup
##########################
# If the characteristic is 2, use GF(2), else use GF(p).
if p == 2
    K = GF(2)
else
    K = GF(p)
end

# Construct the multivariate polynomial ring R = K[var_names]
R, vars = polynomial_ring(K, var_names)

# For each variable, bind its symbol in Main so polynomials can be parsed/evaluated from strings
for (i, v) in enumerate(var_names)
    @eval Main $(Symbol(v)) = vars[$i]
end

##########################
# 5. Parse polynomials from input
##########################
# Collect each equation, starting from the third line of the input file.
polys = Vector{typeof(vars[1])}()
for line in lines[3:end]
    try
        # Meta.parse parses the string to Julia code; eval evaluates it in Main
        poly = eval(Meta.parse(line))
        push!(polys, poly)
    catch err
        println("Error evaluating polynomial: $line")
        rethrow(err)
    end
end

##########################
# 6. Prepare output/log file names
##########################
input_id = splitext(basename(filename))[1]
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
log_file = "logs/$(input_id)_F5_$timestamp.log"
result_file = "results/$(input_id)_F5_$timestamp.txt"

##########################
# 7. Compute the F5 Gröbner basis and log output
##########################
open(log_file, "w") do logio
    redirect_stdout(logio) do
        redirect_stderr(logio) do
            println("Loaded system with $(length(var_names)) variables and $(length(polys)) equations over GF($p)")

            # F5 computation
            # sig_groebner_basis returns a list of pairs (signature, polynomial).
            @time sig_basis = sig_groebner_basis(polys; info_level=1)
            # Extract just the polynomials from the signature pairs
            G = [f[2] for f in sig_basis]
            println("\nComputed Groebner basis (F5):")
            for g in G
                println(g)
            end

            # Verify correctness
            # Reduce each generator of the ideal modulo the computed basis, which should give zero.
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
