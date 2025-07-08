using AlgebraicSolving
import Base.Filesystem: mkpath
using Dates

mkpath("logs")
mkpath("results")

# Read and parse the input file
filename = "data/toy_test.in"
lines = readlines(filename)
var_names = [strip(v) for v in split(strip(lines[1]), ",")]
p = parse(Int, strip(lines[2]))

# Create field and polynomial ring
if p == 2
    K = GF(2)
else
    K = GF(p)
end

R, vars = polynomial_ring(K, var_names)

# Define variables in Main module for eval
for (i, v) in enumerate(var_names)
    @eval Main $(Symbol(v)) = vars[$i]
end

# Initialize the correct type for the polynomials vector
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

I = Ideal(polys)

input_id = splitext(basename(filename))[1]
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")

log_file = "logs/$(input_id)_F4_$timestamp.log"
result_file = "results/$(input_id)_F4_$timestamp.txt"

open(log_file, "w") do logio
    redirect_stdout(logio) do
        redirect_stderr(logio) do
            println("Loaded system with $(length(var_names)) variables and $(length(polys)) equations over GF($p)")
            @time G = groebner_basis(I; info_level=2)
            println("\nComputed Groebner basis (F4):")
            for g in G
                println(g)
            end
            # Verify basis correctness
            for f in getfield(I, :gens)
                r = normal_form(f, I)
                @assert r == 0
            end
            println("\nThe computed Groebner basis is correct!")
            # Save Groebner basis to result file
            open(result_file, "w") do fio
                println(fio, "# Groebner basis (F4) computed for $(filename)")
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
