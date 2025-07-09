using AlgebraicSolving
import Base.Filesystem: mkpath
using Dates

# --- Parse command-line arguments
function parse_args()
    if length(ARGS) < 1
        println("Usage: julia solve_F4_from_file.jl <inputfile> [nthreads]")
        exit(1)
    end
    filename = ARGS[1]
    nthreads = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
    return filename, nthreads
end

filename, nthreads = parse_args()

mkpath("logs")
mkpath("results")

lines = readlines(filename)
var_names = [strip(v) for v in split(strip(lines[1]), ",")]

# --- Field parsing logic supporting both p and p^n formats ---
field_spec = strip(lines[2])
if occursin("^", field_spec)
    base, ext = split(field_spec, "^")
    p = parse(Int, strip(base))
    n = parse(Int, strip(ext))
    K = GF(p, n, "a")  # "a" as primitive element
    field_desc = "GF($p^$n)"
else
    p = parse(Int, field_spec)
    n = 1
    K = GF(p)
    field_desc = "GF($p)"
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
            println("Loaded system with $(length(var_names)) variables and $(length(polys)) equations over $field_desc")
            println("Using $nthreads threads for computation.")
            @time G = groebner_basis(I; info_level=2, nr_thrds=nthreads)
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
                println(fio, "# Field: $field_desc")
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
