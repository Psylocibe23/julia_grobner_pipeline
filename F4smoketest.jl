using AlgebraicSolving

# Define a polynomial ring over GF(2) in 2 variables
R, (x, y) = polynomial_ring(GF(2), ["x", "y"])

# Construct an ideal I = <x^2 + y + 1, x*y + y + 1>
I = Ideal([x^2 + y + 1, x*y + y + 1])

println(typeof(I))
println(fieldnames(typeof(I)))
println(getfield(I, :gens))

# Compute Groebner basis using default (F4 variant)
@time G = groebner_basis(I; info_level=1)

println("Computed Groebner basis:")
for g in G
    println(g)
end

# Check correctness: normal_form of generators should be 0
for f in getfield(I, :gens)
    r = normal_form(f, I)
    @assert r == 0
end

println("Groebner basis correctness check: PASSED")
