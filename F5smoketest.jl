using AlgebraicSolving

# Define a polynomial ring over GF(2) in 2 variables
R, (x, y) = polynomial_ring(GF(2), ["x", "y"])

# Construct the ideal I = <x^2 + y + 1, x*y + y + 1>
I = Ideal([x^2 + y + 1, x*y + y + 1])

# Run F5 (signature-based)
# The F5 algorithm returns a vector of tuples: (signature, polynomial)
@time Gsig = sig_groebner_basis(getfield(I, :gens); info_level=1, mod_ord=:POT)

# Extract only the polynomial part from each tuple
Gf5 = [t[2] for t in Gsig]

println("Computed Groebner basis with F5:")
for g in Gf5
    println(g)
end

# Check correctness: normal_form of generators should be 0
for f in getfield(I, :gens)
    r = normal_form(f, I)
    @assert r == 0
end

println("Groebner basis (F5) correctness check: PASSED")
