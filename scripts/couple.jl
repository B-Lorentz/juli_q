include("ratp.jl")
include("../src/cg_unlimited.jl")
include("../src/q_algebra.jl")
js = []
ms = []
for i in 0:((length(ARGS)-1)÷2 -1)
    push!(js, ratparse(ARGS[2*i + 2]))
    push!(ms, ratparse(ARGS[2*i + 3]))
end
ket, subjs, coeffs = forward(js, ms)

println(alg2f(AlgebraicNumber(sum(coeffs .^ 2))))
io = open(ARGS[1], "w")
print_decomposition(io, ket, subjs, coeffs)
println(io, "___________________________________________________________________________")
py_dectest(io, ket)
close(io)
