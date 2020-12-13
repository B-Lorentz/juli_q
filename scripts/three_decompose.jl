include("ratp.jl")
include("../src/cg_unlimited.jl")
include("../src/q_algebra.jl")
J=ratparse(ARGS[1])
M=ratparse(ARGS[2])
j1=ratparse(ARGS[3])
j2=ratparse(ARGS[4])
j3=ratparse(ARGS[5])

list = three_decompose([j1, j2, j3], J, M)

io = open(ARGS[6], "w")
for tup in list
    println(io, string("\$\$ \\left[", lpp(jK(J, M, 2)), " = ",lpp(tup[1]),  "\\right] \\ ha \\ (J_{12} =", lpp(tup[2]), " )\$\$"))
end
close(io)
