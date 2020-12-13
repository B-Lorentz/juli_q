include("ratp.jl")
include("../src/cg_unlimited.jl")
include("../src/q_algebra.jl")
J=ratparse(ARGS[1])
M=ratparse(ARGS[2])
j1=ratparse(ARGS[3])
j2=ratparse(ARGS[4])

check_decomposition(j1, j2, J, M, ARGS[5])
