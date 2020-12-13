include("ratp.jl")
include("../src/cg_unlimited.jl")
include("../src/q_algebra.jl")
js = []
J=ratparse(ARGS[2])
M=ratparse(ARGS[3])
for i in 1:((length(ARGS)-3))
    push!(js, ratparse(ARGS[i + 3]))

end
println("J = ", J)
println("M = ", M)
println("js=", js)
function all_ms(js)
    mss = (j-> [j]).(-js[1]:js[1])

    for j in js[2:end]
        nmss = []
        for ms in mss
            li = (jx->vcat(ms, jx)).(-j:j)

            nmss = vcat(nmss, li)
        end
        mss = nmss
    end
    mss
end

pp_coup1(io, ket, js, c, M) = println(io, lpp(c), "\\left( ", lpp(ket), " \\right)_{",pp_partials(js),"}")

function get_decomposition(io, J, M, js)
    mss = all_ms(js)
    probs = 0

    recs = []
    for ms in mss
        if sum(ms) == M
            ket, subjs, coeffs = forward(js, ms)
            for cp in zip( coeffs, subjs)
                c, p = cp
                if p[end] == J
                    push!(recs, (ket, p, c))
                    probs = probs + c^2
                end
            end
        end
    end
    print(io, "\$\$ ", lpp(jK(J, M, length(js)+1)), " = ")
    for rec in recs
        ket, p, c = rec
        print(io, " + ")
        pp_coup1(io, ket, p, c/sqrt(probs), M)
    end
    println(io, "\$\$")
    probs
end
io = open(ARGS[1], "w")
println(alg2f(AlgebraicNumber(get_decomposition(io, J, M, js))))
close(io)
