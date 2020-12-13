include("q_algebra.jl")
include("cg.jl")
using Memoize
@memoize function JJdecomp(
    j1::Rational,
    j2::Rational,
    J::Rational,
    hi1 = 0::Integer,
    hi2 = 1::Integer,
)
    m1s, m2s, coeffs = JJlist(J, j1, j2)
    res = ZeroKet()
    for m1m2coe in zip(m1s, m2s, coeffs)
        m1, m2, coe = m1m2coe
        res += coe * jK(j1, m1, hi1) ⊗ jK(j2, m2, hi2)
    end
    res
end

struct PrimitiveBra
    ket::PrimitiveKet
end

dag(a::PrimitiveKet) = PrimitiveBra(a)

function Base.:*(b::PrimitiveBra, k::PrimitiveKet)
    @assert(sameH(b.ket, k), "Can not take inner product between different Hilbert spaces")
    if b.ket == k
        return 1
    else
        return 0
    end
end
Base.:*(b::PrimitiveBra, k::LinKet) = sum(k.coeffs .* ([b] .* k.base_kets))

@memoize function decompose(
    j1::Rational,
    j2::Rational,
    J::Rational,
    M::Rational,
    hi1 = 0::Integer,
    hi2 = 1::Integer,
)
    if J == M
        return JJdecomp(j1, j2, J, hi1, hi2)
    else
        uket = jK(J, M + 1, 0)
        left = apply_op(uket, _J_d)
        dec = decompose(j1, j2, J, M + 1, hi1, hi2)
        right = apply_op(dec, _J_d)

        return (1 / (dag(jK(J, M, 0)) * left)) * right
    end
end

function check_decomposition(j1::Rational, j2::Rational, J::Rational, M::Rational, fname)
    io = open(fname, "w")
    pyprep(io, true)
    deco = decompose(j1, j2, J, M)
    K = jK(J, M, 0)
    for cket in zip(deco.coeffs, deco.base_kets)
        c, ket = cket
        pyCG(c, ket.parts[1], ket.parts[2], K, io, " ")
    end
    println(io, "print('All tests passed: ' , allok)")
    close(io)
end

full_prob(k::LinKet) = sum(alg2f.(k.coeffs) .^ 2)

function three_decompose(js, J::Rational, M::Rational)
    @assert(J in all_Js(js), string("J = ", J, " can not be decomposed to js: ", js))
    jrems = (abs(J - js[end])):(J+js[end])
    results = []

    for jrem in reverse(jrems)
        result = ZeroKet()
        bin = decompose(jrem, js[end], J, M, 0, length(js) - 1)
        #    println(full_prob(bin))
        for c_k in zip(bin.coeffs, bin.base_kets)
            c, k = c_k
            try
                #    deco = full_decompose()
                deco = decompose(js[1], js[2], k.parts[1].j, k.parts[1].m)

                result = result + c * (deco ⊗ k.parts[2])
            catch err
                if ~isa(err, AssertionError)
                    throw(err)
                end
            end
        end
        #    println()
        if ~(result == ZeroKet())
            push!(results, (result, jrem))
        end
    end
    results
end

function all_Js(js)
    if length(js) == 2
        return abs(js[1] - js[2]):(js[1]+js[2])
    else
        Jpos = all_Js(js[1:end-1])
        allJ = []
        for J in Jpos
            allJ = vcat(allJ, (abs(J - js[end]):(J+js[end])))
        end
        return sort(unique(allJ))
    end
end

function great_check(js)
    Js = all_Js(js)
    decons = []
    fulls = []
    for J in Js
        for M = -J:J
            push!()
            push!(decons, full_decompose(js, J, M))
        end
    end
    decons
end

function great_check2(js)
    Js = all_Js(js)
    decons = []
    fulls = []
    J = max(Js)
    for M = -J:J
        push!()
        push!(decons, full_decompose(js, J, M))
    end
    decons
end

function forward(js, ms)
    j_currents = [js[1]]
    dec_currents = [1]
    m_current = ms[1]
    histories = [[]]
    i = 0
    ket = jK(js[1], ms[1], 0)
    for jm in zip(js[2:end], ms[2:end])
        j, m = jm
        ket = ket ⊗ jK(j, m, i+1)
        j_nexts = []
        dec_nexts = []
        next_histories = []
        for thread in zip(j_currents, dec_currents, histories)
            j_, d_, his = thread
            jsums = max(abs(j_ - j), abs(m_current + m)):(j_+j)
            j_nexts = vcat(j_nexts, jsums)
            newket = jK(j, m, i+1)
            coeffs = (J -> CG(jK(j_, m_current, 0), newket, jK(J, m_current + m, 2)) ).(jsums)

            dec_nexts = vcat(dec_nexts, coeffs*d_)#(c-> c*(d_ ⊗ newket)).(coeffs))
            next_histories = vcat(next_histories, (jx->vcat(his, jx)).(jsums) )
        end
        j_currents = j_nexts
        dec_currents = dec_nexts
        m_current = m_current + m
        histories = next_histories
        i += 1
    end
    return ket, histories, dec_currents
end

function pp_partials(js)
    s = "\\left[ "
    i = 2
    for j in js[1:end-1]
        s = string(s, "j_{1...", i, "} = ", lpp(j), ", ")
        i += 1
    end
    s = string(s," \\right]")
    s
end
pp_coup(io, js, c, M) = println(io, lpp(c), " \\left| ", lpp(js[end]), ", ",lpp(M), ", ",pp_partials(js),"\\right\\rangle")
function print_decomposition(io, ket, subjs, coeffs)
    M = sum((k->k.m).(ket.parts))
    print(io, "\$\$ ", lpp(ket), " = ")
    pp_coup(io, subjs[1], coeffs[1], M)
    for jc = zip(subjs[2:end], coeffs[2:end])
        js, c = jc
        print( " + ")
        pp_coup(io, js, c, M)
    end
    println("\$\$")
end

function py_dectest(io::IO, ket::ProductKet)
    println(io, "import sympy")
    println(io, "from sympy.physics.quantum import TensorProduct")
    println(io, "from sympy.physics.quantum.spin import JzKet, couple")
    print(io, "couple(TensorProduct(")
    for part in ket.parts
        print(io, "JzKet(", ket2arg(part), "), ")
    end
    print(io, "))")
end

ket, subjs, coeffs = forward([1, 2, 1//2, 2], [1, 1, 1//2, -1])
println(ket)
println(subjs)
println(coeffs)
println(sum(coeffs .^ 2))
print_decomposition(stdout, ket, subjs, coeffs)
py_dectest(stdout, ket)
#println(great_check([1//1, 1//1,1//1]))
#res = three_decompose( [1//1, 1//1, 1//1], 1//1, 1//1 )
#println(res)
#println(full_prob(res))
#println(decompose(2//1, 1//1, 2//1, 2//1))
#println(JJdecomp(1//1, 1//1, 0//1))
#println(all_Js([3, 13, 5]))
#great_check([1//1,1//1,1//1])
#check_decomposition(8 + 5//2, 3//2, 8//1, 8//1, "de_check.py")
#println(apply_op(jK(5, 4, 0)⊗jK(3//2, -1//2, 1), _J_d))
#println(apply_op(jK(5, 5, 0)⊗jK(3//2, -3//2, 1), _J_d))
