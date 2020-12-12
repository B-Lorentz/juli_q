include("q_algebra.jl")
include("cg.jl")
using Memoize
@memoize function JJdecomp(j1::Rational, j2::Rational, J::Rational, hi1=0::Integer, hi2=1::Integer)
    m1s, m2s, coeffs = JJlist(J, j1, j2)
    res = ZeroKet()
    for m1m2coe in zip(m1s, m2s, coeffs)
        m1, m2, coe =  m1m2coe
        res += coe*jK(j1, m1, hi1)⊗jK(j2, m2, hi2)
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
Base.:*(b::PrimitiveBra, k::LinKet) = sum(k.coeffs .* ( [b] .* k.base_kets))

@memoize function decompose(j1::Rational, j2::Rational, J::Rational, M::Rational, hi1=0::Integer, hi2=1::Integer)
    if J == M
        return JJdecomp(j1, j2, J, hi1, hi2)
    else
        uket = jK(J, M+1, 0)
        left = apply_op(uket, _J_d)
        dec = decompose(j1, j2, J, M+1, hi1, hi2)
        right = apply_op(dec, _J_d)

        return (1/(dag(jK(J, M, 0))*left))*right
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

full_prob(k::LinKet) = sum(alg2f.( k.coeffs) .^2)

function full_decompose(js, J::Rational, M::Rational)
    if length(js) == 2
        return decompose(js[1], js[2], J, M, 0, 1)
    else
        jrems = (abs(J-js[end])):(J+js[end])
        results = ZeroKet()
        sum = 0
        for jrem in reverse(jrems)
            result = ZeroKet()
            bin = decompose(jrem, js[end], J,M, 0, length(js)-1)
        #    println(full_prob(bin))
            for c_k = zip(bin.coeffs, bin.base_kets)
                c, k = c_k
                try
                    deco = full_decompose(js[1:end-1], k.parts[1].j, k.parts[1].m)

                    result = result + c*(deco⊗k.parts[2])
                catch err
                    if ~isa(err, AssertionError)
                        throw(err)
                    end
                end
            end
        #    println()
            results += result
            if result != ZeroKet()
                return result
            end
        end
    #    return (1/sqrt(AlgebraicNumber(sum)))*results
    end
end

function all_Js(js)
    if length(js) == 2
        return abs(js[1]-js[2]):(js[1] + js[2])
    else
        Jpos = all_Js(js[1:end-1])
        allJ =  []
        for J in Jpos
            allJ = vcat(allJ, (abs(J-js[end]):(J + js[end])))
        end
        return sort(unique(allJ))
    end
end

function great_check(js)
    Js = all_Js(js)
    decons = []
    fulls = []
    for J in Js
        for M in -J:J
            push!()
            push!(decons, full_decompose(js, J,M ))
        end
    end
    decons
end

function great_check2(js)
    Js = all_Js(js)
    decons = []
    fulls = []
    J=max(Js)
    for M in -J:J
        push!()
        push!(decons, full_decompose(js, J,M ))
    end
    decons
end

#println(great_check([1//1, 1//1,1//1]))
res = full_decompose( [1//1, 1//1, 1//1], 1//1, 1//1 )
println(res)
println(full_prob(res))
#println(decompose(2//1, 1//1, 2//1, 2//1))
#println(JJdecomp(1//1, 1//1, 0//1))
#println(all_Js([3, 13, 5]))
#great_check([1//1,1//1,1//1])
#check_decomposition(8 + 5//2, 3//2, 8//1, 8//1, "de_check.py")
#println(apply_op(jK(5, 4, 0)⊗jK(3//2, -1//2, 1), _J_d))
#println(apply_op(jK(5, 5, 0)⊗jK(3//2, -3//2, 1), _J_d))
