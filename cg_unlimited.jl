include("q_algebra.jl")
include("cg.jl")
using Memoize
@memoize function JJdecomp(j1::Rational, j2::Rational, J::Rational)
    m1s, m2s, coeffs = JJlist(J, j1, j2)
    res = ZeroKet()
    for m1m2coe in zip(m1s, m2s, coeffs)
        m1, m2, coe =  m1m2coe
        res += coe*jK(j1, m1, 0)⊗jK(j2, m2, 1)
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

@memoize function decompose(j1::Rational, j2::Rational, J::Rational, M::Rational)
    if J == M
        return JJdecomp(j1, j2, J)
    else
        uket = jK(J, M+1, 0)
        left = apply_op(uket, _J_d)
        dec = decompose(j1, j2, J, M+1)
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

function full_decompose(js, J::Rational, M::Rational)
    if length(js) == 2
        return decompose(js[1], js[2], J, M)
    else
        jrems = (J-js[end]):(J+js[end])
        for jrem in jrems
            println(jrem)
            deco = full_decompose(js[1:end-1], jrem, M)
            println(deco)
        end
        return 420
    end
end
println(full_decompose([1//1, 2//1, 1//1], 3//1, 2//1 ))
#check_decomposition(8 + 5//2, 3//2, 8//1, 8//1, "de_check.py")
#println(apply_op(jK(5, 4, 0)⊗jK(3//2, -1//2, 1), _J_d))
#println(apply_op(jK(5, 5, 0)⊗jK(3//2, -3//2, 1), _J_d))
