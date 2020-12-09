include("q_algebra.jl")
#J_d_coeff
#J_u_coeff
function CG(a::EigenKet, b::EigenKet, c::EigenKet)
    J = c.j
    M = c.m
    if ((a.m + b.m) != M) | (J>a.j + b.j) | (J<abs(a.j - b.j))
        return 0
    elseif (J == M)
        m1s, m2s, cs = JJlist(J, a.j, b.j)
        return cs[m1s .== a.m][1]
    else

        _M = M+1

        kl = jK(J, _M, 0)
        c_l = J_d_coeff(kl)

        sum = 0


        try
            k1 = jK(a.j, a.m+1, 0)
            c_r1 = J_d_coeff(k1)
            sum += c_r1*CG(k1, b, kl)
        catch DomainError
        end
        try
            k2 = jK(b.j, b.m+1, 0)
            c_r2 = J_d_coeff(k2)
            sum += c_r2*CG(a, k2, kl)
        catch DomainError
        end
        return sum/c_l
    end
end

function JJlist(J::Rational, j1::Rational, j2::Rational)    #
    @assert((J <= j1 + j2) & ( J >= abs(j1 - j2)), string("J= ", J, " j1 =" ,j1, "j2 = ", j2, " are nonsense"))
    @assert(ishalf(j1), "j1 spin must be half-integer" )
    @assert(ishalf(j2), "j2 spin must be half-integer")
    @assert(ishalf(J), "J spin must be half-integer")
    @assert( (j1+j2).den == J.den , string(J," can't be decomposed into",j1," and ",j2))
    maxm1 = min(j1, J+j2)
    minm1 = max(-j1, J-j2)
    m1s = Array(minm1:(maxm1))
    m2s = J .-m1s

    coeffs = [AlgebraicNumber(1)]
    i = 1
    for m1m2 in zip(m1s[2:end], m2s[2:end])
        m1, m2 = m1m2
        #println(m1," ", m2)
        _m1 = m1
        _m2 = m2+1
        c1 = J_d_coeff(jK(j1, _m1, 0))
        c2 = J_d_coeff(jK(j2, _m2, 0))
        push!(coeffs, -coeffs[i]*c1/c2)
        i += 1
    end
    CSval = coeffs[m1s .== j1]

    (m1s, m2s, coeffs / sqrt(sum(coeffs .^2))*sign(CSval[1]))

end

function randomkets(maxj::Number)::(EigenKet, EigenKet, EigenKet)
    J=rand(0,2*maxj)//2
    diff=rand(1,10)
    j1=rand(0,2*maxj)//2
    println(j1)
    m1=(-j1:j1)[rand(1,2*j1+1)]
    j2=J+diff-j1
    m2=(-j2:j2)[rand(1,2*j2+1)]
    M=m1+m2
    a=jK(j1,m1,0)
    b=jK(j2,m2,0)
    c=jK(J,M,0)
    return (a,b,c)
end

function pythoncheck(N::Integer, fname::String)
    io = open(fname, "w")
    println(io, "from sympy import *")
    println(io, "from sympy.physics.quantum.spin import CG")
    println()
    for _ in 0:N
        kets = randomkets(10)
        cg = CG(kets...)
        println(io, cg)
    end
    close(io)
end
pythoncheck(5,"check.py")
