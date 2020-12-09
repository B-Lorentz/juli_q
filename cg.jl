include("q_algebra.jl")
#J_d_coeff
#J_u_coeff
using Memoize
symRat(x::Rational) = string(" S(",x.num, ")/S(", x.den, ") ")
ket2arg(k::EigenKet) = string( symRat(k.j), ", ", symRat(k.m))

function pyCG(cg, a, b, c, io, logstr)

    println(io, "checky(", round(real(AlgebraicNumber(cg).apprx), digits=9),
     ", CG( ",ket2arg(a), ", ",
            ket2arg(b), ", ",
            ket2arg(c), ").doit().evalf(),'", logstr ,"')")
end
function pyprep(io::IO)
    println(io, "from sympy import *")
    println(io, "from sympy.physics.quantum.spin import CG")
    println(io, "import numpy as np")
    println(io, "def checky(a, b, log):")
    println(io, "   if not (np.allclose(a, float(b)) ):")
    println(io, "       print(f'{a:.6f}', f'{float(b):.6f}', log)")
end
@memoize function CG(a::EigenKet, b::EigenKet, c::EigenKet, io::IO)
    J = c.j
    M = c.m
    cg = 0
    logstr = ""
    if ((a.m + b.m) != M) | (J>a.j + b.j) | (J<abs(a.j - b.j))
        cg = 0
        logstr = string(logstr,"0 brach")
    elseif (J == M)
        m1s, m2s, cs = JJlist(J, a.j, b.j)
        cg = cs[m1s .== a.m][1]
            logstr = string(logstr,"JJ brach")
    else
        _M = M+1
        kl = jK(J, _M, 0)
        c_l = J_d_coeff(kl)

        sum = 0
        try
            k1 = jK(a.j, a.m+1, 0)
            c_r1 = J_d_coeff(k1)
            c1 = CG(k1, b, kl)
            sum += c_r1*c1
        catch e
            if isa(e, DomainError)
                logstr = string(logstr, " de ", a.j," ", a.m+1)
            else
                println(typeof(c_r1), typeof(c1))
                throw(e)
           end
        end
        try
            k2 = jK(b.j, b.m+1, 0)
            c_r2 = J_d_coeff(k2)
            sum += c_r2*CG(a, k2, kl, io)
        catch e
            if isa(e, DomainError)
                logstr = string(logstr, " de ", a.j," ", a.m+1)
            elseif isa(e, MethodError)
                throw(e)
            end
       end
        cg = sum/c_l
        logstr = string(logstr, " main branch")

    end
    pyCG(cg, a, b, c, io, logstr)
    cg
end

@memoize function JJlist(J::Rational, j1::Rational, j2::Rational)    #
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

function randomkets(maxj::Number)

    j1=rand(0:(2*maxj))//2
    m1=rand(-j1:j1)
    j2=rand(0:(2*maxj))//2
    m2=rand(-j2:j2)
    M=m1+m2
    J = rand(abs(M):(j1+j2))
    a=jK(j1,m1,0)
    b=jK(j2,m2,0)
    c=jK(J,M,0)
    return (a,b,c)
end
function pythoncheck(N::Integer, fname::String)
    io = open(fname, "w")
    pyprep(io)

    for _ in 0:N
        println(io, "print('____________________________________')")
        kets = randomkets(5)
        cg = CG(kets..., io)
        #pyCG(cg, kets..., io)
    end
    close(io)
end
pythoncheck(30,"check.py")
#println(jK(5//2, 5//2, 0))
#io = open( "check.py", "w")
#pyprep(io)
#println(CG( jK(9//2 , 3//2, 0) , jK(7//2 ,  -5//2, 0) ,  jK(3//1 ,  -1//1, 0), io) )
#close(io)
#println(randomkets(5))
