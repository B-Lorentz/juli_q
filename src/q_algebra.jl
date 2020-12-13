#tenzorszorzás jel:⊗
using AlgebraicNumbers
struct EigenKet
    j::Rational
    m::Rational
    hil_ind::Integer        #melyik Hilbert-térben van
end
struct ZeroKet
       end
function ishalf(x::Rational)::Bool
    (x.den == 1)|(x.den==2)
end
function jK(j::Number, m::Number, h::Integer)
    j = convert(Rational, j)
    m = convert(Rational, m)
    if (j >= m) & (j >= -m) & ishalf(j) & ishalf(m) &(j.den == m.den)
        return EigenKet(j, m, h)
    else
        throw( DomainError((j, m), "Invalid quantum numbers") )
    end
end
function Base.:(==)(a::EigenKet, b::EigenKet)::Bool
    if ((a.j == b.j) & (a.m == b.m) & (a.hil_ind == b.hil_ind))
        return true
    else
        return false
    end
end

function Base.:>(a::EigenKet, b::EigenKet)::Bool        #definiáljuk, hogy hogyan rendezzük
    if (a.hil_ind > b.hil_ind)
        return true
    elseif (a.hil_ind < b.hil_ind)
        return false
    else
        if (a.j > b.j)
            return true
        elseif (a.j < b.j)
            return false
        else
            if (a.m > b.m)
                return true
            else
                return false
            end
        end
    end
end

function Base.:<(a::EigenKet, b::EigenKet)::Bool        #ez is a rendezéshez
    (~(a > b))&(~(a==b))
end
sameH(a::EigenKet, b::EigenKet) = a.hil_ind == b.hil_ind
diffH(a::EigenKet, b::EigenKet) = ~sameH(a, b)


struct ProductKet           #tenzorszorzata ket-eknek
    parts::Array{EigenKet}
end

Base.:(==)(a::ProductKet, b::ProductKet) = a.parts == b.parts
Base.:(>)(a::ProductKet, b::ProductKet) = a.parts > b.parts
Base.:(<)(a::ProductKet, b::ProductKet) = a.parts < b.parts

sameH(a::ProductKet, b::ProductKet) = all(map(x->sameH(x[1], x[2]), zip(a.parts, b.parts)))

diffH(a::ProductKet, b::ProductKet) = ~any(map(x->sameH(x[1], x[2]), zip(a.parts, b.parts)))
diffH(a::ProductKet, b::EigenKet) = ~any(map(x->sameH(x, b), a.parts))
diffH(b::EigenKet, a::ProductKet) = diffH(a, b)
function hfts(a)            #szép printelés, plusz jelez ha nem félegészek
    if a.den == 1
        return string(a.num)

    elseif a.den == 2
        return string(a.num, "/2")
    else
        throw(DomainError(a, "Not a half-integer"))
    end
end

Base.show(io::IO, a::EigenKet) = print(io, "|", hfts(a.j), " ; ", hfts(a.m), ">_" ,a.hil_ind)

PrimitiveKet = Union{EigenKet,ProductKet}


Base.isless(a::PrimitiveKet,b::PrimitiveKet) = a < b
Base.show(io::IO, a::ZeroKet) = print(io," 0")
function Base.show(io::IO , a::ProductKet)
    print(io,a.parts[1])
    for ket in a.parts[2:end]
        print(io," x ",ket)
    end
end

struct LinKet               #linkombja az EigenKet-eknek
    coeffs::Array{Number}
    base_kets::Array{PrimitiveKet}
end

NonZeroKet = Union{PrimitiveKet,LinKet}

alg2f(x::AlgebraicNumber) = convert(Float64,real(x.apprx))
function Base.show(io::IO, a::LinKet)
    print(io, " ", alg2f(a.coeffs[1]),"(", a.base_kets[1], ")")
    for cket in zip(a.coeffs[2:end], a.base_kets[2:end])
        print(io, " + ", alg2f(cket[1]), "(", cket[2], ") ")
    end
    println()
end

function order(a::LinKet)       #sorbarendezi a LinKet-et, a korábban definiált módon

   p = sortperm(a.base_kets)
   newcoeffs=a.coeffs[p]
   newbase_kets=a.base_kets[p]
   LinKet(newcoeffs,newbase_kets)
end

function Base.:*(x::Number, y::PrimitiveKet)
    if x == 0
        return ZeroKet()
    else
        LinKet([x],[y])
    end
end

function Base.:*(x::Number, y::LinKet)
    if x == 0
        return ZeroKet()
    else
        LinKet(x*y.coeffs,y.base_kets)
    end
end

Base.:*(x::Number, y::ZeroKet) = y

function Base.:+(a::PrimitiveKet, b::PrimitiveKet)
    @assert(sameH(a,b), "Not in same Hilbert space")
    if a == b
        return LinKet([2], [a])
    else
        return LinKet([1,1], sort([a, b]))
    end
end

function q_add(a::LinKet, b::PrimitiveKet,c::Number)
    @assert(sameH(a.base_kets[1],b), "Not in same Hilbert space")
    newcoeffs=[]
    newkets=[]
    isfound=false
    for ipair in zip(a.coeffs,a.base_kets)      #végignézi, hogy van-e már ugyanolyan Ket
        coeff,ket=ipair
        if (ket==b)
            coeff+=c             #ha van ugyanolyan, akkor hozzáad 1-et a coeff-hez
            isfound=true
        end
        if (coeff!=0)
            push!(newkets,ket)
            push!(newcoeffs,coeff)   #kigyűjti a coeff-eket
        end
    end
    #println(newcoeffs, newkets)
    if length(newkets)==0
        return ZeroKet()
    elseif (isfound) #esetleg ez a rész is lefut ha ZeroKet-et kapunk?
        return LinKet(newcoeffs, newkets)
    else
        ans= order(LinKet(vcat(newcoeffs,[c]),vcat(a.base_kets,[b])))
        return ans#LinKet(vcat(newcoeffs,[c]),vcat(newkets,[b]))
    end
end

function _q_add(a::LinKet, b::PrimitiveKet,c::Number)
    @assert(sameH(a.base_kets[1],b), "Not in same Hilbert space")
    if length(a.coeffs) == 1
        if a.base_kets[1] == b
            return (a.coeffs[1] + c)*b
        else
            return order(LinKet(vcat(a.coeffs,[c]),vcat(a.base_kets,[b])))
        end
    else
        if a.base_kets[1] == b
            return (a.coeffs[1] + c)*b
        end
    end
end

Base.:+(a::LinKet, b::PrimitiveKet) = q_add(a,b,1)

Base.:+( b::PrimitiveKet, a::LinKet) = a+b
Ket = Union{ZeroKet, NonZeroKet}
Base.:+(a::ZeroKet,b::NonZeroKet) = b
Base.:+(b::NonZeroKet,a::ZeroKet) = b
Base.:+(a::ZeroKet, b::ZeroKet) = a

function Base.:+(a::LinKet, b::LinKet)
   @assert(sameH(a.base_kets[1], b.base_kets[1]), "Not in same Hilbert space")
   ans=a
   for ipair in zip(b.coeffs,b.base_kets)
       #println("bfore:" , ans)
       ans=q_add(ans,ipair[2],ipair[1])
       #println("after:", ans)
   end
   return ans
end

lpp(x::Integer) = string(x)
function lpp(a::Integer, b::Integer)
    if b == 1
        return string(a)
    else
        return string("\\frac{",a,"}{",b,"}")
    end
end
lpp(a::Rational) = lpp(a.num, a.den)
function lpp(x::AlgebraicNumber)

    if imag(x.apprx)==0
        if length(x.coeff) == 2
            lpp(-x.coeff[1], x.coeff[2])
        elseif length(x.coeff) == 3
            if x.coeff[2] == 0
                sgn = ""
                if real(x.apprx)<0
                    sgn = "-"
                end
                return string(sgn,"\\sqrt{", lpp(abs(x.coeff[1]), abs(x.coeff[3])),"}")
            else
                print(x.coeff)
                return ""
            end
        else
            return string(x)
        end
    else
        return string(x)
    end
end

function lpp(a::EigenKet)
    string( "\\left|", lpp(a.j), ",", lpp(a.m), "\\right\\rangle_" ,a.hil_ind)
end
function lpp(vec::LinKet)
    ans=string(lpp(vec.coeffs[1]),"\\ ",lpp(vec.base_kets[1]))
    for ck in zip(vec.coeffs[2:end],vec.base_kets[2:end])
        coeff, ket = ck
        ans=string(ans," + ",lpp(coeff), "\\ ",lpp(ket))
    end
    return ans
end
lpp(x::ZeroKet)=string(0)

function lpp(x::ProductKet)
    ans=string(lpp(x.parts[1]))
    for ket in x.parts[2:end]
        ans = string(ans, " \\otimes ", lpp(ket))
    end
    return ans
end

function flpp(io::IOStream, kets::Array{Ket})
    for ket in kets
        println(io, "\$\$",lpp(ket),"\$\$")
    end
end
function ⊗(a::EigenKet,b::EigenKet)
    @assert(diffH(a, b), "Can't ⊗, in the same Hilbert space")
    ProductKet([a,b])
end

⊗(a::ZeroKet, b::Ket) = a
⊗(a::Ket, b::ZeroKet) = b


function ⊗(a::ProductKet, b::EigenKet)
    @assert(diffH(a, b), "Can't ⊗, partly in the same Hilbert space")
    ProductKet(vcat(a.parts, [b]))
end

function ⊗(a::EigenKet, b::ProductKet)
    @assert(diffH(a, b), "Can't ⊗, partly in the same Hilbert space")
    ProductKet(vcat([a], b.parts))
end

function ⊗(a::LinKet,b::PrimitiveKet)
    newbase=[]
    for ket in a.base_kets
        push!(newbase,ket⊗b)
    end
    LinKet(a.coeffs,newbase)
end

function ⊗(a::PrimitiveKet,b::LinKet)
    newbase=[]
    for ket in b.base_kets
        push!(newbase,a⊗ket)
    end
    LinKet(b.coeffs,newbase)
end

function ⊗(a::LinKet, b::LinKet)
    parts = []
    for cak in zip(a.coeffs, a.base_kets)
        ca, ck = cak
        parts = vcat(parts, ca*(ck⊗b))
    end
    sum(parts)
end

function J_d_coeff(x::EigenKet)
    j=AlgebraicNumber(x.j)
    m=AlgebraicNumber(x.m)
    sqrt(j*(j+1)-m*(m-1))
end

function J_u_coeff(x::EigenKet)
    j=AlgebraicNumber(x.j)
    m=AlgebraicNumber(x.m)
    sqrt(j*(j+1)-m*(m+1))
end

function _J_d(a::EigenKet)
    if -a.m == a.j
        return ZeroKet()
    else
        return J_d_coeff(a)*EigenKet(a.j, a.m-1, a.hil_ind)
    end
end
function _J_u(a::EigenKet)
    if a.m == a.j
        return ZeroKet()
    else
        return J_u_coeff(a)*EigenKet(a.j, a.m+1, a.hil_ind)
    end
end
apply_op(a::EigenKet, op) = op(a)
function apply_op(a::LinKet, op)
    result = ZeroKet()
    for ipair in zip(a.coeffs, a.base_kets)      #végignézi, hogy van-e már ugyanolyan Ket
        coeff,ket=ipair
        result=result+ coeff * apply_op(ket, op)
    end
    result
end
function apply_op(a::ProductKet, op)
    if length(a.parts) == 1
        return apply_op(a.parts[1], op)
    else
        return apply_op(a.parts[1], op)⊗ProductKet(a.parts[2:end]) + a.parts[1]⊗apply_op(ProductKet(a.parts[2:end]), op)
    end
end
