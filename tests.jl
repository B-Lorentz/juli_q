include("q_algebra.jl")
a_ = jK(3, 3, 2)
b_ = jK(3, -1, 2)
c_ = jK(1, 0, 3)
d_ = jK(3/2,1/2,3)
#println(a_ + b_)
#println((c_+d_)⊗(a_+b_))
#io = open("myfile.txt", "w")
#println(io, apply_op(a_, _J_d))
#flpp(io, convert(Array{Ket},[ apply_op(a_⊗c_, _J_d)]))
#println(io)
#close(io)
#-2 , 0*x , 3*x**2
function check_eq(a::Ket, b::Ket)
    println(a)
    println(b)
    @assert(typeof(a+(-1)*b) == ZeroKet)
    println()
end
check_eq(2*a_,a_+a_)
check_eq((1//4 + 1//5)*b_,(1//4)*b_ +(1//5)*b_)
check_eq(-3*(a_⊗c_), a_⊗(-3*c_))
check_eq(apply_op(jK(1, -1, 0), _J_d), ZeroKet())

# TODO: Ez nem fut le, de nem az assert, hanem valami üreslista-kezelés miatt az összedásban
#check_eq(apply_op(jK(1/2, 1/2, 0)⊗jK(1, 0, 1), _J_d),
#jK(1/2, -1/2, 0)⊗jK(1, 0, 1) + sqrt(AlgebraicNumber(2))*(jK(1/2, 1/2, 0)⊗jK(1,-1, 1)))
check_eq(a_ + b_, a_ + b_)
