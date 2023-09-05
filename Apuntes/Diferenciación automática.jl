mutable struct Dual{T <: Real}
    fun :: T
    der :: T
end

function Dual(x::T,y::S) where {T<:Real, S<:Real}
    return Dual(promote(x,y)...)
end

function fun(x::Dual)
    return x.fun
end

function der(x::Dual)
    return x.der
end

function Dual(c::Real)
    return Dual(c,0)
end

function dual(x)
    return Dual(x,1)
end

import Base.:+, Base.:-, Base.:*, Base.:/

function +(x::Dual,y::Dual)
    Dual(x.fun+y.fun,x.der+y.der)
end

function -(x::Dual,y::Dual)
    Dual(x.fun-y.fun,x.der-y.der)
end

function *(x::Dual,y::Dual)
    Dual(x.fun*y.fun,x.fun*y.der+x.der*y.fun)
end

function /(x::Dual,y::Dual)
    Dual(x.fun/y.fun,(y.fun*x.der-x.fun*y.der)/y.fun^2)
end
