module DifAutom

    ### 1

    struct Dual{T <: Real}
        fun :: T
        der :: T
    end

    function Dual(x::T,y::S) where {T<:Real, S<:Real}
        return Dual(promote(x,y)...)
    end

    function Dual(x::T,y::S) where {T<:Integer, S<:Integer}
        x, y, _ = promote(x,y,1.0)
        return Dual(x,y)
    end

    function fun(x::Dual)
        return x.fun
    end

    function der(x::Dual)
        return x.der
    end

    ### 2

    function Dual(c::Real)
    return Dual(c,0)
    end

    function dual(x::Real)
        return Dual(x,1)
    end

    import Base.:+, Base.:-, Base.:*, Base.:/, Base.==, Base.:inv

    function +(x::Dual,y::Dual)
        return Dual(x.fun+y.fun,x.der+y.der)
    end

    function +(x::Dual)
        return x
    end

    function -(x::Dual,y::Dual)
        return Dual(x.fun-y.fun,x.der-y.der)
    end

    function -(x::Dual)
        return Dual(-x.fun,-x.der)
    end

    function *(x::Dual,y::Dual)
        return Dual(x.fun*y.fun,x.fun*y.der+x.der*y.fun)
    end

    function /(x::Dual,y::Dual)
        return Dual(x.fun/y.fun,(y.fun*x.der-x.fun*y.der)/y.fun^2)
    end

    function ==(x::Dual,y::Dual)
        return (x.fun == y.fun) && (x.der == y.der)
    end

    function +(x::Real,y::Dual)
        return Dual(x+y.fun,y.der)
    end

    function -(x::Real,y::Dual)
        return Dual(x-y.fun,-y.der)
    end

    function *(x::Real,y::Dual)
        return Dual(x*y.fun,x*y.der)
    end

    function /(x::Real,y::Dual)
        return Dual(x/y.fun,(-x*y.der)/y.fun^2)
    end

    function ==(x::Real,y::Dual)
        return Dual(x) == y
    end

    function +(x::Dual,y::Real)
        return Dual(x.fun+y,x.der)
    end

    function -(x::Dual,y::Real)
        return Dual(x.fun-y,x.der)
    end

    function *(x::Dual,y::Real)
        return y*x
    end

    function /(x::Dual,y::Real)
        return Dual(x.fun/y,(y*x.der)/y^2)
    end

    function ==(x::Dual,y::Real)
        return x == Dual(y)
    end

    function inv(x::Dual)
       return Dual(1/x.fun,-x.der/x.fun^2) 
    end

    ### 3

    import Base.Math.:^
    function ^(x::Dual,y::Real)
        return Dual(x.fun^y,y*x.fun^(y-1)*x.der)
    end

    import Base.Math.:sqrt
    function sqrt(x::Dual)
        return Dual(sqrt(x.fun),x.der*1/(2*sqrt(x.fun)))
    end

    import Base.Math.:sin
    function sin(x::Dual)
        return Dual(sin(x.fun),x.der*cos(x.fun))
    end

    import Base.Math.:cos
    function cos(x::Dual)
        return Dual(cos(x.fun),-x.der*sin(x.fun))
    end

    import Base.Math.:tan
    function tan(x::Dual)
        return Dual(tan(x.fun),x.der*sec(x.fun)^2)
    end

    import Base.Math.:asin
    function asin(x::Dual)
        return Dual(asin(x.fun),x.der*1/sqrt(1-x.fun^2))
    end

    import Base.Math.:acos
    function acos(x::Dual)
        return Dual(acos(x.fun),-x.der*1/sqrt(1-x.fun^2))
    end

    import Base.Math.:atan
    function atan(x::Dual)
        return Dual(atan(x.fun),x.der*1/(1+x.fun^2))
    end

    import Base.Math.:sinh
    function sinh(x::Dual)
        return Dual(sinh(x.fun),x.der*cosh(x.fun))
    end

    import Base.Math.:cosh
    function cosh(x::Dual)
        return Dual(cosh(x.fun),x.der*sinh(x.fun))
    end

    import Base.Math.:tanh
    function tanh(x::Dual)
        return Dual(tanh(x.fun),x.der*sech(x.fun)^2)
    end

    import Base.Math.:asinh
    function asinh(x::Dual)
        return Dual(asinh(x.fun),x.der/cosh(asinh(x.fun)))
    end

    import Base.Math.:acosh
    function acosh(x::Dual)
        return Dual(acosh(x.fun),x.der/sinh(acosh(x.fun)))
    end

    import Base.Math.:atanh
    function atanh(x::Dual)
        return Dual(atanh(x.fun),x.der*cosh(atanh(x.fun))^2)
    end

    import Base.Math.:exp
    function exp(x::Dual)
        return Dual(exp(x.fun),x.der*exp(x.fun))
    end

    import Base.Math.:log
    function log(x::Dual)
        return Dual(log(x.fun),x.der/x.fun)
    end

end
