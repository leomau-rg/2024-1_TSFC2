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

using Test

@testset "Creación de Duales" begin
    @test typeof(Dual(1.0, 2.0)) == Dual{Float64}
    @test typeof(Dual(1.0, big(2.0))) == Dual{BigFloat}
    @test typeof(Dual(1, 2)) == Dual{Float64}
    @test typeof(Dual(big(1), 2)) == Dual{BigFloat}
    @test typeof(Dual(1, 2//1)) == Dual{Rational{Int}}
    @test typeof(Dual(1//1, 2//1)) == Dual{Rational{Int}}
    @test typeof(Dual(1//1, big(2)//1)) == Dual{Rational{BigInt}}
    @test typeof(Dual(1.0, big(2)//1)) == Dual{BigFloat}
    @test typeof(Dual(Int128(1))) == Dual{Float64}

    @test Dual(1.0, 2.0) == Dual{Float64}(1.0, 2.0)
    @test Dual(1.0, big(2.0)) == Dual{BigFloat}(1.0, 2.0)
    @test Dual(1, 2) == Dual{Float64}(1, 2)
    @test Dual(1.0, 2//1) == Dual{Rational{Int}}(1, 2)
    @test Dual(1//1, 2//1) == Dual{Rational{Int}}(1, 2)
    @test Dual(1//1, big(2)//1) == Dual(big(1)//1, 2//1)

    @test Dual(1) == Dual(1.0) == Dual(1.0, 0.0)
    @test 1.0 == Dual(1)
    @test Dual(big(1)) == Dual(big(1.0))
    @test Dual(big(1)) == big(1.0)
    @test !(Dual(1,2)==1)
    @test Dual(1//1) == Dual(big(1), 0)

end

a = Dual(1)
b = dual(2)
c = Dual(0.5, 0.5)

@testset "Operaciones aritméticas con Duales" begin
    @test +a == a
    @test -b == Dual(-b.fun, -b.der)
    @test a + b == Dual(3, 1)
    @test a - b == Dual(-1, -1)
    @test 2.5 * a == Dual(2.5)
    @test c == dual(1) * 0.5
    @test a * c == c
    @test b * c == Dual(1.0, 1.5)
    @test a/a == 1.0
    @test inv(a) == a/1.0
    @test b/b == 1.0
    @test inv(b) == Dual(0.5, -1/4)
    @test c/b == Dual(1/4, 1/8)
    @test b/c == inv(Dual(1/4, 1/8))
end

@testset "Funciones con Duales" begin
    @test a^2 == a
    @test b^2 == Dual(4, 4)
    @test sqrt(Dual(4, 4)) == Dual(4,4)^0.5 == b
    @test c^2 == Dual(1/4, 1/2)

    # Aquí usamos ≈ (en lugar de ==) porque puede variar de implementación a implementación
    @test fun(b^2.5) ≈ fun(Dual(sqrt(fun(b)^5), 2.5*fun(b)^1.5))
    @test der(b^2.5) ≈ der(Dual(sqrt(fun(b)^5), 2.5*fun(b)^1.5))
    b1 = b^2.5
    b2 = sqrt(b^5)
    @test fun(b1) ≈ fun(b2) && der(b1) ≈ der(b2)

    let
        f(x) = (3x^2-8x+5)/(7x^3-1)
        u = dual(2//1)
        fu = f(u)
        @test fun(fu) == 1//55
        @test der(fu) == 136//55^2
    end

    u = dual(0.5)
    @test exp(u) == Dual(exp(fun(u)), exp(fun(u))*der(u)) == Dual(1.6487212707001282, 1.6487212707001282)
    @test log(u) == Dual(log(fun(u)), der(u)/fun(u)) == Dual(- 0.6931471805599453, 2.0)
    @test sin(u) == Dual(sin(fun(u)), cos(fun(u))*der(u)) == Dual(0.479425538604203, 0.8775825618903728)
    @test cos(u) == Dual(cos(fun(u)), -sin(fun(u))*der(u)) == Dual(0.8775825618903728, -0.479425538604203)
    fu = tan(u)
    @test fu == Dual(tan(fun(u)), (sec(fun(u)))^2*der(u))
    @test der(fu) ≈ 1.2984464104095248
    @test asin(u) == Dual(asin(fun(u)), der(u)/sqrt(1-(fun(u))^2)) == Dual(0.5235987755982989, 1.1547005383792517)
    fu = acos(u)
    @test fu == Dual(acos(fun(u)), -der(u)/sqrt(1-(fun(u))^2))
    @test der(fu) ≈ - 1.1547005383792517
    @test atan(u) == Dual(atan(fun(u)), der(u)/(1+(fun(u))^2)) == Dual(0.4636476090008061, 0.8)
    @test atanh(u) == Dual(atanh(fun(u)), der(u)/(1-(fun(u))^2)) == Dual(0.5493061443340549, 1.3333333333333333)

    u = dual(2.0)
    @test sinh(u) == Dual(sinh(fun(u)), cosh(fun(u))*der(u)) == Dual(3.6268604078470186, 3.7621956910836314)
    @test cosh(u) == Dual(cosh(fun(u)), sinh(fun(u))*der(u)) == Dual(3.7621956910836314, 3.6268604078470186)
    fu = tanh(u)
    #@test fu == Dual(tanh(fun(u)), (sech(fun(u)))^2*der(u))
    @test fun(fu) ≈ fun(Dual(tanh(fun(u)), (sech(fun(u)))^2*der(u))) && der(fu) ≈ der(Dual(tanh(fun(u)), (sech(fun(u)))^2*der(u)))
    @test der(fu) ≈ 0.07065082485316443
    @test asinh(u) == Dual(asinh(fun(u)), der(u)/sqrt((fun(u))^2+1)) == Dual(1.4436354751788103, 0.4472135954999579)
    fu = acosh(u)
    @test fu == Dual(acosh(fun(u)), der(u)/sqrt((fun(u))^2-1))
    @test der(fu) ≈ 0.5773502691896258

    u = Dual(1.0, 2.0)
    @test exp(u) == Dual(exp(fun(u)), exp(fun(u))*der(u)) == Dual(2.718281828459045, 5.43656365691809)
    @test log(u) == Dual(log(fun(u)), der(u)/fun(u)) == Dual(0.0, 2.0)
    @test sin(u) == Dual(sin(fun(u)), cos(fun(u))*der(u)) == Dual(0.8414709848078965, 1.0806046117362795)
    @test cos(u) == Dual(cos(fun(u)), -sin(fun(u))*der(u)) == Dual(0.5403023058681398, - 1.682941969615793)
    @test tan(u) == Dual(tan(fun(u)), (sec(fun(u)))^2*der(u))
    # @test der(tan(u)) ≈ 6.85103764162952
    @test sinh(u) == Dual(sinh(fun(u)), cosh(fun(u))*der(u)) == Dual(1.1752011936438014, 3.0861612696304874)
    @test cosh(u) == Dual(cosh(fun(u)), sinh(fun(u))*der(u)) == Dual(1.5430806348152437, 2.3504023872876028)
    #@test tanh(u) == Dual(tanh(fun(u)), (sech(fun(u)))^2*der(u)) == Dual(0.7615941559557649, 0.8399486832280523)
    @test fun(tanh(u)) ≈ fun(Dual(tanh(fun(u)), (sech(fun(u)))^2*der(u))) ≈ fun(Dual(0.7615941559557649, 0.8399486832280523)) && der(tanh(u)) ≈ der(Dual(tanh(fun(u)), (sech(fun(u)))^2*der(u))) ≈ der(Dual(0.7615941559557649, 0.8399486832280523))
    @test asinh(u) == Dual(asinh(fun(u)), der(u)/sqrt((fun(u))^2+1)) == Dual(0.881373587019543, 1.414213562373095)

    u = Dual(0.0, 2.0)
    @test asin(u) == Dual(asin(fun(u)), der(u)/sqrt(1-(fun(u))^2)) == Dual(0.0, 2.0)
    @test acos(u) == Dual(acos(fun(u)), -der(u)/sqrt(1-(fun(u))^2)) == Dual(1.5707963267948966, -2.0)
    @test atan(u) == Dual(atan(fun(u)), der(u)/(1+(fun(u))^2)) == Dual(0.0, 2.0)
    @test atanh(u) == Dual(atanh(fun(u)), der(u)/(1-(fun(u))^2)) == Dual(0.0, 2.0)

    u = Dual(3.0, 2.0)
    @test acosh(u) == Dual(acosh(fun(u)), der(u)/sqrt((fun(u))^2-1)) == Dual(1.762747174039086, 0.7071067811865475)
end

"""
    paso_newton(f, x0)
Implementa un iterado del método de Newton para f(x) usando x0 como
punto de inicio. Devuelve el siguiente iterado.
"""
function paso_newton(f, x0)
    dualf = f(dual(x0))
    return x0 - fun(dualf) / der(dualf)
end

"""
    newton(f, x0, n)

Implementa n pasos del método de Newton para f(x) usando x0 como
punto de inicio. Devuelve el último iterado.
"""
function newton(f, x0, n)
    x1 = x0
    for i = 1:n
        x1 = paso_newton(f, x1)
    end
    return x1
end

@testset "Método de Newton en 1d" begin
    @test newton(x->x^2-2, 10.0, 20) == sqrt(2)
    @test newton(x->x^2-2, -2.0, 19) == -sqrt(2)

    let
        f(x) = 2*cos(x) - x + x^2/10
        x1 = newton(f, 10.0, 20)
        @test @show(abs(f(x1))) < 1.e-14

        x2 = newton(f, 1.0, 20)
        @test @show(abs(f(x2))) < 1.e-14
    end
end




