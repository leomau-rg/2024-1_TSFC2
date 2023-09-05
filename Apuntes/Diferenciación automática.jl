mutable struct Dual
    funcion::Float64
    derivada::Float64
end

# Definición de la suma:
import Base.:+
Base.:+(x::Dual,y::Dual) = Dual(x.funcion+y.funcion,x.derivada+y.derivada)

Dual(5.0,2.0) + Dual(1.0,1.0)


