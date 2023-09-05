using Plots

f(x) = 3*x^3 - 2

function derivada_derecha(f,h,x₀)
    (f(x₀+h)-f(x₀)) / h
end

function derivada_simetrica(f,h,x₀)
    (f(x₀+h)-f(x₀-h)) / (2*h)
end

function derivada_imaginaria(f,h,x₀)
    imag((f(x₀+h*im))) / h
    
end

ies = Float64[]
for i = 1:120
    push!(ies,10.0^(-0.125 * i))
end

println()

plot(title = "Precisión de diferentes derivadas",
    yscale = :log10, 
    xscale = :log10, 
    xlims = (1e-16,1), 
    ylims = (1e-16,1),
    legend = :left, 
    xlabel= "h", 
    ylabel = "Error")

plot!(ies,abs.(derivada_simetrica.(f,ies,1).-9), label = "Simétrica",)
plot!(ies,abs.(derivada_derecha.(f,ies,1).-9), label = "Derecha")


plot!(ies,abs.(derivada_imaginaria.(f,ies,1).-9), label = "Compleja")


