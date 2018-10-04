
using SpecialFunctions
using Plots

B(a, b) = gamma(a) * gamma(b) / gamma(a + b)
pdf_ψ(x) = 1 / B(1, 1) * exp(x) / (1 + exp(x))^2

plot(pdf_ψ, legend = false)
savefig("julia-scripts/week-03_plot1.png")

pdf_ψ(x) = 1 / gamma(1) * exp(x - exp(x))
