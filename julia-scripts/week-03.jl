
using SpecialFunctions
using Plots

B(a, b) = gamma(a) * gamma(b) / gamma(a + b)
pdf_ψ(x) = 1 / B(1, 1) * exp(x) / (1 + exp(x))^2

plot(pdf_ψ, legend = false, xlab = "\\phi", ylab = "pdf (\\phi)")
savefig("julia-scripts/week-03_plot1.png")

pdf_ψ(x) = 1 / gamma(1) * exp(x - exp(x))

plot(pdf_ψ, legend = false, xlab = "\\psi", ylab = "pdf (\\psi)")
savefig("julia-scripts/week-03_plot2.png")
