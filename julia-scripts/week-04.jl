# Exercise 2 -----

using Distributions

function generatevalues(nval, F, G)
    result = Vector{Float64}(undef, nval)
    ind = 1
    while ind ≤ nval
        x = rand(G)
        u = rand(Uniform(0, M * pdf(G, x)))
        if u < pdf(F, x)
            result[ind] = x
            ind += 1
        end
    end
    return result
end

y = generatevalues(1000, Normal(), Uniform())

using KernelDensity, Plots

k = kde(y)
plot(k.x, k.density)
