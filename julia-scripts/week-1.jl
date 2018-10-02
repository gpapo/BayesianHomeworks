
using LinearAlgebra

binomial(5, 1)
binomial(5, 3)

b = Rational[1//2 - 8//30, 1//4 - 2//15]
A = Rational[2//5 4//5; 1//10 3//5]
ω₂⁵, ω₄⁵ = A \ b

b = Rational[1//30 - 1//30 + 1 * 1//4 - 1//2, 1//30 - 1//30]
A = Rational[6 -4; 1 -2]
ω₃, ω₄ = A \ b
