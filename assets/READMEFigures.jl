using Fortuna, Distributions
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type=:svg)

# Define the distribution types and moments of the marginal random variables:
RVTypesₓ = ["Gamma", "Gumbel"]
μₓ = [10, 15]
σₓ = [1.5, 2.5]

# Generate each marginal random variable and place it into a random vector:
X = Vector{Distribution}(undef, length(RVTypesₓ))
for i in eachindex(X)
    X[i] = generaterv(RVTypesₓ[i], "Moments", [μₓ[i], σₓ[i]])
end

# Define the correlation between the marginal random variables:
ρˣ = [1 0.75; 0.75 1]

# Perform Nataf transformation of the correlated marginal random variables:
NatafObject = NatafTransformation(X, ρˣ)

# Generate samples of the correlated marginal random variables:
NumSamples = 2 * 10^3
XSamples, ZSamples, USamples = samplerv(NatafObject, NumSamples)

# Plot:
begin
    F = Figure(resolution=72 .* (6, 6), fontsize=14, fonts=(; regular=texfont()))
    A = Axis(F[1, 1],
        title=L"Nataf Transformation$$",
        xlabel=L"$x_1$, $z_1$, $u_1$", ylabel=L"$x_2$, $z_2$, $u_2$")
    scatter!(XSamples, label=L"$X$",
        markersize=3)
    scatter!(ZSamples, label=L"$Z$",
        markersize=3)
    scatter!(USamples, label=L"$U$",
        markersize=3)
    limits!(-10, 30, -10, 30)
    axislegend(A, position=:rt)
    F
end