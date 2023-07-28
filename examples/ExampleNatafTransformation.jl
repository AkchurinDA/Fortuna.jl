#=
Author: Damir Akchurin
Date:   July 14, 2023
=#

# Preamble:
using Fortuna, Distributions

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

# Output the distorted correlation matrix:
display(NatafObject.ρᶻ)

# Output the lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix:
display(NatafObject.L)

# Output the inverse of lower triangular matrix of the Cholesky decomposition of the distorted correlation matrix:
display(NatafObject.L⁻¹)

# Generate samples of the correlated marginal random variables:
NumSamples = 10^3
XSamples, ZSamples, USamples = samplerv(NatafObject, NumSamples)

# Plot generated samples and the joint PDF at these samples:
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type=:svg)

begin
    F = Figure(resolution=72 .* (9, 9), fonts=(; regular=texfont()))

    A = Axis(F[1, 1],
        xlabel=L"$x_1$, $z_1$, $u_1$", ylabel=L"$x_2$, $z_2$, $u_2$",
        xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
        limits=(-10, 30, -10, 30),
        aspect=1)

    # X:
    scatter!(XSamples, label=L"$X$",
        color=(:red, 0.25),
        markersize=5)
    # Z:
    scatter!(ZSamples, label=L"$Z$",
        color=(:blue, 0.25),
        markersize=5)
    # U:
    scatter!(USamples, label=L"$U$",
        color=(:green, 0.25),
        markersize=5)

    axislegend(position=:rt)
    F
end