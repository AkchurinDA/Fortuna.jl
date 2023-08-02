using Fortuna, Distributions
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type=:svg)

# Define the distribution types and moments of the marginal random variables:
RVTypesₓ = ["Gamma", "Gumbel"]
μₓ = [5, 7.5]
σₓ = [1, 2.5]

# Generate each marginal random variable and place it into a random vector:
X = Vector{Distribution}(undef, length(RVTypesₓ))
for i in eachindex(X)
    X[i] = generaterv(RVTypesₓ[i], "Moments", [μₓ[i], σₓ[i]])
end

# Define the correlation between the marginal random variables:
ρˣ = [1 0.9; 0.9 1]

# Perform Nataf transformation of the correlated marginal random variables:
NatafObject = NatafTransformation(X, ρˣ)

# Generate samples of the correlated marginal random variables:
NumSamples = 10^3
XSamples, ZSamples, USamples = samplerv(NatafObject, NumSamples)

# Plot:
begin
    F = Figure(resolution=72 .* (18, 6), fontsize=14, fonts=(; regular=texfont()))
    A = Axis(F[1, 1],
        title=L"Samples of correlated non-normal random variables $X$",
        xlabel=L"$x_1$", ylabel=L"$x_2$",
        aspect=1)
    scatter!(XSamples, label=L"$X$",
        color=:black,
        markersize=3)
    limits!(-10, 20, -10, 20)
    A = Axis(F[1, 2],
        title=L"Samples of correlated standard normal random variables $Z$",
        xlabel=L"$z_1$", ylabel=L"$z_2$",
        aspect=1)
    scatter!(ZSamples, label=L"$Z$",
        color=:black,
        markersize=3)
    limits!(-10, 20, -10, 20)
    A = Axis(F[1, 3],
        title=L"Samples of uncorrelated standard normal random variables $U$",
        xlabel=L"$u_1$", ylabel=L"$u_2$",
        aspect=1)
    scatter!(USamples, label=L"$U$",
        color=:black,
        markersize=3)
    limits!(-10, 20, -10, 20)
    F
end