using Fortuna, Distributions
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type=:svg)

# Define the distribution types and moments of the marginal random variables:
RVTypesₓ = ["Gamma", "Gumbel"]
μₓ = [2.5, 0]
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
NumSamples = 1000
XSamples, ZSamples, USamples = samplerv(NatafObject, NumSamples)

# Plot:
# XPoints = Observable(Point2f[(XSamples[1, 1], XSamples[1, 2])])
# ZPoints = Observable(Point2f[(ZSamples[1, 1], ZSamples[1, 2])])
# UPoints = Observable(Point2f[(USamples[1, 1], USamples[1, 2])])
XPoints = Observable(Point2f[(NaN, NaN)])
ZPoints = Observable(Point2f[(NaN, NaN)])
UPoints = Observable(Point2f[(NaN, NaN)])

F = Figure(resolution=72 .* (18, 6), fontsize=16, fonts=(; regular=texfont()))
A = Axis(F[1, 1],
    title=L"Correlated non-normal random variables $\vec{X}$",
    xlabel=L"$x_1$", ylabel=L"$x_2$",
    xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
    xminorticksvisible=true, yminorticksvisible=true,
    xminorgridvisible=true, yminorgridvisible=true,
    aspect=1)
scatter!(XPoints, label=L"$X$",
    color=:darkgreen,
    markersize=3)
limits!(-10, 10, -10, 10)
A = Axis(F[1, 2],
    title=L"Correlated standard normal random variables $\vec{Z}$",
    xlabel=L"$z_1$", ylabel=L"$z_2$",
    xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
    xminorticksvisible=true, yminorticksvisible=true,
    xminorgridvisible=true, yminorgridvisible=true,
    aspect=1)
scatter!(ZPoints, label=L"$Z$",
    color=:darkblue,
    markersize=3)
limits!(-10, 10, -10, 10)
A = Axis(F[1, 3],
    title=L"Uncorrelated standard normal random variables $\vec{U}$",
    xlabel=L"$u_1$", ylabel=L"$u_2$",
    xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
    xminorticksvisible=true, yminorticksvisible=true,
    xminorgridvisible=true, yminorgridvisible=true,
    aspect=1)
scatter!(UPoints, label=L"$U$",
    color=:darkred,
    markersize=3)
limits!(-10, 10, -10, 10)

record(F, "assets/NatafTransformation.mp4", 1:NumSamples, framerate=200) do i
    XPoints[] = push!(XPoints[], Point2f(XSamples[i, :]))
    ZPoints[] = push!(ZPoints[], Point2f(ZSamples[i, :]))
    UPoints[] = push!(UPoints[], Point2f(USamples[i, :]))
end