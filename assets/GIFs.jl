using Fortuna, Distributions
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type=:svg)

# # Define the distribution types and moments of the marginal random variables:
# RVTypesₓ = ["Gamma", "Gumbel"]
# μₓ = [2.5, 0]
# σₓ = [1, 2.5]

# # Generate each marginal random variable and place it into a random vector:
# X = Vector{Distribution}(undef, length(RVTypesₓ))
# for i in eachindex(X)
#     X[i] = generaterv(RVTypesₓ[i], "Moments", [μₓ[i], σₓ[i]])
# end

# # Define the correlation between the marginal random variables:
# ρˣ = [1 0.9; 0.9 1]

# # Perform Nataf transformation of the correlated marginal random variables:
# NatafObject = NatafTransformation(X, ρˣ)

# # Generate samples of the correlated marginal random variables:
# NumSamples = 1000
# XSamples, ZSamples, USamples = samplerv(NatafObject, NumSamples)

# # Plot:
# # XPoints = Observable(Point2f[(XSamples[1, 1], XSamples[1, 2])])
# # ZPoints = Observable(Point2f[(ZSamples[1, 1], ZSamples[1, 2])])
# # UPoints = Observable(Point2f[(USamples[1, 1], USamples[1, 2])])
# XPoints = Observable(Point2f[(NaN, NaN)])
# ZPoints = Observable(Point2f[(NaN, NaN)])
# UPoints = Observable(Point2f[(NaN, NaN)])

# F = Figure(resolution=72 .* (18, 6), fontsize=16, fonts=(; regular=texfont()))
# A = Axis(F[1, 1],
#     title=L"Correlated non-normal random variables $\vec{X}$",
#     xlabel=L"$x_1$", ylabel=L"$x_2$",
#     xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
#     xminorticksvisible=true, yminorticksvisible=true,
#     xminorgridvisible=true, yminorgridvisible=true,
#     aspect=1)
# scatter!(XPoints, label=L"$X$",
#     color=:darkgreen,
#     markersize=3)
# limits!(-10, 10, -10, 10)
# A = Axis(F[1, 2],
#     title=L"Correlated standard normal random variables $\vec{Z}$",
#     xlabel=L"$z_1$", ylabel=L"$z_2$",
#     xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
#     xminorticksvisible=true, yminorticksvisible=true,
#     xminorgridvisible=true, yminorgridvisible=true,
#     aspect=1)
# scatter!(ZPoints, label=L"$Z$",
#     color=:darkblue,
#     markersize=3)
# limits!(-10, 10, -10, 10)
# A = Axis(F[1, 3],
#     title=L"Uncorrelated standard normal random variables $\vec{U}$",
#     xlabel=L"$u_1$", ylabel=L"$u_2$",
#     xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
#     xminorticksvisible=true, yminorticksvisible=true,
#     xminorgridvisible=true, yminorgridvisible=true,
#     aspect=1)
# scatter!(UPoints, label=L"$U$",
#     color=:darkred,
#     markersize=3)
# limits!(-10, 10, -10, 10)

# record(F, "assets/NatafTransformation.gif", 1:NumSamples, framerate=200) do i
#     XPoints[] = push!(XPoints[], Point2f(XSamples[i, :]))
#     ZPoints[] = push!(ZPoints[], Point2f(ZSamples[i, :]))
#     UPoints[] = push!(UPoints[], Point2f(USamples[i, :]))
# end

# Define a random vector of correlated marginal distributions:
X₁ = generaterv("Normal", "Moments", [10, 2])
X₂ = generaterv("Normal", "Moments", [20, 5])
X = [X₁, X₂]
ρˣ = [1 0.5; 0.5 1]

# Perform Nataf transformation of the correlated marginal random variables:
NatafObject = NatafTransformation(X, ρˣ)

# Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
G(x::Vector) = x[1]^2 - 2 * x[2]

# Define reliability problems:
Problem = ReliabilityProblem(X, ρˣ, G)

# Perform the reliability analysis using MCFOSM:
Solution = analyze(Problem, FORM(iHLRF(ϵ₁=10e-3, ϵ₂=10e-3)))

# Point on a limit state function:
x = range(0, 20, 100)
GPointsX = hcat(x, (1 / 2) .* x .^ 2)
GPointsU = transformsamples(NatafObject, GPointsX, "X2U")

# Plot:
XPoints = Observable(Point2f[(NaN, NaN)])
UPoints = Observable(Point2f[(NaN, NaN)])

F = Figure(resolution=72 .* (12, 6), fontsize=16, fonts=(; regular=texfont()))
A = Axis(F[1, 1],
    title=L"Space of non-normal random variables $\vec{X}$",
    xlabel=L"$x_1$", ylabel=L"$x_2$",
    xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
    xminorticksvisible=true, yminorticksvisible=true,
    xminorgridvisible=true, yminorgridvisible=true,
    aspect=1)
band!(GPointsX[:, 1], 15 * ones(100,), GPointsX[:, 2], label=L"Safe domain$$",
    color=(:darkgreen, 0.25))
band!(GPointsX[:, 1], GPointsX[:, 2], 25 * ones(100,), label=L"Failure domain$$",
    color=(:darkred, 0.25))
lines!(GPointsX, label=L"Failure boundary$$",
    color=:black)
scatterlines!(XPoints, label=L"Design point $\vec{x}^*$",
    color=:black,
    markersize=6)
limits!(5, 15, 15, 25)
axislegend(A, position=:rt)
A = Axis(F[1, 2],
    title=L"Space of standard normal random variables $\vec{U}$",
    xlabel=L"$u_1$", ylabel=L"$u_2$",
    xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
    xminorticksvisible=true, yminorticksvisible=true,
    xminorgridvisible=true, yminorgridvisible=true,
    aspect=1)
band!(vcat(GPointsU[1:55, 1], collect(range(GPointsU[55, 1], 5, 100))), -5 * ones(155,), vcat(GPointsU[1:55, 2], 5 * ones(100,)), label=L"Safe domain$$",
    color=(:darkgreen, 0.25))
band!(GPointsU[1:55, 1], GPointsU[1:55, 2], 5 * ones(55,), label=L"Failure domain$$",
    color=(:darkred, 0.25))
lines!(GPointsU, label=L"Failure boundary$$",
    color=:black)
scatterlines!(UPoints, label=L"Design point $\vec{u}^*$",
    color=:black,
    markersize=6)
limits!(-5, 5, -5, 5)
axislegend(A, position=:rt)

record(F, "assets/FORM.gif", 1:size(Solution.x, 2), framerate=1) do i
    XPoints[] = push!(XPoints[], Point2f(Solution.x[:, i]))
    UPoints[] = push!(UPoints[], Point2f(Solution.u[:, i]))
end