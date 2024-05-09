using Fortuna
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

X₁  = randomvariable("Gamma", "M", [10, 1.5])
X₂  = randomvariable("Gamma", "M", [15, 2.5])
X   = [X₁, X₂]

ρˣ = [1 -0.75; -0.75 1]

NatafObject = NatafTransformation(X, ρˣ)

XSamples, ZSamples, USamples = rand(NatafObject, 10000, LHS())

begin
    F = Figure(size = 72 .* (18, 6), fonts = (; regular = texfont()), fontsize = 14)

    A = Axis(F[1, 1], 
        title = L"$X$-space",
        xlabel = L"$x_{1}$", ylabel = L"$x_{2}$",
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        limits = (0, 20, 5, 25),
        aspect = 1)
    
    scatter!(XSamples[1, :], XSamples[2, :],
        color = (:crimson, 0.5), 
        strokecolor = (:black, 0.5), strokewidth = 0.25,
        markersize = 6)
    
    A = Axis(F[1, 2], 
        title = L"$Z$-space",
        xlabel = L"$z_{1}$", ylabel = L"$z_{2}$",
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        limits = (-6, +6, -6, +6),
        aspect = 1)
    
    scatter!(ZSamples[1, :], ZSamples[2, :],
        color = (:steelblue, 0.5),
        strokecolor = (:black, 0.5), strokewidth = 0.25,
        markersize = 6)

    A = Axis(F[1, 3], 
        title = L"$U$-space",
        xlabel = L"$u_{1}$", ylabel = L"$u_{2}$",
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        limits = (-6, +6, -6, +6),
        aspect = 1)
    
    scatter!(USamples[1, :], USamples[2, :],
        color = (:forestgreen, 0.5),
        strokecolor = (:black, 0.5), strokewidth = 0.25,
        markersize = 6)

    display(F)
end

save("docs/src/assets/Plots (Examples)/NatafTransformation-1.svg", F)

xRange₁ = range(0, 20, 500)
xRange₂ = range(5, 25, 500)
fSamples = [pdf(NatafObject, [x₁, x₂]) for x₁ in xRange₁, x₂ in xRange₂]

begin
    F = Figure(size = 72 .* (6, 6), fonts = (; regular = texfont()), fontsize = 14)

    A = Axis(F[1, 1],
        xlabel = L"$x_{1}$", ylabel = L"$x_{2}$",
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        limits = (0, 20, 5, 25),
        aspect = 1)
    
    contour!(xRange₁, xRange₂, fSamples,
        levels = 25,
        colormap = (:turbo, 0.75))

    display(F)
end

save("docs/src/assets/Plots (Examples)/NatafTransformation-2.svg", F)
