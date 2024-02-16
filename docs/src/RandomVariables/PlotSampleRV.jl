using Fortuna, Random
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

Random.seed!(1)

X = generaterv("Gamma", "M", [10, 1.5])

XSamples = samplerv(X, 5000, ITS())

begin
    F = Figure(size = 72 .* (8, 6), fonts = (; regular = texfont()), fontsize = 14)
    A = Axis(F[1, 1], 
            xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
            xminorticksvisible = true, yminorticksvisible = true,
            xgridvisible = true, ygridvisible = true,
            xminorgridvisible = true, yminorgridvisible = true,
            xlabel = L"$x$", ylabel = "PDF",
            limits = (5.0, 17.5, 0, nothing))

    hist!(XSamples,
        color = :steelblue,
        strokecolor = :black, strokewidth = 0.25,
        bins = 25, normalization = :pdf)
    
    display(F)
end

save("Sample-RVariable.svg", F)

X₁  = generaterv("Gamma", "M", [10, 1.5])
X₂  = generaterv("Gumbel", "M", [15, 2.5])
X   = [X₁, X₂]

XSamples = samplerv(X, 5000, ITS())

begin
    F = Figure(size = 72 .* (6, 6), fonts = (; regular = texfont()), fontsize = 14)
    A = Axis(F[1, 1], 
            xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
            xminorticksvisible = true, yminorticksvisible = true,
            xgridvisible = true, ygridvisible = true,
            xminorgridvisible = true, yminorgridvisible = true,
            xlabel = L"$x_{1}$", ylabel = L"$x_{2}$",
            limits = (5.0, 17.5, 7.5, 30.0),
            aspect = 1)

    scatter!(XSamples[:, 1], XSamples[:, 2],
        color = (:steelblue, 0.75),
        strokecolor = (:black, 0.75), strokewidth = 0.25,
        markersize = 6)
    
    display(F)
end

save("Sample-RVector-U.svg", F)

X₁  = generaterv("Gamma", "M", [10, 1.5])
X₂  = generaterv("Gumbel", "M", [15, 2.5])
X   = [X₁, X₂]

ρˣ = [1 0.90; 0.90 1]

TransformationObject = NatafTransformation(X, ρˣ)

XSamples, ZSamples, USamples = samplerv(TransformationObject, 5000, ITS())

begin
    F = Figure(size = 72 .* (6, 6), fonts = (; regular = texfont()), fontsize = 14)
    A = Axis(F[1, 1], 
            xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
            xminorticksvisible = true, yminorticksvisible = true,
            xgridvisible = true, ygridvisible = true,
            xminorgridvisible = true, yminorgridvisible = true,
            xlabel = L"$x_{1}$", ylabel = L"$x_{2}$",
            limits = (5.0, 17.5, 7.5, 30.0),
            aspect = 1)

    scatter!(XSamples[:, 1], XSamples[:, 2],
        color = (:steelblue, 0.75),
        strokecolor = (:black, 0.75), strokewidth = 0.25,
        markersize = 6)
    
    display(F)
end

save("Sample-RVector-C.svg", F)