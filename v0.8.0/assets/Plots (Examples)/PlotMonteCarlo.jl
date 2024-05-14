using Fortuna
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

X₁  = randomvariable("Normal", "M", [0, 1])
X₂  = randomvariable("Normal", "M", [0, 1])
X   = [X₁, X₂]

ρˣ = [1 0; 0 1]

NatafObject = NatafTransformation(X, ρˣ)

β               = 3
g(x::Vector)    = β * sqrt(2) - x[1] - x[2]

Problem = ReliabilityProblem(X, ρˣ, g)
Solution = solve(Problem, MC())

xRange₁ = range(-3, +6, 500)
xRange₂ = range(-3, +6, 500)
gSamples = [g([x₁, x₂]) for x₁ in xRange₁, x₂ in xRange₂]

begin
    F = Figure(size = 72 .* (6, 6), fonts = (; regular = texfont()), fontsize = 14)

    A = Axis(F[1, 1],
        xlabel = L"$x_{1}$", ylabel = L"$x_{2}$",
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        limits = (minimum(xRange₁), maximum(xRange₁), minimum(xRange₂), maximum(xRange₂)),
        aspect = 1)

    contour!(xRange₁, xRange₂, gSamples, label = L"$g(\vec{\mathbf{x}}) = 0$",
        levels = [0],
        color = (:black, 0.25))

    contourf!(xRange₁, xRange₂, gSamples,
        levels = [0],
        extendhigh = (:green, 0.25), extendlow = (:red, 0.25))

    scatter!(Solution.Samples[1, 1:100:end], Solution.Samples[2, 1:100:end],
        color = (:steelblue, 0.5), 
        strokecolor = (:black, 0.5), strokewidth = 0.25,
        markersize = 6)
    
    text!(4.5, 5.5, text = L"$g(\vec{\mathbf{x}}) \leq 0$",
        color = :black, 
        align = (:center, :center), fontsize = 12)

    text!(4.5, 5.0, text = "(Failure domain)",
        color = :black, 
        align = (:center, :center), fontsize = 12)

    text!(4.5, -2.0, text = L"$g(\vec{\mathbf{x}}) > 0$",
        color = :black, 
        align = (:center, :center), fontsize = 12)

    text!(4.5, -2.5, text = "(Safe domain)",
        color = :black, 
        align = (:center, :center), fontsize = 12)
    
    display(F)
end

save("docs/src/assets/Plots (Examples)/MonteCarlo-1.svg", F)