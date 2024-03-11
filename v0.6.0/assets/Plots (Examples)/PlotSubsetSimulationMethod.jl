using Fortuna
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

X₁  = randomvariable("Normal", "M", [0, 1])
X₂  = randomvariable("Normal", "M", [0, 1])
X   = [X₁, X₂]

ρˣ  = [1 0; 0 1]

NatafObject = NatafTransformation(X, ρˣ)

β               = 3
g(x::Vector)    = β * sqrt(2) - x[1] - x[2]

Problem = ReliabilityProblem(X, ρˣ, g)
Solution = solve(Problem, SSM())

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
        limits = (-3, +6, -3, +6),
        aspect = 1)

    contour!(xRange₁, xRange₂, gSamples,
        levels = [0],
        color = (:black, 0.25))

    contourf!(xRange₁, xRange₂, gSamples,
        levels = [0],
        extendhigh = (:green, 0.25), extendlow = (:red, 0.25))

    for i in eachindex(Solution.CSubset)
        contour!(xRange₁, xRange₂, gSamples,
            levels = [Solution.CSubset[i]],
            color = (:black, 0.25))

        scatter!(Solution.XSamplesSubset[i][1, 1:100:end], Solution.XSamplesSubset[i][2, 1:100:end],
            alpha = 0.5, 
            strokecolor = (:black, 0.5), strokewidth = 0.25,
            markersize = 6)
    end

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

    tooltip!(1.0, 1.0, L"\Omega_{f, 1}",
        color = :black, 
        offset = 0, outline_linewidth = 1,
        fontsize = 12)

    tooltip!(1.7, 1.7, L"\Omega_{f, 2}",
        color = :black, 
        offset = 0, outline_linewidth = 1,
        fontsize = 12)

    tooltip!(2.4, 2.4, L"\Omega_{f, 3}",
        color = :black, 
        offset = 0, outline_linewidth = 1,
        fontsize = 12)


    display(F)
end

save("docs/src/assets/Plots (Examples)/SubsetSimulationMethod-1.svg", F)