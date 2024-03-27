using Fortuna
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

# Define random vector:
X₁  = randomvariable("Normal", "M", [0, 1])
X₂  = randomvariable("Normal", "M", [0, 1])
X   = [X₁, X₂]

# Define correlation matrix:
ρˣ  = [1 0; 0 1]

# Define limit state function:
a = 5.50
b = 0.02
c = 5 / 6
d = π / 3
g(x::Vector) = a - x[2] + b * x[1] ^ 3 + c * sin(d * x[1])

# Define reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Solve reliability problem using Subset Simulation Method:
Solution = solve(Problem, SSM())

# Plot:
xRange₁ = range(-9, +9, 500)
xRange₂ = range(-9, +9, 500)
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

    contour!(xRange₁, xRange₂, gSamples,
        levels = [0],
        color = (:black, 0.25))

    contourf!(xRange₁, xRange₂, gSamples,
        levels = [0],
        extendhigh = (:green, 0.25), extendlow = (:red, 0.25))

    for i in eachindex(Solution.CSubset)
        contour!(xRange₁, xRange₂, gSamples,
            levels = [Solution.CSubset[i]],
            color = (:black, 0.25),
            linestyle = :dash)

        scatter!(Solution.XSamplesSubset[i][1, 1:100:end], Solution.XSamplesSubset[i][2, 1:100:end],
            alpha = 0.5, 
            strokecolor = (:black, 0.5), strokewidth = 0.25,
            markersize = 6)
    end

    text!(-6.0, 6.0, text = L"$g(\vec{\mathbf{x}}) \leq 0$",
        color = :black, 
        align = (:center, :bottom), fontsize = 12)

    text!(-6.0, 6.0, text = "(Failure domain)",
        color = :black, 
        align = (:center, :top), fontsize = 12)

    text!(6.0, -6.0, text = L"$g(\vec{\mathbf{x}}) > 0$",
        color = :black, 
        align = (:center, :bottom), fontsize = 12)

    text!(6.0, -6.0, text = "(Safe domain)",
        color = :black, 
        align = (:center, :top), fontsize = 12)
    
    display(F)
end