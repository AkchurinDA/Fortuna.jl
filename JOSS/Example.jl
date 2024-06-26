# Preamble:
using Fortuna

# Define the random vector:
U₁ = randomvariable("Normal", "M", [0, 1])
U₂ = randomvariable("Normal", "M", [0, 1])
U = [U₁, U₂]

# Define the correlation matrix:
ρ = [1 0; 0 1]

# Define the limit state function:
g(u::Vector) = 0.5 * (u[1] - 2) ^ 2 - 1.5 * (u[2] - 5) ^ 3 - 3

# Define the reliability problem:
Problem = ReliabilityProblem(U, ρ, g)

# Solve the reliability problem using the FORM:
Solution = solve(Problem, FORM())
println("Geometric reliability index: ", Solution.β)
println("Failure probability: ", Solution.PoF)
# Geometric reliability index: 3.932419
# Failure probability: 4.204761E-5

# Plot the failure domain:
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

uRange₁ = -6:0.01:6
uRange₂ = -6:0.01:6

NatafObject = NatafTransformation(U, ρ)
fValues = [pdf(NatafObject, [u₁, u₂]) for u₁ in uRange₁, u₂ in uRange₂]
fValues = fValues ./ maximum(fValues)
gValues = [g([u₁, u₂]) for u₁ in uRange₁, u₂ in uRange₂]

begin
    F = Figure(size = 72 .* (12, 6), fonts = (; regular = texfont()))

    A = Axis(F[1, 1],
        xlabel = L"u_1",
        ylabel = L"u_2",
        xticks = -6:3:6,
        yticks = -6:3:6,
        xminorticks = IntervalsBetween(5),
        yminorticks = IntervalsBetween(5),
        xminorgridvisible = true,
        yminorgridvisible = true,
        limits = (-6, 6, -6, 6),
        aspect = 1)

    contourf!(A, uRange₁, uRange₂, fValues, 
        levels   = 0:0.1:1, 
        colormap = cgrad([:transparent, :teal]))

    contour!(A, uRange₁, uRange₂, fValues, 
        levels    = 0:0.1:1, 
        colormap  = cgrad([:transparent, :black]),
        linewidth = 0.5)

    contourf!(A, uRange₁, uRange₂, gValues,
        levels     = [0],
        extendlow  = (  :red, 0.1),
        extendhigh = (:green, 0.1))

    contour!(A, uRange₁, uRange₂, gValues,
        levels    = [0],
        color     = :black,
        linewidth = 0.5)

    text!(A, (+2, +5), 
        text = L"Failure domain, \\ $g(u_1, u_2) \leq 0$",
        align = (:center, :center))

    text!(A, (-2, -5), 
        text = L"Safe domain, \\ $g(u_1, u_2) > 0$",
        align = (:center, :center))

    text!(A, (+3, -3), 
        text = L"Joint PDF, \\ $f_{\vec{U}}(u_1, u_2)$",
        align = (:center, :center))

    arc!(A, (0, -3), 3, π / 12, π / 4,
        color     = :black,
        linewidth = 0.5)

    display(F)
end

save("JOSS/Example.pdf", F)