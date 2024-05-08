import Makie

export reliabilityplot
export reliabilityplot!

Makie.@recipe(ReliabilityPlot, Problem, Solution) do scene
    Makie.Attributes(
        scale     = 1.5,
        meshsize  = 1000,
        safecolor = :green,
        failcolor = :red)
end

function Makie.plot!(P::ReliabilityPlot)
    # Extract the problem parameters:
    Problem = P[:Problem]
    X       = Problem[].X
    g       = Problem[].g

    # Extract the solution parameters:
    Solution = P[:Solution]
    x        = Solution[].x
    β        = Solution[].β

    # Extract the plotting attributes:
    scale     = P.attributes.scale[]
    meshsize  = P.attributes.meshsize[]
    safecolor = P.attributes.safecolor[]
    failcolor = P.attributes.failcolor[]

    # Define the range of interest:
    xS₁, xS₂ = x[1, 1  ], x[2, 1  ]
    xF₁, xF₂ = x[1, end], x[2, end]
    Δx₁, Δx₂ = xS₂ - xS₁, xF₂ - xF₁
    xRange₁  = collect(range(xS₁ - scale * Δx₁, xF₁ + scale * Δx₁, length = meshsize))
    xRange₂  = collect(range(xS₂ - scale * Δx₂, xF₂ + scale * Δx₂, length = meshsize))

    # Compute the limit state function:
    gValue = [g([x₁, x₂]) for x₁ in xRange₁, x₂ in xRange₂]

    # Plot the limit state function:
    Makie.contour!(xRange₁, xRange₂, gValue;
        levels = [0],
        color = :black)

    # Return the plot:
    return P
end