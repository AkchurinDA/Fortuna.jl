import Makie

export reliabilityplot
export reliabilityplot!

Makie.@recipe(ReliabilityPlot, Problem, Solution) do scene
    Makie.Attributes()
end

function Makie.plot!(P::ReliabilityPlot)
    # Return the plot:
    return P
end