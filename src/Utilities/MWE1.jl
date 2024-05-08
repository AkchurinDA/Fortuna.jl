using Makie, CairoMakie

@recipe(ReliabilityPlot2) do scene
    Attributes(
        meshes = 1000,
    )
end

function Makie.plot!(p::ReliabilityPlot2)
    contour!(p, rand(10, 10))
end

reliabilityplot2(rand(10))