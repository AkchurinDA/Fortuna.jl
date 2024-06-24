using CairoMakie, MathTeXEngine
using Makie.GeometryBasics
CairoMakie.activate!(type = :svg)

# Integrand evaluations:
begin
    F = Figure(size = 72 .* (8, 6), fonts = (; regular = texfont()), fontsize = 18)

    A = Axis(F[1, 1],
        xlabel = L"Number of dimensions, $N$",
        ylabel = "Required number of integrand evaluations",
        yscale = log10,
        limits = (1, 10, 1, nothing),
        aspect = 4 / 3)

    N = 1:10
    scatterlines!(A, N, 3 .^ N, 
        label = "1-point Gauss-Krondrod quadrature rule\napplied in each dimensions")
    scatterlines!(A, N, 2 .^ N .+ 2 .* N .^ 2 .+ 2 .* N .+ 1,
        label = "Genz-Malik cubature rule")

    axislegend(A, position = :rb)

    display(F)
end

save("temp/IntegrandEvaluations.pdf", F)

# Subdivision of the initial integration region:
begin
    F = Figure(size = 72 .* (12, 6), fonts = (; regular = texfont()), fontsize = 18)

    A = Axis(F[1, 1],
        title  = L"Integration region$$",
        limits = (-1.1, +1.1, -1.1, +1.1),
        aspect = 1)

    hidedecorations!(A)
    hidespines!(A)

    arc!(A, (0, 0), 1, 0, 2 * π,
        color = :black)

    scatterlines!(A, [0, 0.5 * cos(3 * π / 4)], [0, 0.5 * sin(3 * π / 4)],
        color = :black)

    text!(A, (0.5 * cos(3 * π / 4), 0.5 * sin(3 * π / 4)), 
        text = L"$(t, \varphi_1)$",
        align = (:right, :bottom),
        color = :black)

    A = Axis(F[1, 2],
        title  = L"Subdivided integration region$$",
        limits = (-1.1, +1.1, -1.1, +1.1),
        aspect = 1)

    hidedecorations!(A)
    hidespines!(A)

    # Define the number of subdivisions along each dimension:
    N = [2, 8]

    # Define the lower and upper bounds of each integration region:
    Δ           = [1, 2 * π] ./ N
    tLower      = collect(range(0,     1 - Δ[1], step = Δ[1]))
    φLower      = collect(range(0, 2 * π - Δ[2], step = Δ[2]))
    LowerBounds = vec(collect(Iterators.product(tLower, φLower)))
    tUpper      = collect(range(Δ[1],     1, step = Δ[1]))
    φUpper      = collect(range(Δ[2], 2 * π, step = Δ[2]))
    UpperBounds = vec(collect(Iterators.product(tUpper, φUpper)))

    for i in 1:prod(N)
        tLower, φLower = LowerBounds[i]
        tUpper, φUpper = UpperBounds[i]
        tMidpoint      = (tLower + tUpper) / 2
        φMidpoint      = (φLower + φUpper) / 2

        EDGE1 = Vector{Tuple{Float64, Float64}}(undef, 100)
        EDGE2 = Vector{Tuple{Float64, Float64}}(undef, 100)
        EDGE3 = Vector{Tuple{Float64, Float64}}(undef, 100)
        EDGE4 = Vector{Tuple{Float64, Float64}}(undef, 100)

        tRANGE1 = collect(range(tLower, tUpper, 100))
        tRANGE2 = collect(range(tUpper, tLower, 100))
        φRANGE1 = collect(range(φLower, φUpper, 100))
        φRANGE2 = collect(range(φUpper, φLower, 100))
        for j in 1:100
            EDGE1[j] = (    tLower * cos(φRANGE1[j]),     tLower * sin(φRANGE1[j]))
            EDGE2[j] = (tRANGE1[j] *     cos(φUpper), tRANGE1[j] *     sin(φUpper))
            EDGE3[j] = (    tUpper * cos(φRANGE2[j]),     tUpper * sin(φRANGE2[j]))
            EDGE4[j] = (tRANGE2[j] *     cos(φLower), tRANGE2[j] *     sin(φLower))
        end
        
        poly!(A, Polygon(Point2f[vcat(EDGE1, EDGE2, EDGE3, EDGE4)...]),
            color       = tMidpoint,
            colorrange  = (0, 1),
            colormap    = cgrad([:deepskyblue, :white]),
            strokecolor = :black,
            strokewidth = 0.5)

        text!(A, (tMidpoint * cos(φMidpoint), tMidpoint * sin(φMidpoint)),
            text     = "$(i)",
            align    = (:center, :center),
            color    = :black,
            fontsize = 12)
    end

    arc!(A, (0, 0), 1, 0, 2 * π,
        color = :black)

    display(F)
end

save("temp/IntegrationRegionSubdivision.pdf", F)