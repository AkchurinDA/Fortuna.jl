using Fortuna, Random
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

N(PoF, V) = (1 - PoF) / (PoF * V ^ 2)

VRange          = collect(range(0.01, 1, 1000))
PoFRange        = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]
PoFRangeText    = ["0.1", "0.01", "0.001", "0.0001", "0.00001", "0.000001"]

begin
    F = Figure(size = 72 .* (6, 6), fonts = (; regular = texfont()), fontsize = 14)

    A = Axis(F[1, 1], 
            xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(10),
            xminorticksvisible = true, yminorticksvisible = true,
            xminorgridvisible = true, yminorgridvisible = true,
            xlabel = L"$V_{P_{f}}$", ylabel = L"$N$",
            yscale = log10,
            limits = (0, 1, nothing, nothing),
            aspect = 1)

    for i in eachindex(PoFRange)
        lines!(VRange, N.(PoFRange[i], VRange),
            color = :black)

        text!(VRange[end], N(PoFRange[i], VRange[end]),
            text = L"P_{f} = %$(PoFRangeText[i])",
            align = (:right, :top),
            color = :black, fontsize = 12)
    end
    
    display(F)
end

save("RecommendedN.svg", F)