using Fortuna, Random
using CairoMakie, MathTeXEngine
CairoMakie.activate!(type = :svg)

Random.seed!(1)

X₁  = generaterv("Normal", "M", [0, 1])
X₂  = generaterv("Normal", "M", [0, 1])
X   = [X₁, X₂]
ρˣ  = [1 0; 0 1]
NatafObject = NatafTransformation(X, ρˣ)

β               = 3
g(x::Vector)    = β * sqrt(2) - x[1] - x[2]

Problem = ReliabilityProblem(X, ρˣ, g)

q           = MvNormal([β / sqrt(2), β / sqrt(2)], [1 -0.5; -0.5 1])
NumSamples  = vcat(1000:500:9000, 10000:5000:90000, 100000:50000:900000, 1000000)
PoFValues   = Matrix{Float64}(undef, length(NumSamples), 2)
for i in eachindex(NumSamples)
    display(i)
    PoFValues[i, 1] = analyze(Problem, MCS(NumSamples[i], ITS())).PoF
    PoFValues[i, 2] = analyze(Problem, IS(q, NumSamples[i])).PoF
end

PoF = cdf(Normal(), -β)
begin
    F = Figure(size = 72 .* (8, 6), fonts = (; regular = texfont()), fontsize = 14)

    A = Axis(F[1, 1],
        xlabel = L"$N$", ylabel = L"$P_{f}$",
        xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        xscale = log10,
        limits = (NumSamples[1], NumSamples[end], 0, 0.003),
        aspect = 4 / 3)

    hlines!(cdf(Normal(), -β), label = L"Exact solution$$",
        color = :grey, 
        linestyle = :dash, linewidth = 1)

    scatterlines!(NumSamples, PoFValues[:, 1], label = L"MCS()$$",
        color = :steelblue,
        linestyle = :solid, linewidth = 1,
        markersize = 6)
    
    scatterlines!(NumSamples, PoFValues[:, 2], label = L"IS()$$",
        color = :crimson,
        linestyle = :solid, linewidth = 1,
        markersize = 6)

    axislegend(position = :rt, fontsize = 12)

    display(F)
end

save("ImportanceSampling-1.svg", F)

xRange₁     = range(-3, +6, 500)
xRange₂     = range(-3, +6, 500)
gSamples    = [g([x₁, x₂]) for x₁ in xRange₁, x₂ in xRange₂]
fSamples    = [jointpdf(NatafObject, [x₁, x₂]) for x₁ in xRange₁, x₂ in xRange₂]
qSamples    = [pdf(q, [x₁, x₂]) for x₁ in xRange₁, x₂ in xRange₂]

begin
    F = Figure(size = 72 .* (6, 6), fonts = (; regular = texfont()), fontsize = 14)

    A = Axis(F[1, 1],
        xlabel = L"$x_{1}$", ylabel = L"$x_{2}$",
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        limits = (minimum(xRange₁), maximum(xRange₁), minimum(xRange₂), maximum(xRange₂)),
        aspect = 1)

    contour!(xRange₁, xRange₂, gSamples, label = L"$g(\mathbf{\vec{x}}) = 0$",
        levels = [0],
        color = :black)
    
    Samples = analyze(Problem, MCS(500, ITS())).Samples
    contour!(xRange₁, xRange₂, fSamples, label = L"$f_{\mathbf{\vec{X}}}(\mathbf{\vec{x}})$",
        levels = 10,
        color = (:grey, 0.5))
    scatter!(Samples[:, 1], Samples[:, 2], label = L"$f_{\mathbf{\vec{X}}}(\mathbf{\vec{x}})$",
        color = (:steelblue, 0.5),
        markersize = 6)
    
    Samples = analyze(Problem, IS(q, 500)).Samples
    contour!(xRange₁, xRange₂, qSamples, label = L"$q_{\mathbf{\vec{X}}}(\mathbf{\vec{x}})$",
        levels = 10,
        color = (:grey, 0.5))
    scatter!(Samples[:, 1], Samples[:, 2], label = L"$q_{\mathbf{\vec{X}}}(\mathbf{\vec{x}})$",
        color = (:crimson, 0.5),
        markersize = 6)

    axislegend(position = :rt, fontsize = 12, merge = true)

    display(F)
end

save("ImportanceSampling-2.svg", F)