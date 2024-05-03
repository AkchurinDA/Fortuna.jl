@testset "Random Variables" begin
    # Exponential:
    Moments = zeros(10, 2)
    for i in 1:10
        Moments[i, 1] = i
        Moments[i, 2] = i
    end

    RandomVariables = [randomvariable("Exponential", "M", Moments[i, :]) for i in axes(Moments, 1)]

    @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))

    # Frechet:
    μ = 1:10
    σ = 0.1:0.1:1
    Moments = zeros(length(μ) * length(σ), 2)
    for i in eachindex(μ)
        for j in eachindex(σ)
            Moments[(i - 1) * length(σ) + j, 1] = μ[i]
            Moments[(i - 1) * length(σ) + j, 2] = σ[j]
        end
    end

    RandomVariables = [randomvariable("Frechet", "M", Moments[i, :]) for i in axes(Moments, 1)]

    @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
end