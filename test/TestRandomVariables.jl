@testset "Random Variables" begin
    @testset "Exponential" begin
        Moments = zeros(10, 2)
        for i in 1:10
            Moments[i, 1] = i
            Moments[i, 2] = i
        end

        RandomVariables = [randomvariable("Exponential", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end

    @testset "Frechet" begin
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

    @testset "Gamma" begin
        μ = 1:10
        σ = 0.1:0.1:1
        Moments = zeros(length(μ) * length(σ), 2)
        for i in eachindex(μ)
            for j in eachindex(σ)
                Moments[(i - 1) * length(σ) + j, 1] = μ[i]
                Moments[(i - 1) * length(σ) + j, 2] = σ[j]
            end
        end

        RandomVariables = [randomvariable("Gamma", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end

    @testset "Gumbel" begin
        μ = 1:10
        σ = 0.1:0.1:1
        Moments = zeros(length(μ) * length(σ), 2)
        for i in eachindex(μ)
            for j in eachindex(σ)
                Moments[(i - 1) * length(σ) + j, 1] = μ[i]
                Moments[(i - 1) * length(σ) + j, 2] = σ[j]
            end
        end

        RandomVariables = [randomvariable("Gumbel", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end

    @testset "LogNormal" begin
        μ = 1:10    
        σ = 0.1:0.1:1
        Moments = zeros(length(μ) * length(σ), 2)
        for i in eachindex(μ)
            for j in eachindex(σ)
                Moments[(i - 1) * length(σ) + j, 1] = μ[i]
                Moments[(i - 1) * length(σ) + j, 2] = σ[j]
            end
        end

        RandomVariables = [randomvariable("LogNormal", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end
    
    @testset "Normal" begin
        μ = 1:10
        σ = 0.1:0.1:1
        Moments = zeros(length(μ) * length(σ), 2)
        for i in eachindex(μ)
            for j in eachindex(σ)
                Moments[(i - 1) * length(σ) + j, 1] = μ[i]
                Moments[(i - 1) * length(σ) + j, 2] = σ[j]
            end
        end

        RandomVariables = [randomvariable("Normal", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end

    @testset "Poisson" begin
        Moments = zeros(10, 2)
        for i in 1:10
            Moments[i, 1] = i
            Moments[i, 2] = sqrt(i)
        end

        RandomVariables = [randomvariable("Poisson", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end

    @testset "Uniform" begin
        μ = 1:10
        σ = 0.1:0.1:1
        Moments = zeros(length(μ) * length(σ), 2)
        for i in eachindex(μ)
            for j in eachindex(σ)
                Moments[(i - 1) * length(σ) + j, 1] = μ[i]
                Moments[(i - 1) * length(σ) + j, 2] = σ[j]
            end
        end

        RandomVariables = [randomvariable("Uniform", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end

    @testset "Weibull" begin
        μ = 1:10
        σ = 0.1:0.1:1
        Moments = zeros(length(μ) * length(σ), 2)
        for i in eachindex(μ)
            for j in eachindex(σ)
                Moments[(i - 1) * length(σ) + j, 1] = μ[i]
                Moments[(i - 1) * length(σ) + j, 2] = σ[j]
            end
        end

        RandomVariables = [randomvariable("Weibull", "M", Moments[i, :]) for i in axes(Moments, 1)]

        @test isapprox(hcat(mean.(RandomVariables), std.(RandomVariables)), Moments, rtol = 10 ^ (-9))
    end
end