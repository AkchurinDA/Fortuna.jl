@testset "Sampling: Uncorrelated random variables" begin
    # Generate a random vector X with uncorrelated marginal random variables X₁ and X₂:
    X₁ = generaterv("Gamma", "M", [10, 1.5])
    X₂ = generaterv("Gumbel", "M", [15, 2.5])
    X = [X₁, X₂]

    # Generate samples:
    NumSamples = 5 * 10^6
    XSamplesITS = samplerv(X, NumSamples, ITS())
    XSamplesLHS = samplerv(X, NumSamples, LHS())

    # Test the results:
    @test isapprox(mean(XSamplesITS[:, 1]), 10, rtol=0.01)
    @test isapprox(mean(XSamplesITS[:, 2]), 15, rtol=0.01)
    @test isapprox(std(XSamplesITS[:, 1]), 1.5, rtol=0.01)
    @test isapprox(std(XSamplesITS[:, 2]), 2.5, rtol=0.01)
    @test isapprox(mean(XSamplesLHS[:, 1]), 10, rtol=0.01)
    @test isapprox(mean(XSamplesLHS[:, 2]), 15, rtol=0.01)
    @test isapprox(std(XSamplesLHS[:, 1]), 1.5, rtol=0.01)
    @test isapprox(std(XSamplesLHS[:, 2]), 2.5, rtol=0.01)
end

@testset "Sampling: Correlated random variables" begin
    # Generate a random vector X of correlated marginal distributions:
    X₁ = generaterv("Gamma", "M", [10, 1.5])
    X₂ = generaterv("Gumbel", "M", [15, 2.5])
    X = [X₁, X₂]
    ρˣ = [1 0.75; 0.75 1]

    # Perform Nataf transformation of the correlated marginal random variables:
    NatafObject = NatafTransformation(X, ρˣ)

    # Generate samples:
    NumSamples = 5 * 10^6
    XSamples, ZSamples, USamples = samplerv(NatafObject, NumSamples)

    # Test the results:
    @test isapprox(mean(XSamples[:, 1]), 10, rtol=0.01)
    @test isapprox(mean(XSamples[:, 2]), 15, rtol=0.01)
    @test isapprox(std(XSamples[:, 1]), 1.5, rtol=0.01)
    @test isapprox(std(XSamples[:, 2]), 2.5, rtol=0.01)
    @test isapprox(cor(XSamples, dims=1), [1 0.75; 0.75 1], rtol=0.01)
end