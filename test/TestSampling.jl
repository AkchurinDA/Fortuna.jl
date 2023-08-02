@testset "Sampling: Correlated random variables" begin
    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Gamma", "Moments", [10, 1.5])
    X₂ = generaterv("Gumbel", "Moments", [15, 2.5])
    X = [X₁, X₂]
    ρˣ = [1 0.75; 0.75 1]

    # Perform Nataf transformation of the correlated marginal random variables:
    NatafObject = NatafTransformation(X, ρˣ)

    # Generate samples:
    NumSamples = 10^6
    XSamples, ZSamples, USamples = samplerv(NatafObject, NumSamples)

    # Test the results:
    @test isapprox(mean(XSamples[:, 1]), 10, rtol=0.01)
    @test isapprox(mean(XSamples[:, 2]), 15, rtol=0.01)
    @test isapprox(std(XSamples[:, 1]), 1.5, rtol=0.01)
    @test isapprox(std(XSamples[:, 2]), 2.5, rtol=0.01)
    @test isapprox(cor(XSamples, dims=1), [1 0.75; 0.75 1], rtol=0.01)
end