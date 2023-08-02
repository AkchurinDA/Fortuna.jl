@testset "Nataf Transformation: Distorted Correlation Matrix #1 (Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Normal", "Moments", [0, 1])
    X₂ = generaterv("Normal", "Moments", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)
    ρᶻ = NatafObject.ρᶻ

    # Test the results:
    @test ρᶻ == ρˣ
end

@testset "Nataf Transformation: Distorted Correlation Matrix #2 (Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Uniform", "Parameters", [0, 1])
    X₂ = generaterv("Uniform", "Parameters", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)
    ρᶻ = NatafObject.ρᶻ

    # Test the results:
    @test ρᶻ == ρˣ
end

@testset "Nataf Transformation: Distorted Correlation Matrix #3 (Non-Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Normal", "Moments", [0, 1])
    X₂ = generaterv("Normal", "Moments", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0.8; 0.8 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)
    ρᶻ = NatafObject.ρᶻ

    # Test the results:
    @test isapprox(ρᶻ, [1 0.8; 0.8 1], rtol=10^(-6))
end

@testset "Nataf Transformation: Distorted Correlation Matrix #4 (Non-Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Uniform", "Parameters", [0, 1])
    X₂ = generaterv("Uniform", "Parameters", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0.8; 0.8 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)
    ρᶻ = NatafObject.ρᶻ

    # Test the results:
    @test isapprox(ρᶻ, [1 0.8134732861515996; 0.8134732861515996 1], rtol=10^(-6))
end

@testset "Nataf Transformation: Sampling" begin
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