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