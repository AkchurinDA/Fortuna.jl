@testset "Nataf Transformation: Distorted Correlation Matrix #1" begin
    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Normal", "Moments", [0, 1])
    X₂ = generaterv("Normal", "Moments", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)
    ρᶻ = NatafObject.ρᶻ

    # Test the results:
    @test ρᶻ == [1 0; 0 1]
end

@testset "Nataf Transformation: Distorted Correlation Matrix #2" begin
    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Uniform", "Parameters", [0, 1])
    X₂ = generaterv("Uniform", "Parameters", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)
    ρᶻ = NatafObject.ρᶻ

    # Test the results:
    @test ρᶻ == [1 0; 0 1]
end
