@testset "Reliability Analysis: MCFOSM" begin
    # Example from "Structural and System Reliability" by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Normal", "Moments", [10, 2])
    X₂ = generaterv("Normal", "Moments", [20, 5])
    X = [X₁, X₂]
    ρˣ = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G₁(x::Vector) = x[1]^2 - 2 * x[2]
    G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

    # Perform the reliability analysis using MCFOSM:
    β₁ = MCFOSM(Problem₁)
    β₂ = MCFOSM(Problem₂)

    # Test the results:
    @test round(β₁, digits=2) == 1.66
    @test round(β₂, digits=2) == 4.29
end

@testset "Reliability Analysis: FORM" begin
    # Example from "Structural and System Reliability" by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Normal", "Moments", [10, 2])
    X₂ = generaterv("Normal", "Moments", [20, 5])
    X = [X₁, X₂]
    ρˣ = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G₁(x::Vector) = x[1]^2 - 2 * x[2]
    G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

    # Perform the reliability analysis using FORM:
    β₁, _, _ = FORM(Problem₁)
    β₂, _, _ = FORM(Problem₂)

    # Test the results:
    @test round(β₁, digits=2) == 2.11
    @test round(β₂, digits=2) == 2.11
end