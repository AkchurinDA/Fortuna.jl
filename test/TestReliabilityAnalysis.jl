@testset "Reliability Analysis: FORM - MCFOSM" begin
    # Example from "Structural and System Reliability" book by Armen Der Kiureghian

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
    β₁ = analyze(Problem₁, FORM(MCFOSM()))
    β₂ = analyze(Problem₂, FORM(MCFOSM()))

    # Test the results:
    @test round(β₁, digits=2) == 1.66
    @test round(β₂, digits=2) == 4.29
end

@testset "Reliability Analysis: FORM - HLRF #1" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Normal", "Parameters", [200, 20])
    X₂ = generaterv("Normal", "Parameters", [150, 10])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G(x::Vector) = x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using FORM:
    β, x, u = analyze(Problem, FORM(HLRF()))

    # Test the results:
    @test isapprox(β, 2.236067977499917, rtol=10^(-9))
    @test isapprox(x[:, end], [160, 160], rtol=10^(-9))
    @test isapprox(u[:, end], [-2, 1], rtol=10^(-9))
end

@testset "Reliability Analysis: FORM - HLRF #2" begin
    # Example from "Structural and System Reliability" book by Armen Der Kiureghian

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
    β₁, _, _ = analyze(Problem₁, FORM(HLRF()))
    β₂, _, _ = analyze(Problem₂, FORM(HLRF()))

    # Test the results:
    @test round(β₁, digits=2) == 2.11
    @test round(β₂, digits=2) == 2.11
end

@testset "Reliability Analysis: FORM - iHLRF #1" begin
    # Example from "Structural and System Reliability" book by Armen Der Kiureghian

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
    β₁, _, _ = analyze(Problem₁, FORM(iHLRF()))
    β₂, _, _ = analyze(Problem₂, FORM(iHLRF()))

    # Test the results:
    @test round(β₁, digits=2) == 2.11
    @test round(β₂, digits=2) == 2.11
end

@testset "Reliability Analysis: FORM - iHLRF #2" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = generaterv("Normal", "Parameters", [200, 20])
    X₂ = generaterv("Normal", "Parameters", [150, 10])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G(x::Vector) = x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using FORM:
    β, x, u = analyze(Problem, FORM(iHLRF()))

    # Test the results:
    @test isapprox(β, 2.236067977499917, rtol=10^(-9))
    @test isapprox(x[:, end], [160, 160], rtol=10^(-9))
    @test isapprox(u[:, end], [-2, 1], rtol=10^(-9))
end