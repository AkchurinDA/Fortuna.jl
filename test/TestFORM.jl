@testset "Reliability Analysis: FORM - MCFOSM" begin
    # Example 5.1 (p. 110) from "Structural and System Reliability" book by Armen Der Kiureghian

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
    @test isapprox(β₁, 1.66, rtol=0.01)
    @test isapprox(β₂, 4.29, rtol=0.01)
end

@testset "Reliability Analysis: FORM - HLRF #1" begin
    # Example 5.2 (p. 118) from "Structural and System Reliability" book by Armen Der Kiureghian

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
    β₁, _, x₁, u₁ = analyze(Problem₁, FORM(HLRF()))
    β₂, _, x₂, u₂ = analyze(Problem₂, FORM(HLRF()))

    # Test the results:
    @test isapprox(β₁, 2.11, rtol=0.01)
    @test isapprox(β₂, 2.11, rtol=0.01)
    @test isapprox(x₁[:, end], [6.14, 18.9], rtol=0.01)
    @test isapprox(x₂[:, end], [6.14, 18.9], rtol=0.01)
    @test isapprox(u₁[:, end], [-1.928, 0.852], rtol=0.01)
    @test isapprox(u₂[:, end], [-1.928, 0.852], rtol=0.01)
end

@testset "Reliability Analysis: FORM - HLRF #2" begin
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
    β, PoF, x, u = analyze(Problem, FORM(HLRF()))

    # Test the results:
    @test isapprox(β, 2.236067977499917, rtol=10^(-9))
    @test isapprox(PoF, 0.012673659338729965, rtol=10^(-9))
    @test isapprox(x[:, end], [160, 160], rtol=10^(-9))
    @test isapprox(u[:, end], [-2, 1], rtol=10^(-9))
end

@testset "Reliability Analysis: FORM - HLRF #3" begin
    # Example 6.5 (p. 147) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁ = generaterv("Normal", "Moments", [250, 250 * 0.3])
    M₂ = generaterv("Normal", "Moments", [125, 125 * 0.3])
    P = generaterv("Gumbel", "Moments", [2500, 2500 * 0.2])
    Y = generaterv("Weibull", "Moments", [40000, 40000 * 0.1])
    X = [M₁, M₂, P, Y]
    ρˣ = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    a = 0.190
    s₁ = 0.030
    s₂ = 0.015
    G(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using curve-fitting SORM:
    β, PoF, x, u = analyze(Problem, FORM(HLRF()))
    @test isapprox(β, 2.47, rtol=0.01)
    @test isapprox(PoF, 0.00682, rtol=0.01)
    @test isapprox(x[:, end], [341, 170, 3223, 31770], rtol=0.01)
    @test isapprox(u[:, end], [1.210, 0.699, 0.941, -1.80], rtol=0.01)
    # Note: There is a typo in the book for this example. The last coordinate of the design point in U-space must be -1.80.
end

@testset "Reliability Analysis: FORM - iHLRF #1" begin
    # Example 5.2 (p. 118) from "Structural and System Reliability" book by Armen Der Kiureghian

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
    β₁, _, x₁, u₁ = analyze(Problem₁, FORM(iHLRF()))
    β₂, _, x₂, u₂ = analyze(Problem₂, FORM(iHLRF()))

    # Test the results:
    @test isapprox(β₁, 2.11, rtol=0.01)
    @test isapprox(β₂, 2.11, rtol=0.01)
    @test isapprox(x₁[:, end], [6.14, 18.9], rtol=0.01)
    @test isapprox(x₂[:, end], [6.14, 18.9], rtol=0.01)
    @test isapprox(u₁[:, end], [-1.928, 0.852], rtol=0.01)
    @test isapprox(u₂[:, end], [-1.928, 0.852], rtol=0.01)
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
    β, PoF, x, u = analyze(Problem, FORM(iHLRF()))

    # Test the results:
    @test isapprox(β, 2.236067977499917, rtol=10^(-9))
    @test isapprox(PoF, 0.012673659338729965, rtol=10^(-9))
    @test isapprox(x[:, end], [160, 160], rtol=10^(-9))
    @test isapprox(u[:, end], [-2, 1], rtol=10^(-9))
end

@testset "Reliability Analysis: FORM - iHLRF #3" begin
    # Example 6.5 (p. 147) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁ = generaterv("Normal", "Moments", [250, 250 * 0.3])
    M₂ = generaterv("Normal", "Moments", [125, 125 * 0.3])
    P = generaterv("Gumbel", "Moments", [2500, 2500 * 0.2])
    Y = generaterv("Weibull", "Moments", [40000, 40000 * 0.1])
    X = [M₁, M₂, P, Y]
    ρˣ = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    a = 0.190
    s₁ = 0.030
    s₂ = 0.015
    G(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using curve-fitting SORM:
    β, PoF, x, u = analyze(Problem, FORM(iHLRF()))
    @test isapprox(β, 2.47, rtol=0.01)
    @test isapprox(PoF, 0.00682, rtol=0.01)
    @test isapprox(x[:, end], [341, 170, 3223, 31770], rtol=0.01)
    @test isapprox(u[:, end], [1.210, 0.699, 0.941, -1.80], rtol=0.01)
end