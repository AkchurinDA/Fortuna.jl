@testset "Reliability Analysis: FORM - MCFOSM" begin
    # Example 5.1 (p. 110) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁  = generaterv("Normal", "M", [10, 2])
    X₂  = generaterv("Normal", "M", [20, 5])
    X   = [X₁, X₂]
    ρˣ  = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G₁(x::Vector) = x[1]^2 - 2 * x[2]
    G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

    # Perform the reliability analysis using MCFOSM:
    Solution₁ = analyze(Problem₁, FORM(MCFOSM()))
    Solution₂ = analyze(Problem₂, FORM(MCFOSM()))

    # Test the results:
    @test isapprox(Solution₁.β, 1.66, rtol = 0.01)
    @test isapprox(Solution₂.β, 4.29, rtol = 0.01)
end

@testset "Reliability Analysis: FORM - HLRF #1" begin
    # Example 5.2 (p. 118) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁  = generaterv("Normal", "M", [10, 2])
    X₂  = generaterv("Normal", "M", [20, 5])
    X   = [X₁, X₂]
    ρˣ  = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G₁(x::Vector) = x[1]^2 - 2 * x[2]
    G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

    # Perform the reliability analysis using FORM:
    Solution₁ = analyze(Problem₁, FORM(HLRF()))
    Solution₂ = analyze(Problem₂, FORM(HLRF()))

    # Test the results:
    @test isapprox(Solution₁.β, 2.11, rtol = 0.01)
    @test isapprox(Solution₂.β, 2.11, rtol = 0.01)
    @test isapprox(Solution₁.x[:, end], [6.14, 18.9], rtol = 0.01)
    @test isapprox(Solution₂.x[:, end], [6.14, 18.9], rtol = 0.01)
    @test isapprox(Solution₁.u[:, end], [-1.928, 0.852], rtol = 0.01)
    @test isapprox(Solution₂.u[:, end], [-1.928, 0.852], rtol = 0.01)
end

@testset "Reliability Analysis: FORM - HLRF #2" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁  = generaterv("Normal", "Parameters", [200, 20])
    X₂  = generaterv("Normal", "Parameters", [150, 10])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G(x::Vector) = x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using FORM:
    Solution = analyze(Problem, FORM(HLRF()))

    # Test the results:
    @test isapprox(Solution.β, 2.236067977499917, rtol = 10^(-9))
    @test isapprox(Solution.PoF, 0.012673659338729965, rtol = 10^(-9))
    @test isapprox(Solution.x[:, end], [160, 160], rtol = 10^(-9))
    @test isapprox(Solution.u[:, end], [-2, 1], rtol = 10^(-9))
end

@testset "Reliability Analysis: FORM - HLRF #3" begin
    # Example 6.5 (p. 147) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁  = generaterv("Normal", "M", [250, 250 * 0.3])
    M₂  = generaterv("Normal", "M", [125, 125 * 0.3])
    P   = generaterv("Gumbel", "M", [2500, 2500 * 0.2])
    Y   = generaterv("Weibull", "M", [40000, 40000 * 0.1])
    X   = [M₁, M₂, P, Y]
    ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    a   = 0.190
    s₁  = 0.030
    s₂  = 0.015
    G(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using curve-fitting SORM:
    Solution = analyze(Problem, FORM(HLRF()))
    @test isapprox(Solution.β, 2.47, rtol = 0.01)
    @test isapprox(Solution.PoF, 0.00682, rtol = 0.01)
    @test isapprox(Solution.x[:, end], [341, 170, 3223, 31770], rtol = 0.01)
    @test isapprox(Solution.u[:, end], [1.210, 0.699, 0.941, -1.80], rtol = 0.01)
    # Note: There is a typo in the book for this example. The last coordinate of the design point in U-space must be -1.80.
end

@testset "Reliability Analysis: FORM - iHLRF #1" begin
    # Example 5.2 (p. 118) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁  = generaterv("Normal", "M", [10, 2])
    X₂  = generaterv("Normal", "M", [20, 5])
    X   = [X₁, X₂]
    ρˣ  = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G₁(x::Vector) = x[1]^2 - 2 * x[2]
    G₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, G₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, G₂)

    # Perform the reliability analysis using FORM:
    Solution₁ = analyze(Problem₁, FORM(iHLRF()))
    Solution₂ = analyze(Problem₂, FORM(iHLRF()))

    # Test the results:
    @test isapprox(Solution₁.β, 2.11, rtol = 0.01)
    @test isapprox(Solution₂.β, 2.11, rtol = 0.01)
    @test isapprox(Solution₁.x[:, end], [6.14, 18.9], rtol = 0.01)
    @test isapprox(Solution₂.x[:, end], [6.14, 18.9], rtol = 0.01)
    @test isapprox(Solution₁.u[:, end], [-1.928, 0.852], rtol = 0.01)
    @test isapprox(Solution₂.u[:, end], [-1.928, 0.852], rtol = 0.01)
end

@testset "Reliability Analysis: FORM - iHLRF #2" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁  = generaterv("Normal", "Parameters", [200, 20])
    X₂  = generaterv("Normal", "Parameters", [150, 10])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    G(x::Vector) = x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using FORM:
    Solution = analyze(Problem, FORM(iHLRF()))

    # Test the results:
    @test isapprox(Solution.β, 2.236067977499917, rtol = 10^(-9))
    @test isapprox(Solution.PoF, 0.012673659338729965, rtol = 10^(-9))
    @test isapprox(Solution.x[:, end], [160, 160], rtol = 10^(-9))
    @test isapprox(Solution.u[:, end], [-2, 1], rtol = 10^(-9))
end

@testset "Reliability Analysis: FORM - iHLRF #3" begin
    # Example 6.5 (p. 147) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁  = generaterv("Normal", "M", [250, 250 * 0.3])
    M₂  = generaterv("Normal", "M", [125, 125 * 0.3])
    P   = generaterv("Gumbel", "M", [2500, 2500 * 0.2])
    Y   = generaterv("Weibull", "M", [40000, 40000 * 0.1])
    X   = [M₁, M₂, P, Y]
    ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    a   = 0.190
    s₁  = 0.030
    s₂  = 0.015
    G(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, G)

    # Perform the reliability analysis using curve-fitting SORM:
    Solution = analyze(Problem, FORM(iHLRF()))
    @test isapprox(Solution.β, 2.47, rtol = 0.01)
    @test isapprox(Solution.PoF, 0.00682, rtol = 0.01)
    @test isapprox(Solution.x[:, end], [341, 170, 3223, 31770], rtol = 0.01)
    @test isapprox(Solution.u[:, end], [1.210, 0.699, 0.941, -1.80], rtol = 0.01)
end