@testset "Reliability Analysis: SSM #1" begin
    # Define a random vector of uncorrelated marginal distributions:
    X₁ = generaterv("Normal", "M", [0, 1])
    X₂ = generaterv("Normal", "M", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define a limit state function:
    g(x::Vector) = 3.5 - (x[1] + x[2]) / sqrt(2)

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using SSM:
    Solution = analyze(Problem, SSM())

    # Test the results:
    @test isapprox(Solution.PoF, 0.000232629, rtol=0.10)
end

@testset "Reliability Analysis: SSM #2" begin
    # Define a random vector of uncorrelated marginal distributions:
    X₁ = generaterv("Exponential", "P", 1)
    X₂ = generaterv("Exponential", "P", 1)
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define a limit state function:
    g(x::Vector) = 10 - x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using SSM:
    Solution = analyze(Problem, SSM())

    # Test the results:
    @test isapprox(Solution.PoF, 0.000499399, rtol=0.10)
end

@testset "Reliability Analysis: SSM #3" begin
    # Define a random vector of uncorrelated marginal distributions:
    X₁ = generaterv("Normal", "M", [0, 1])
    X₂ = generaterv("Normal", "M", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define a limit state function:
    g(x::Vector) = minimum([3.2 + (x[1] + x[2]) / sqrt(2), 0.1 * (x[1] - x[2])^2 - (x[1] + x[2]) / sqrt(2) + 2.5])

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using SSM:
    Solution = analyze(Problem, SSM())

    # Test the results:
    @test isapprox(Solution.PoF, 0.0049, rtol=0.10)
end

@testset "Reliability Analysis: SSM #4" begin
    # Define a random vector of uncorrelated marginal distributions:
    X₁ = generaterv("Normal", "M", [0, 1])
    X₂ = generaterv("Normal", "M", [0, 1])
    X = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define a limit state function:
    g(x::Vector) = 5 - x[2] - 0.5 * (x[1] - 0.1)^2

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using SSM:
    Solution = analyze(Problem, SSM())

    # Test the results:
    @test isapprox(Solution.PoF, 0.00301, rtol=0.10)
end