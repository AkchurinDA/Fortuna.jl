@testset "Reliability Analysis: SSM #1" begin
    # Define a list of reliability indices of interest:
    βList = 1:6

    # Set an RNG seed:
    Random.seed!(123)

    for i in eachindex(βList)
        # Define a random vector of uncorrelated marginal distributions:
        X₁  = randomvariable("Normal", "M", [0, 1])
        X₂  = randomvariable("Normal", "M", [0, 1])
        X   = [X₁, X₂]
        ρˣ  = [1 0; 0 1]

        # Define a limit state function:
        g(x::Vector) = βList[i] * sqrt(2) - x[1] - x[2]

        # Define reliability problems:
        Problem = ReliabilityProblem(X, ρˣ, g)

        # Perform the reliability analysis using SSM:
        Solution = solve(Problem, SSM())

        # Test the results:
        @test isapprox(Solution.PoF, cdf(Normal(), -βList[i]), rtol = 0.05)
    end
end

@testset "Reliability Analysis: SSM #2" begin
    # https://www.researchgate.net/publication/370230768_Structural_reliability_analysis_by_line_sampling_A_Bayesian_active_learning_treatment

    # Define random vector:
    X₁  = randomvariable("Normal", "M", [0, 1])
    X₂  = randomvariable("Normal", "M", [0, 1])
    X   = [X₁, X₂]

    # Define correlation matrix:
    ρˣ  = [1 0; 0 1]

    # Define limit state function:
    a = 5.50
    b = 0.02
    c = 5 / 6
    d = π / 3
    g(x::Vector) = a - x[2] + b * x[1] ^ 3 + c * sin(d * x[1])

    # Define reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Solve reliability problem using Subset Simulation Method:
    Solution = solve(Problem, SSM())

    # Test the results:
    @test isapprox(Solution.PoF, 3.53 * 10 ^ (-7), rtol = 0.05)
end