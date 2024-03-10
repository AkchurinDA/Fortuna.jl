@testset "Reliability Problems: IS" begin
    # Define a random vector of correlated marginal distributions:
    X₁  = randomvariable("Normal", "M", [0, 1])
    X₂  = randomvariable("Normal", "M", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Define a limit state function:
    g(x::Vector) = 4 * sqrt(2) - x[1] - x[2]

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using Importance Sampling method:
    q = MvNormal([4 / sqrt(2), 4 / sqrt(2)], [2.5 0; 0 2.5])
    Solution = solve(Problem, IS(q, 10 ^ 6))

    # Test the results:
    PoF = Solution.PoF
    @test isapprox(PoF, cdf(Normal(), -4), rtol = 0.01)
end