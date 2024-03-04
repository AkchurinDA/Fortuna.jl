@testset "Reliability Problems: IS" begin
    # Define a random vector of correlated marginal distributions:
    X₁  = generaterv("Normal", "M", [0, 1])
    X₂  = generaterv("Normal", "M", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Define a limit state function:
    g(x::Vector) = 4 * sqrt(2) - x[1] - x[2]

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using Monte Carlo simulations:
    Solution = analyze(Problem, IS(MvNormal([4 / sqrt(2), 4 / sqrt(2)], [2.5 0; 0 2.5]), 10 ^ 6))
    @test isapprox(Solution.PoF, cdf(Normal(0, 1), -4), rtol = 0.05)
end