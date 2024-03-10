@testset "Monte Carlo #1" begin
    # Define a list of reliability indices of interest:
    βList = 1:3

    for i in eachindex(βList)
        # Define a random vector of correlated marginal distributions:
        X₁  = randomvariable("Normal", "M", [0, 1])
        X₂  = randomvariable("Normal", "M", [0, 1])
        X   = [X₁, X₂]
        ρˣ  = [1 0; 0 1]

        # Define a limit state function:
        g(x::Vector) = βList[i] * sqrt(2) - x[1] - x[2]

        # Define a reliability problem:
        Problem = ReliabilityProblem(X, ρˣ, g)

        # Perform the reliability analysis using Importance Sampling method:
        Solution = solve(Problem, MC())

        # Test the results:
        @test isapprox(Solution.PoF, cdf(Normal(), -βList[i]), rtol = 0.05)
    end
end

@testset "Monte Carlo #2" begin
    # Define a random vector of correlated marginal distributions:
    M₁  = randomvariable("Normal", "M", [250, 250 * 0.3])
    M₂  = randomvariable("Normal", "M", [125, 125 * 0.3])
    P   = randomvariable("Gumbel", "M", [2500, 2500 * 0.2])
    Y   = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
    X   = [M₁, M₂, P, Y]
    ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    a   = 0.190
    s₁  = 0.030
    s₂  = 0.015
    g(x) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using Monte Carlo simulations:
    Solution = solve(Problem, MC())

    # Test the results:
    @test isapprox(Solution.PoF, 0.00931, rtol = 0.05)
end