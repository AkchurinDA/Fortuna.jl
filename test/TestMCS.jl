@testset "Reliability Problems: MCS" begin
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

    # Perform the reliability analysis using Monte Carlo simulations:
    Solution = analyze(Problem, MCS())
    @test isapprox(Solution.PoF, 0.00931, rtol = 0.05)
end