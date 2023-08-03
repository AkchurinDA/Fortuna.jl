@testset "Reliability Analysis: SORM - CF" begin
    # Examples 6.5 (p. 147) and 7.1 (p. 188) from "Structural and System Reliability" book by Armen Der Kiureghian

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
    Solution = analyze(Problem, SORM(CF()))

    # Test the results:
    @test isapprox(Solution.β₂, [2.35, 2.35], rtol=0.01)
    @test isapprox(Solution.PoF₂, [0.00960, 0.00914], rtol=0.01)
    @test isapprox(Solution.κ, [-0.155, -0.0399, 0], rtol=0.01)
end