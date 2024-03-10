@testset "Sensitivity Problems" begin
    # Example 6.7 (p. 161) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁  = randomvariable("Normal", "M", [250, 250 * 0.3])
    M₂  = randomvariable("Normal", "M", [125, 125 * 0.3])
    P   = randomvariable("Gumbel", "M", [2500, 2500 * 0.2])
    Y   = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
    X   = [M₁, M₂, P, Y]
    ρˣ  = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    g(x, θ) = 1 - x[1] / (θ[1] * x[4]) - x[2] / (θ[2] * x[4]) - (x[3] / (θ[3] * x[4]))^2

    # Define parameters of the limit state function:
    s₁  = 0.030
    s₂  = 0.015
    a   = 0.190
    θ   = [s₁, s₂, a]

    # Define a reliability problem:
    Problem = SensitivityProblem(X, ρˣ, g, θ)

    # Perform the sensitivity analysis:
    Solution = solve(Problem)

    # Test the results:
    @test isapprox(Solution.∇β, [36.8, 73.6, 9.26], rtol = 0.01)
    @test isapprox(Solution.∇PoF, [-0.700, -1.400, -0.176], rtol = 0.01)
end