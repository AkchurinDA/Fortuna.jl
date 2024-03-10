@testset "Reliability Analysis: SORM - CF" begin
    # Examples 6.5 (p. 147) and 7.2 (p. 188) from "Structural and System Reliability" book by Armen Der Kiureghian

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
    g(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using curve-fitting SORM:
    Solution = solve(Problem, SORM(CF()))

    # Test the results:
    @test isapprox(Solution.β₂, [2.35, 2.35], rtol = 0.01)
    @test isapprox(Solution.PoF₂, [0.00960, 0.00914], rtol = 0.01)
    @test isapprox(Solution.κ, [-0.155, -0.0399, 0], rtol = 0.01)
end

@testset "Reliability Analysis: SORM - PF" begin
    # Examples 6.5 (p. 147) and 7.7 (p. 196) from "Structural and System Reliability" book by Armen Der Kiureghian

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
    g(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using point-fitting SORM:
    Solution = solve(Problem, SORM(PF()))

    # Test the results:
    @test isapprox(Solution.β₂, [2.36, 2.36], rtol = 0.01)
    @test isapprox(Solution.PoF₂, [0.00913, 0.00913], rtol = 0.05)
    @test isapprox(Solution.FittingPoints⁻[1, :], [-2.47, +2.27], rtol = 0.05)
    @test isapprox(Solution.FittingPoints⁻[2, :], [-2.47, +2.43], rtol = 0.05)
    @test isapprox(Solution.FittingPoints⁻[3, :], [-2.47, +2.05], rtol = 0.05)
    @test isapprox(Solution.FittingPoints⁺[1, :], [+2.47, +2.34], rtol = 0.05)
    @test isapprox(Solution.FittingPoints⁺[2, :], [+2.47, +2.44], rtol = 0.05)
    @test isapprox(Solution.FittingPoints⁺[3, :], [+2.47, +2.13], rtol = 0.05)
    @test isapprox(Solution.κ₁[1, :], [-0.0630, -0.0405], rtol = 0.01)
    @test isapprox(Solution.κ₁[2, :], [-0.0097, -0.0120], rtol = 0.01)
    @test isapprox(Solution.κ₁[3, :], [-0.1380, -0.1110], rtol = 0.01)
end