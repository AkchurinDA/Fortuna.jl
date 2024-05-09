@testset "SORM #1 - CF" begin
    # Examples 6.5 (p. 147) and 7.2 (p. 188) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁ = randomvariable("Normal",  "M", [250,   250   * 0.3])
    M₂ = randomvariable("Normal",  "M", [125,   125   * 0.3])
    P  = randomvariable("Gumbel",  "M", [2500,  2500  * 0.2])
    Y  = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
    X  = [M₁, M₂, P, Y]
    ρˣ = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    a            = 0.190
    s₁           = 0.030
    s₂           = 0.015
    g(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using curve-fitting SORM:
    Solution = solve(Problem, SORM(CF()))

    # Test the results:
    @test isapprox(Solution.β₂,   [2.35, 2.35],         rtol = 1E-2)
    @test isapprox(Solution.PoF₂, [0.00960, 0.00914],   rtol = 1E-2)
    @test isapprox(Solution.κ,    [-0.155, -0.0399, 0], rtol = 1E-2)
end

@testset "SORM #2 - CF" begin
    # https://www.researchgate.net/publication/370230768_Structural_reliability_analysis_by_line_sampling_A_Bayesian_active_learning_treatment

    # Define random vector:
    m   = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.05])
    k₁  = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.10])
    k₂  = randomvariable("LogNormal", "M", [0.2, 0.2 * 0.10])
    r   = randomvariable("LogNormal", "M", [0.5, 0.5 * 0.10])
    F₁  = randomvariable("LogNormal", "M", [0.4, 0.4 * 0.20])
    t₁  = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.20])
    X   = [m, k₁, k₂, r, F₁, t₁]

    # Define correlation matrix:
    ρˣ  = Matrix(1.0 * I, 6, 6)

    # Define limit state function:
    g(x::Vector) = 3 * x[4] - abs(((2 * x[5]) / (x[2] + x[3])) * sin((x[6] / 2) * sqrt((x[2] + x[3]) / x[1])))

    # Define reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Solve reliability problem using Curve-Fitting Method:
    Solution = solve(Problem, SORM(CF()))

    # Test the results:
    @test isapprox(Solution.PoF₂[1], 4.08 * 10 ^ (-8), rtol = 5E-2)
    @test isapprox(Solution.PoF₂[2], 4.08 * 10 ^ (-8), rtol = 5E-2)
end

@testset "SORM #3 - PF" begin
    # Examples 6.5 (p. 147) and 7.7 (p. 196) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁ = randomvariable("Normal", "M",  [250,   250   * 0.3])
    M₂ = randomvariable("Normal", "M",  [125,   125   * 0.3])
    P  = randomvariable("Gumbel", "M",  [2500,  2500  * 0.2])
    Y  = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
    X  = [M₁, M₂, P, Y]
    ρˣ = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    a            = 0.190
    s₁           = 0.030
    s₂           = 0.015
    g(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4]))^2

    # Define a reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using point-fitting SORM:
    Solution = solve(Problem, SORM(PF()))

    # Test the results:
    @test isapprox(Solution.β₂,                   [2.36, 2.36],       rtol = 5E-2)
    @test isapprox(Solution.PoF₂,                 [0.00913, 0.00913], rtol = 5E-2)
    @test isapprox(Solution.FittingPoints⁻[1, :], [-2.47, +2.27],     rtol = 5E-2)
    @test isapprox(Solution.FittingPoints⁻[2, :], [-2.47, +2.43],     rtol = 5E-2)
    @test isapprox(Solution.FittingPoints⁻[3, :], [-2.47, +2.05],     rtol = 5E-2)
    @test isapprox(Solution.FittingPoints⁺[1, :], [+2.47, +2.34],     rtol = 5E-2)
    @test isapprox(Solution.FittingPoints⁺[2, :], [+2.47, +2.44],     rtol = 5E-2)
    @test isapprox(Solution.FittingPoints⁺[3, :], [+2.47, +2.13],     rtol = 5E-2)
    @test isapprox(Solution.κ₁[1, :],             [-0.0630, -0.0405], rtol = 5E-2)
    @test isapprox(Solution.κ₁[2, :],             [-0.0097, -0.0120], rtol = 5E-2)
    @test isapprox(Solution.κ₁[3, :],             [-0.1380, -0.1110], rtol = 5E-2)
end

@testset "SORM #4 - PF" begin
    # https://www.researchgate.net/publication/370230768_Structural_reliability_analysis_by_line_sampling_A_Bayesian_active_learning_treatment

    # Define random vector:
    m  = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.05])
    k₁ = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.10])
    k₂ = randomvariable("LogNormal", "M", [0.2, 0.2 * 0.10])
    r  = randomvariable("LogNormal", "M", [0.5, 0.5 * 0.10])
    F₁ = randomvariable("LogNormal", "M", [0.4, 0.4 * 0.20])
    t₁ = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.20])
    X  = [m, k₁, k₂, r, F₁, t₁]

    # Define correlation matrix:
    ρˣ  = Matrix(1.0 * I, 6, 6)

    # Define limit state function:
    g(x::Vector) = 3 * x[4] - abs(((2 * x[5]) / (x[2] + x[3])) * sin((x[6] / 2) * sqrt((x[2] + x[3]) / x[1])))

    # Define reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Solve reliability problem using Point-Fitting Method:
    Solution = solve(Problem, SORM(PF()))
    
    # Test the results:
    @test isapprox(Solution.PoF₂[1], 4.08 * 10 ^ (-8), rtol = 5E-2)
    @test isapprox(Solution.PoF₂[2], 4.08 * 10 ^ (-8), rtol = 5E-2)
end