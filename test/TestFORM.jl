@testset "FORM #1 - MCFOSM" begin
    # Example 5.1 (p. 110) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁ = randomvariable("Normal", "M", [10, 2])
    X₂ = randomvariable("Normal", "M", [20, 5])
    X  = [X₁, X₂]
    ρˣ = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    g₁(x::Vector) = x[1]^2 - 2 * x[2]
    g₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

    # Perform the reliability analysis using MCFOSM:
    Solution₁ = solve(Problem₁, FORM(MCFOSM()))
    Solution₂ = solve(Problem₂, FORM(MCFOSM()))

    # Test the results:
    @test isapprox(Solution₁.β, 1.66, rtol = 1E-2)
    @test isapprox(Solution₂.β, 4.29, rtol = 1E-2)
end

@testset "FORM #2 - RF" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = randomvariable("Normal", "M", [200, 20])
    X₂ = randomvariable("Normal", "M", [150, 10])
    X  = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    g(x::Vector) = x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using iHLRF:
    Solution = solve(Problem, FORM(RF()))

    # Test the results:
    @test isapprox(Solution.β,         2.236067977499917, rtol = 1E-9)
    @test isapprox(Solution.x[:, end], [160, 160],        rtol = 1E-9)
    @test isapprox(Solution.u[:, end], [-2, 1],           rtol = 1E-9)
end

@testset "FORM #3 - HLRF" begin
    # Example 5.2 (p. 118) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁ = randomvariable("Normal", "M", [10, 2])
    X₂ = randomvariable("Normal", "M", [20, 5])
    X  = [X₁, X₂]
    ρˣ = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    g₁(x::Vector) = x[1]^2 - 2 * x[2]
    g₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

    # Perform the reliability analysis using HLRF:
    Solution₁ = solve(Problem₁, FORM(HLRF()))
    Solution₂ = solve(Problem₂, FORM(HLRF()))

    # Test the results:
    @test isapprox(Solution₁.β,         2.11,            rtol = 1E-2)
    @test isapprox(Solution₂.β,         2.11,            rtol = 1E-2)
    @test isapprox(Solution₁.x[:, end], [6.14, 18.9],    rtol = 1E-2)
    @test isapprox(Solution₂.x[:, end], [6.14, 18.9],    rtol = 1E-2)
    @test isapprox(Solution₁.u[:, end], [-1.928, 0.852], rtol = 1E-2)
    @test isapprox(Solution₂.u[:, end], [-1.928, 0.852], rtol = 1E-2)
end

@testset "FORM #4 - HLRF" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = randomvariable("Normal", "M", [200, 20])
    X₂ = randomvariable("Normal", "M", [150, 10])
    X  = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    g(x::Vector) = x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using HLRF:
    Solution = solve(Problem, FORM(HLRF()))

    # Test the results:
    @test isapprox(Solution.β,         2.236067977499917,    rtol = 1E-9)
    @test isapprox(Solution.PoF,       0.012673659338729965, rtol = 1E-9)
    @test isapprox(Solution.x[:, end], [160, 160],           rtol = 1E-9)
    @test isapprox(Solution.u[:, end], [-2, 1],              rtol = 1E-9)
end

@testset "FORM #5 - HLRF" begin
    # Example 6.5 (p. 147) from "Structural and System Reliability" book by Armen Der Kiureghian

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

    # Perform the reliability analysis using HLRF:
    Solution = solve(Problem, FORM(HLRF()))

    # Test the results:
    @test isapprox(Solution.β,         2.47,                          rtol = 1E-2)
    @test isapprox(Solution.PoF,       0.00682,                       rtol = 1E-2)
    @test isapprox(Solution.x[:, end], [341, 170, 3223, 31770],       rtol = 1E-2)
    @test isapprox(Solution.u[:, end], [1.210, 0.699, 0.941, -1.80],  rtol = 1E-2)
    @test isapprox(Solution.γ,         [0.269, 0.269, 0.451, -0.808], rtol = 1E-2)
    # Note: There is a typo in the book for this example. The last coordinate of the design point in U-space must be -1.80.
end

@testset "FORM #6 - HLRF" begin
    # https://www.researchgate.net/publication/370230768_Structural_reliability_analysis_by_line_sampling_A_Bayesian_active_learning_treatment

    # Define random vector:
    m  = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.05])
    k₁ = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.10])
    k₂ = randomvariable("LogNormal", "M", [0.2, 0.2 * 0.10])
    r  = randomvariable("LogNormal", "M", [0.5, 0.5 * 0.10])
    F₁ = randomvariable("LogNormal", "M", [0.4, 0.4 * 0.20])
    t₁ = randomvariable("LogNormal", "M", [1.0, 1.0 * 0.20])
    X  = [m, k₁, k₂, r, F₁, t₁]
    ρˣ = Matrix(1.0 * I, 6, 6)

    # Define limit state function:
    g(x::Vector) = 3 * x[4] - abs(((2 * x[5]) / (x[2] + x[3])) * sin((x[6] / 2) * sqrt((x[2] + x[3]) / x[1])))

    # Define reliability problem:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using HLRF:
    Solution = solve(Problem, FORM(HLRF()))
    
    # Test the results:
    @test isapprox(Solution.PoF, 4.88 * 10 ^ (-8), rtol = 1E-2)
end

@testset "FORM #7 - iHLRF" begin
    # Example 5.2 (p. 118) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    X₁ = randomvariable("Normal", "M", [10, 2])
    X₂ = randomvariable("Normal", "M", [20, 5])
    X  = [X₁, X₂]
    ρˣ = [1 0.5; 0.5 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    g₁(x::Vector) = x[1]^2 - 2 * x[2]
    g₂(x::Vector) = 1 - 2 * x[2] / x[1]^2

    # Define reliability problems:
    Problem₁ = ReliabilityProblem(X, ρˣ, g₁)
    Problem₂ = ReliabilityProblem(X, ρˣ, g₂)

    # Perform the reliability analysis using iHLRF:
    Solution₁ = solve(Problem₁, FORM(iHLRF()))
    Solution₂ = solve(Problem₂, FORM(iHLRF()))

    # Test the results:
    @test isapprox(Solution₁.β,         2.11,            rtol = 1E-2)
    @test isapprox(Solution₂.β,         2.11,            rtol = 1E-2)
    @test isapprox(Solution₁.x[:, end], [6.14, 18.9],    rtol = 1E-2)
    @test isapprox(Solution₂.x[:, end], [6.14, 18.9],    rtol = 1E-2)
    @test isapprox(Solution₁.u[:, end], [-1.928, 0.852], rtol = 1E-2)
    @test isapprox(Solution₂.u[:, end], [-1.928, 0.852], rtol = 1E-2)
end

@testset "FORM #8 - iHLRF" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector of correlated marginal distributions:
    X₁ = randomvariable("Normal", "M", [200, 20])
    X₂ = randomvariable("Normal", "M", [150, 10])
    X  = [X₁, X₂]
    ρˣ = [1 0; 0 1]

    # Define two equivalent limit state functions to demonstrate the invariance problem of the MCFOSM method:
    g(x::Vector) = x[1] - x[2]

    # Define reliability problems:
    Problem = ReliabilityProblem(X, ρˣ, g)

    # Perform the reliability analysis using iHLRF:
    Solution = solve(Problem, FORM(iHLRF()))

    # Test the results:
    @test isapprox(Solution.β,         2.236067977499917,    rtol = 1E-9)
    @test isapprox(Solution.PoF,       0.012673659338729965, rtol = 1E-9)
    @test isapprox(Solution.x[:, end], [160, 160],           rtol = 1E-9)
    @test isapprox(Solution.u[:, end], [-2, 1],              rtol = 1E-9)
end

@testset "FORM #9 - iHLRF" begin
    # Example 6.5 (p. 147) from "Structural and System Reliability" book by Armen Der Kiureghian

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

    # Perform the reliability analysis using iHLRF:
    Solution = solve(Problem, FORM(iHLRF()))

    # Test the results:
    @test isapprox(Solution.β,         2.47,                          rtol = 1E-2)
    @test isapprox(Solution.PoF,       0.00682,                       rtol = 1E-2)
    @test isapprox(Solution.x[:, end], [341, 170, 3223, 31770],       rtol = 1E-2)
    @test isapprox(Solution.u[:, end], [1.210, 0.699, 0.941, -1.80],  rtol = 1E-2)
    @test isapprox(Solution.γ,         [0.269, 0.269, 0.451, -0.808], rtol = 1E-2)
    # Note: There is a typo in the book for this example. The last coordinate of the design point in U-space must be -1.80.
end

@testset "FORM #10 - iHLRF" begin
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

    # Perform the reliability analysis using iHLRF:
    Solution = solve(Problem, FORM(iHLRF()))

    # Test the results:
    @test isapprox(Solution.PoF, 4.88 * 10 ^ (-8), rtol = 1E-2)
end