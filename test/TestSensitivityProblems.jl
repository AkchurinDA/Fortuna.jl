@testset "Sensitivity Problems (Type I)" begin
    # Example 6.7 (p. 161) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define a random vector of correlated marginal distributions:
    M₁ = randomvariable("Normal",  "M", [250,   250   * 0.3])
    M₂ = randomvariable("Normal",  "M", [125,   125   * 0.3])
    P  = randomvariable("Gumbel",  "M", [2500,  2500  * 0.2])
    Y  = randomvariable("Weibull", "M", [40000, 40000 * 0.1])
    X  = [M₁, M₂, P, Y]
    ρˣ = [1 0.5 0.3 0; 0.5 1 0.3 0; 0.3 0.3 1 0; 0 0 0 1]

    # Define a limit state function:
    g(x::Vector, θ::Vector) = 1 - x[1] / (θ[1] * x[4]) - x[2] / (θ[2] * x[4]) - (x[3] / (θ[3] * x[4])) ^ 2

    # Define parameters of the limit state function:
    s₁ = 0.030
    s₂ = 0.015
    a  = 0.190
    Θ  = [s₁, s₂, a]

    # Define a sensitivity problem:
    Problem = SensitivityProblemTypeI(X, ρˣ, g, Θ)

    # Perform the sensitivity analysis:
    Solution = solve(Problem)

    # Test the results:
    @test isapprox(Solution.∇β,   [+36.80, +73.60, +9.260], rtol = 0.01)
    @test isapprox(Solution.∇PoF, [-0.700, -1.400, -0.176], rtol = 0.01)
end

@testset "Sensitivity Problems (Type II)" begin
    # Example 6.7 (p. 161) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define the random vector as a function of its moments:
    function XFunction(Θ::Vector)
        M₁ = randomvariable("Normal",  "M", [Θ[1], Θ[2]])
        M₂ = randomvariable("Normal",  "M", [Θ[3], Θ[4]])
        P  = randomvariable("Gumbel",  "P", [Θ[5], Θ[6]])
        Y  = randomvariable("Weibull", "P", [Θ[7], Θ[8]])

        return [M₁, M₂, P, Y]
    end

    # Define the correlation matrix:
    ρˣ = [
        1.0 0.5 0.3 0.0
        0.5 1.0 0.3 0.0
        0.3 0.3 1.0 0.0
        0.0 0.0 0.0 1.0]

    # Define the parameters of the random vector:
    Θ = [
        250,  250   * 0.30,
        125,  125   * 0.30,
        2257, 1 / 0.00257,
        12.2, 41700]

    # Define the limit state function:
    a            = 0.190
    s₁           = 0.030
    s₂           = 0.015
    g(x::Vector) = 1 - x[1] / (s₁ * x[4]) - x[2] / (s₂ * x[4]) - (x[3] / (a * x[4])) ^ 2

    # Define a sensitivity problem:
    Problem  = SensitivityProblemTypeII(XFunction, ρˣ, g, Θ)

    # Perform the sensitivity analysis:
    Solution = solve(Problem)

    # Test the results:
    @test isapprox(Solution.∇β,   1E-3 * [-3.240, -3.920, -6.480, -7.840, -0.546, -1.334, +88.8, +0.0951], rtol = 0.05)
    @test isapprox(Solution.∇PoF, 1E-4 * [+0.617, +0.746, +1.230, +1.490, +0.104, -0.254, -16.9, -0.0181], rtol = 0.05)
end