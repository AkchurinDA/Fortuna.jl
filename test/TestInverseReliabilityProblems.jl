@testset "Inverse Reliability Problems" begin
    # Example 6.12 (p. 174) from "Structural and System Reliability" book by Armen Der Kiureghian

    # Define the random vector:
    X₁ = randomvariable("Normal", "M", [0, 1])
    X₂ = randomvariable("Normal", "M", [0, 1])
    X₃ = randomvariable("Normal", "M", [0, 1])
    X₄ = randomvariable("Normal", "M", [0, 1])
    X  = [X₁, X₂, X₃, X₄]

    # Define the correlation matrix:
    ρˣ = Matrix{Float64}(1.0 * I, 4, 4)

    # Define the limit state function:
    g(x::Vector, θ::Real) = exp(-θ * (x[1] + 2 * x[2] + 3 * x[3])) - x[4] + 1.5

    # Define the target reliability index:
    β = 2

    # Define an inverse reliability problem:
    Problem = InverseReliabilityProblem(X, ρˣ, g, β)

    # Perform the inverse reliability analysis:
    Solution = solve(Problem, 0.1, x₀ = [0.2, 0.2, 0.2, 0.2])

    # Test the results:
    @test isapprox(Solution.x[:, end], [+0.218, +0.436, +0.655, +1.826], atol = 1E-3)
    @test isapprox(Solution.θ[end], 0.367, atol = 1E-3)
end