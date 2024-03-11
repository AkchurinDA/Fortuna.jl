@testset "Reliability Analysis: SSM #1" begin
    # Define a list of reliability indices of interest:
    βList = 1:6

    # Set an RNG seed:
    Random.seed!(123)

    for i in eachindex(βList)
        # Define a random vector of uncorrelated marginal distributions:
        X₁  = randomvariable("Normal", "M", [0, 1])
        X₂  = randomvariable("Normal", "M", [0, 1])
        X   = [X₁, X₂]
        ρˣ  = [1 0; 0 1]

        # Define a limit state function:
        g(x::Vector) = βList[i] * sqrt(2) - x[1] - x[2]

        # Define reliability problems:
        Problem = ReliabilityProblem(X, ρˣ, g)

        # Perform the reliability analysis using SSM:
        Solution = solve(Problem, SSM())

        # Test the results:
        @test isapprox(Solution.PoF, cdf(Normal(), -βList[i]), rtol = 0.05)
    end
end