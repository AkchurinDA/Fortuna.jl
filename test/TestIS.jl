@testset "Importance Sampling #1" begin
    # Set an RNG seed:
    Random.seed!(123)

    # Define a list of reliability indices of interest:
    βList = 1:6

    for i in eachindex(βList)
        # Define a random vector of correlated marginal distributions:
        X₁ = randomvariable("Normal", "M", [0, 1])
        X₂ = randomvariable("Normal", "M", [0, 1])
        X  = [X₁, X₂]
        ρˣ = [1 0; 0 1]

        # Define a limit state function:
        g(x::Vector) = βList[i] * sqrt(2) - x[1] - x[2]

        # Define a reliability problem:
        Problem = ReliabilityProblem(X, ρˣ, g)

        # Perform the reliability analysis using Importance Sampling method:
        ProposalPDF = MvNormal([βList[i] / sqrt(2), βList[i] / sqrt(2)], [1 0; 0 1])
        Solution    = solve(Problem, IS(ProposalPDF, 10 ^ 6))

        # Test the results:
        @test isapprox(Solution.PoF, cdf(Normal(), -βList[i]), rtol = 5E-2)
    end
end