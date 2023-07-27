# Mean-Centered First-Order Second-Moment Method:
function MCFOSM(Problem::ReliabilityProblem)
    # Extract data:
    g = Problem.g
    X = Problem.X
    ρˣ = Problem.ρˣ

    # Compute the means of marginal distrbutions:
    Mˣ = mean.(X)

    # Convert the correlation matrix into covariance matrix:
    σˣ = std.(X)
    Dˣ = diagm(σˣ)
    Σˣ = Dˣ * ρˣ * Dˣ

    # Compute gradient of the limit state function and evaluate it at the means of the marginal distributions:
    ∇g = gradient(g, Mˣ)

    # Compute the reliability index:
    β = g(Mˣ) / sqrt(transpose(∇g) * Σˣ * ∇g)

    return β
end