using Fortuna
using Distributions

# Define the problem dimensionality:
N = 2

# Define the random vector:
X  = [randomvariable("Normal", "M", [0, 1]) for _ in 1:N]

# Define the correlation matrix:
ρˣ = Matrix(1.0 * I, N, N)

# Define the limit state function:
β = 6
g(x::AbstractVector) = β * sqrt(N) - sum(x)

# Set up the counter:
const Counter = Ref(0)
CountFunctionCalls(F) = x -> begin
    Counter[] += 1
    F(x)
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, CountFunctionCalls(g))

# Solve the reliability problem using the DCM:
DCMSolution = solve(Problem, DCM(NumInitialDivisions = vcat(4, 12), ϵ = 0.10), RL = β)
println("Failure probability: ", DCMSolution.I)
println("Number of function calls: ", Counter[])
println("Error:", abs(DCMSolution.I - cdf(Normal(), -β)) / cdf(Normal(), -β) * 100)