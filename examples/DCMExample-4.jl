using Fortuna

# Define the random vector:
X₁ = randomvariable("Normal", "P", [0, 1])
X₂ = randomvariable("Normal", "P", [0, 1])
X  = [X₁, X₂]

# Define the correlation matrix:
ρˣ = Matrix(1.0 * I, 2, 2)

# Define the limit state function:
g(x::AbstractVector) = 0.5 * (x[1] - 2) ^ 2 - 1.5 * (x[2] - 5) ^ 3 - 3

# Set up the counter:
const Counter = Ref(0)
CountFunctionCalls(F) = x -> begin
    Counter[] += 1
    F(x)
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, CountFunctionCalls(g))

# Solve the reliability problem using the DCM:
DCMSolution = solve(Problem, DCM(NumInitialDivisions = [12, 12], ϵ = 0.01))
println("Failure probability: ", DCMSolution.I)
println("Number of function calls: ", Counter[])
println("Error:", abs(DCMSolution.I - 2.85E-5) / 2.85E-5 * 100) 