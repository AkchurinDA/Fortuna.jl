using Fortuna

# Define the random vector:
X₁ = randomvariable("Normal", "M", [0, 1])
X₂ = randomvariable("Normal", "M", [0, 1])
X  = [X₁, X₂]

# Define the correlation matrix:
ρˣ = Matrix(1.0 * I, 2, 2)

# Define the limit state function:
g(x::AbstractVector) = 10 - (x[1] ^ 2 - 5 * cospi(2 * x[1]) + x[2] ^ 2 - 5 * cospi(2 * x[2]))

# Set up the counter:
const Counter = Ref(0)
CountFunctionCalls(F) = x -> begin
    Counter[] += 1
    F(x)
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, CountFunctionCalls(g))

# Solve the reliability problem using the DCM:
DCMSolution = solve(Problem, DCM(NumInitialDivisions = [4, 12], ϵ = 0.01), RU = 6)
println("Failure probability: ", DCMSolution.I)
println("Number of function calls: ", Counter[])
println("Error:", abs(DCMSolution.I - 7.31E-2) / 7.31E-2)