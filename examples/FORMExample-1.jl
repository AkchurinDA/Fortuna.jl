using Fortuna

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

# Perform the reliability analysis using FORM:
Solution = solve(Problem, FORM())
println("β   = $(Solution.β)  ")
println("PoF = $(Solution.PoF)")