using Fortuna

# Define random vector:
X₁  = randomvariable("Normal",      "M",    [1260,  1260 * 0.05])
X₂  = randomvariable("LogNormal",   "M",    [300,    300 * 0.10])
X₃  = randomvariable("Normal",      "M",    [770,    770 * 0.05])
X₄  = randomvariable("LogNormal",   "M",    [0.35,  0.35 * 0.10])
X₅  = randomvariable("LogNormal",   "M",    [30,      30 * 0.15])
X₆  = randomvariable("Normal",      "M",    [400,    400 * 0.05])
X₇  = randomvariable("LogNormal",   "M",    [80,      80 * 0.20])
X  = [X₁, X₂, X₃, X₄, X₅, X₆, X₇]

# Define correlation matrix:
ρˣ  = Matrix(1.0 * I, 7, 7)

# Define limit state function:
g(x::Vector) = x[1] * x[2] * x[3] - (x[1] ^ 2 * x[2] ^ 2 * x[4]) / (x[5] * x[6]) - x[7]

# Define reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Solve reliability problem using Subset Simulation Method:
Solution = solve(Problem, SSM())
Solution.PoF₂