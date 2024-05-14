# Preamble:
using Fortuna
using PyCall

# Load Python modules:
ops = pyimport("openseespy.opensees")

# Define the random vector:
t₁ = randomvariable("Normal", "M", [   20,    20 * 0.01])
t₂ = randomvariable("Normal", "M", [   50,    50 * 0.01])
P₁ = randomvariable("Gumbel", "M", [20000, 20000 * 0.01])
P₂ = randomvariable("Gumbel", "M", [20000, 20000 * 0.01])
P₃ = randomvariable("Gumbel", "M", [20000, 20000 * 0.01])
X  = [t₁, t₂, P₁, P₂, P₃]

# Define the correlation matrix:
ρˣ = [
    1.0 0.5 0.0 0.0 0.0
    0.5 1.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0 0.0
    0.0 0.0 0.0 1.0 0.0
    0.0 0.0 0.0 0.0 1.0]

# Define the limit state function:
function g(x::Vector)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    ops.node(1,     0,    0)
    ops.node(2,  4000, 2000)
    ops.node(3,  8000, 4000)
    ops.node(4,  8000,    0)
    ops.node(5, 12000, 2000)
    ops.node(6, 16000,    0)

    ops.fix(1, 1, 1)
    ops.fix(6, 0, 1)

    ops.uniaxialMaterial("Elastic", 1, 210000)

    A = x[1] * x[2]
    ops.element("Truss", 1, 1, 2, A, 1)
    ops.element("Truss", 2, 1, 4, A, 1)
    ops.element("Truss", 3, 2, 4, A, 1)
    ops.element("Truss", 4, 2, 3, A, 1)
    ops.element("Truss", 5, 3, 4, A, 1)
    ops.element("Truss", 6, 3, 5, A, 1)
    ops.element("Truss", 7, 4, 5, A, 1)
    ops.element("Truss", 8, 4, 6, A, 1)
    ops.element("Truss", 9, 5, 6, A, 1)

    ops.timeSeries("Linear", 1)
    
    ops.pattern("Plain", 1, 1)
    
    ops.load(2, 0, -x[3])
    ops.load(3, 0, -x[4])
    ops.load(5, 0, -x[5])

    ops.system("BandSPD")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1)
    ops.algorithm("Linear")
    ops.analysis("Static")

    ops.analyze(1)

    Δ = -ops.nodeDisp(2, 2)

    return 10 - Δ
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Solve the reliability problem using the FORM:
Solution = solve(Problem, FORM(), Differentiation = :Numeric)
println("β   = $(Solution.β)  ")
println("PoF = $(Solution.PoF)")
