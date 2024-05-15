# http://www1.coe.neu.edu/~jfhajjar/home/Denavit%20and%20Hajjar%20-%20Geometric%20Nonlinearity%20in%20OpenSees%20-%20Report%20No.%20NEU-CEE-2013-02%202013.pdf
# Figure 5 (Major-axis bending)

# Preamble:
using Fortuna
using InstantFrame

# Define the random vector:
X1 = randomvariable("Normal", "M", [29000, 0.05 * 29000])
X2 = randomvariable("Normal", "M", [  110, 0.05 *   110])
X  = [X1, X2]

# Define the correlation matrix:
ρˣ = [1 0; 0 1]

# Define the limit state function:
function g(x::Vector)
    Material = InstantFrame.Material(
        names = ["Steel"], 
        E     = [x[1]], 
        ν     = [0.3], 
        ρ     = [492 / 32.170 / 12 ^ 4 / 1000])
        
    CrossSection = InstantFrame.CrossSection(
        names = ["Beam"], 
        A     = [9.12], 
        Iy    = [37.1], 
        Iz    = [x[2]], 
        J     = [0.001])

    Connection = InstantFrame.Connection(
        names     = ["Rigid"], 
        stiffness = (ux = [Inf], uy = [Inf], uz = [Inf], rx = [Inf], ry = [Inf], rz = [Inf]))

    Node = InstantFrame.Node(
        numbers     = [1, 2], 
        coordinates = [
            (  0, 0, 0), 
            (180, 0, 0)])

    Element = InstantFrame.Element(
        types         = [""],
        numbers       = [1], 
        nodes         = [(1, 2)], 
        orientation   = [0], 
        connections   = [("Rigid", "Rigid")], 
        cross_section = ["Beam"], 
        material      = ["Steel"])

    Support = InstantFrame.Support(
        nodes     = [1, 2], 
        stiffness = (uX = [Inf, 0], uY = [Inf, 0], uZ = [Inf, Inf], rX = [Inf, Inf], rY=[Inf, Inf], rZ=[Inf, 0]))

    UniformLoad = InstantFrame.UniformLoad(
        labels     = ["Test"], 
        elements   = [1], 
        magnitudes = (qX = [0, 0], qY = [0, 0], qZ = [0, 0], mX = [0, 0], mY = [0, 0], mZ = [0, 0]))

    PointLoad = InstantFrame.PointLoad(
        labels     = ["Test"], 
        nodes      = [2], 
        magnitudes = (FX = [-50], FY=[1], FZ=[0], MX=[0], MY=[0], MZ=[0]))

    Model = InstantFrame.solve(Node, CrossSection, Material, Connection, Element, Support, UniformLoad, PointLoad, analysis_type = "second order", solution_tolerance = 1E-6)

    Δ = Model.solution.displacements[2][2]

    return 1 - Δ
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Solve the reliability problem:
Solution = solve(Problem, FORM(), Differentiation = :Numeric)
Solution.β   # 2.784322708521028
Solution.PoF # 0.002681981827787

# Redefine the limit state function:
H = 1
P = 50
L = 180
function g(x::Vector)
    α = sqrt((P * L ^ 2) / (x[1] * x[2]))
    Δ = (H * L ^ 3) / (3 * x[1] * x[2]) * ((3 * (tan(α) - α)) / α ^ 3)
    return 1 - Δ
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Solve the reliability problem:
Solution = solve(Problem, FORM())
Solution.β   # 2.778811787854018
Solution.PoF # 0.002727906330022