# Preamble:
using Fortuna
using PyCall

# Load OpenSeesPy package:
ops = pyimport("openseespy.opensees")

# Define the random variables:
X₁ = randomvariable("Normal", "M", [29000, 0.05 * 29000]) # Young's modulus
X₂ = randomvariable("Normal", "M", [  110, 0.05 *   110]) # Moment of inertia about major axis
X  = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0; 0 1]

# Define the FE model of the cantilever beam:
function CantileverBeam(x::Vector)
    # Remove any previous models:
    ops.wipe()

    # Define the model parameters:
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    # Define the nodes:
    ops.node( 1,  0 * 18, 0)
    ops.node( 2,  1 * 18, 0)
    ops.node( 3,  2 * 18, 0)
    ops.node( 4,  3 * 18, 0)
    ops.node( 5,  4 * 18, 0)
    ops.node( 6,  5 * 18, 0)
    ops.node( 7,  6 * 18, 0)
    ops.node( 8,  7 * 18, 0)
    ops.node( 9,  8 * 18, 0)
    ops.node(10,  9 * 18, 0)
    ops.node(11, 10 * 18, 0)

    # Define the boundary conditions:
    ops.fix(1, 1, 1, 1)

    # Define the material properties:
    ops.uniaxialMaterial("Elastic", 1, 29000)

    # Define the cross-sectional properties:
    A = 9.12
    ops.section("Elastic",  1, x[1], A, x[2])

    # Define the transformation:
    ops.geomTransf("PDelta", 1)

    # Define the elements:
    ops.element("elasticBeamColumn",  1,  1,  2,  1, 1)
    ops.element("elasticBeamColumn",  2,  2,  3,  1, 1)
    ops.element("elasticBeamColumn",  3,  3,  4,  1, 1)
    ops.element("elasticBeamColumn",  4,  4,  5,  1, 1)
    ops.element("elasticBeamColumn",  5,  5,  6,  1, 1)
    ops.element("elasticBeamColumn",  6,  6,  7,  1, 1)
    ops.element("elasticBeamColumn",  7,  7,  8,  1, 1)
    ops.element("elasticBeamColumn",  8,  8,  9,  1, 1)
    ops.element("elasticBeamColumn",  9,  9, 10,  1, 1)
    ops.element("elasticBeamColumn", 10, 10, 11,  1, 1)

    # Define the loads:
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(11,   0, -1, 0)
    ops.load(11, -50,  0, 0)

    # Define the solver parameters:
    ops.system("BandSPD")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.algorithm("Linear")

    # Solve:
    ops.integrator("LoadControl", 0.01)
    ops.analysis("Static")
    ops.analyze(100)

    # Get the vertical displacement at the free end:
    Δ = -ops.nodeDisp(11, 2)

    # Return the result:
    return Δ 
end

# Define the limit state function:
g(x::Vector) = 1 - CantileverBeam(x)

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using the FORM:
Solution = solve(Problem, FORM(), Differentiation = :Numeric)
println("FORM:")
println("β: $(Solution.β)")
println("PoF: $(Solution.PoF)")

# Perform the reliability analysis using the SORM:
Solution = solve(Problem, SORM(), Differentiation = :Numeric)
println("SORM:")
println("β: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β: $(Solution.β₂[2]) (Breitung)")
println("PoF: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF: $(Solution.PoF₂[2]) (Breitung)")

# Preamble:
using Surrogates

# Define the training points:
LowerBound = [
    29000 - 6 * 0.05 * 29000,
      110 - 6 * 0.05 *   110]
UpperBound = [
    29000 + 6 * 0.05 * 29000,
      110 + 6 * 0.05 *   110]
XTrain = sample(50, LowerBound, UpperBound, SobolSample())
YTrain = [CantileverBeam([x...]) for x in XTrain]

# Fit a Kriging surrogate model to the training points:
CantileverBeamSurrogate = Kriging(XTrain, YTrain, LowerBound, UpperBound)

# Define the limit state function using the surrogate model:
gSurrogate(x::Vector) = 1 - CantileverBeamSurrogate(x)

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, gSurrogate)

# Perform the reliability analysis using the FORM:
Solution = solve(Problem, FORM())
println("FORM:")
println("β: $(Solution.β)")
println("PoF: $(Solution.PoF)")

# Perform the reliability analysis using the SORM:
Solution = solve(Problem, SORM())
println("SORM:")
println("β: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β: $(Solution.β₂[2]) (Breitung)")
println("PoF: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF: $(Solution.PoF₂[2]) (Breitung)")