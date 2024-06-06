# Advanced Usage

## Finite Element Models

The real power of `Fortuna.jl` package comes from the ability to define limit state functions using finite element (FE) models of complex systems. For example, if the roof drift of a frame subjected to wind loads ``\Delta`` must not exceed a certain limit ``\Delta_{0}``, then the limit state function is given by ``g(\vec{X}) = \Delta_{0} - \Delta(\vec{X})``, where ``\Delta`` must be computed from the FE model of this frame with a proper consideration of all random variables involved in the problem formulation ``\vec{X}``. 

`Fortuna.jl` package allows you to solve reliability problems for FE models built using any Julia package with such capabilities, e.g., [`InstantFrame.jl`](https://github.com/runtosolve/InstantFrame.jl) and [`ONSAS.jl`](https://github.com/ONSAS/ONSAS.jl) packages for structural system problems, [`PowerSystems.jl`](https://github.com/NREL-Sienna/PowerSystems.jl) package for power system problems, [`Gridap.jl`](https://github.com/gridap/Gridap.jl) package for more general FE problems, and many more. Moreover, thanks to the great Julia-to-Python infrastructure developed by the [`JuliaPy`](https://github.com/JuliaPy) organization, you can also build FE models using Python packages with such capabilities. The most prominent example of such Python package is, of course, [`OpenSeesPy`](https://github.com/zhuminjie/OpenSeesPy) which serves as the standard to simulate the performance of structural and geotechnical systems subjected to earthquakes.

### `OpenSeesPy`

#### Installation

To be able to work with [`OpenSeesPy`](https://openseespydoc.readthedocs.io/en/latest/) package directly from Julia you need to install two packages that allow Julia to talk to Python: 

| Package | Description |
| :--- | :--- |
| [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl) | Let's you call Python functions directly from Julia |
| [`Conda.jl`](https://github.com/JuliaPy/Conda.jl) | Let's you install Python packages from `conda` and `pip` package managers that can be then called directly from Julia using `PyCall.jl` package |

To install `PyCall.jl` and `Conda.jl` packages, type `]` in Julia REPL to enter the built-in Julia package manager and execute the following command:

```
pkg> add PyCall Conda
```

Now you need to install `OpenSeesPy` package through `Conda.jl` package. To do that, run the following commands in Julia REPL:

```julia
using Conda
Conda.pip_interop(true)            # Activate "pip"
Conda.pip("install", "openseespy") # Install OpenSeesPy using "pip"
```

Now you can work with `OpenSeesPy` package directly from Julia!

!!! note
    If you are experiencing any problems installing OpenSeesPy package please refer to Conda.jl package's [documentation](https://github.com/JuliaPy/Conda.jl).

#### Example

Consider a cantilever beam subjected to simultaneous axial and transverse loading based on the example provided in [Denavit:2013](@citet) with the only difference that Young's modulus ``E = X_{1}`` and moment of inertia about major axis ``I = X_{2}`` are uncorrelated normally-distributed random variables. The cross-sectional area of the beam ``A`` is ``9.12 \text{ in.}^{2}``.

```@raw html
<img src="../assets/Examples-OpenSees-1.png" class="center" style="max-height:350px; border-radius:2.5px;"/>
```

Let's define the limit state function as 

```math
g(\vec{X}) = \Delta_{0} - \Delta(\vec{X})
```

where ``\Delta(\vec{X})`` is the downward deflection at the free end of the beam, which must not exceed the deflection limit ``\Delta_{0}`` of ``1 \text{ in.}``. The goal of this example is to find the reliability indices ``\beta`` and probabilities of failure ``P_{f}`` using First- and Second-Order Reliability Methods.

```@setup 1
# Install OpenSeesPy:
using Conda
Conda.pip_interop(true)
Conda.pip("install", "openseespy")

# Force Julia to use its own Python distribution via Conda.jl:
using Pkg
ENV["PYTHON"] = ""
Pkg.build("PyCall")
```

```@example 1
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
    ops.node( 1,   0, 0)
    ops.node( 2,  18, 0)
    ops.node( 3,  36, 0)
    ops.node( 4,  54, 0)
    ops.node( 5,  72, 0)
    ops.node( 6,  90, 0)
    ops.node( 7, 108, 0)
    ops.node( 8, 126, 0)
    ops.node( 9, 144, 0)
    ops.node(10, 162, 0)
    ops.node(11, 180, 0)

    # Define the boundary conditions:
    ops.fix(1, 1, 1, 1)

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
Solution = solve(Problem, FORM(), diff = :numeric)
println("FORM:")
println("β: $(Solution.β)")
println("PoF: $(Solution.PoF)")

# Perform the reliability analysis using the SORM:
Solution = solve(Problem, SORM(), diff = :numeric)
println("SORM:")
println("β: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β: $(Solution.β₂[2]) (Breitung)")
println("PoF: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF: $(Solution.PoF₂[2]) (Breitung)")
```

Observe that `Differentiation = :Numeric` keyword argument was used in `solve()` function. This forces `Fortuna.jl` package to use numerical differentiation to evaluate the gradients and Hessians of the limit state function since the finite element (FE) model of the cantilever beam, defined using `OpenSeesPy` package, is not automatically differentiable. Note that you would still obtain a solution if you didn't use `Differentiation = :Numeric` keyword argument, but each time `Fortuna.jl` package needs to differentiate the limit state function, it would (1) attempt to differentiate it automatically, and after it fails, (2) it would differentiate it numerically, significantly slowing down the reliability analysis.

### `InstantFrame.jl`

#### Installation

To install `InstantFrame.jl` package, type `]` in Julia REPL to enter the built-in Julia package manager and execute the following command:

```
pkg> add InstantFrame
```

#### Example

Let's solve previously defined problem using `InstantFrame.jl` package.

```@example 1
# Preamble:
using Fortuna
using InstantFrame

# Define the random variables:
X₁ = randomvariable("Normal", "M", [29000, 0.05 * 29000]) # Young's modulus
X₂ = randomvariable("Normal", "M", [  110, 0.05 *   110]) # Moment of inertia about major-axis
X  = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0; 0 1]

# Define the FE model of the cantilever beam:
function CantileverBeam(x::Vector)
    # Define the material properties:
    Material = InstantFrame.Material(
        names = ["Steel"], 
        E     = [x[1]], 
        ν     = [0.3], 
        ρ     = [0])
    
    # Define the cross-sectional properties:
    CrossSection = InstantFrame.CrossSection(
        names = ["Beam"], 
        A     = [9.12], 
        Iy    = [0], 
        Iz    = [x[2]], 
        J     = [0])

    # Define the element-to-element connections:
    Connection = InstantFrame.Connection(
        names     = ["Rigid"], 
        stiffness = (ux = [Inf], uy = [Inf], uz = [Inf], rx = [Inf], ry = [Inf], rz = [Inf]))

    # Define the nodes:
    Node = InstantFrame.Node(
        numbers     = 1:11, 
        coordinates = [
            (  0, 0, 0), 
            ( 18, 0, 0), 
            ( 36, 0, 0), 
            ( 54, 0, 0), 
            ( 72, 0, 0), 
            ( 90, 0, 0), 
            (108, 0, 0), 
            (126, 0, 0), 
            (144, 0, 0), 
            (162, 0, 0), 
            (180, 0, 0)])

    # Define the elements:
    Element = InstantFrame.Element(
        types         = ["", "", "", "", "", "", "", "", "", ""],
        numbers       = 1:10, 
        nodes         = [
            ( 1,  2), 
            ( 2,  3), 
            ( 3,  4), 
            ( 4,  5), 
            ( 5,  6), 
            ( 6,  7), 
            ( 7,  8), 
            ( 8,  9), 
            ( 9, 10), 
            (10, 11)], 
        orientation   = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
        connections   = [
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid"), 
            ("Rigid", "Rigid")], 
        cross_section = [ "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam",  "Beam"], 
        material      = ["Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel", "Steel"])

    # Define the boundary conditions:
    Support = InstantFrame.Support(
        nodes     = 1:11, 
        stiffness = (
            uX = [Inf,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], 
            uY = [Inf,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], 
            uZ = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf],
            rX = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf], 
            rY = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf], 
            rZ = [Inf,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0]))

    # Define the distributed loads:
    UniformLoad = InstantFrame.UniformLoad(nothing)

    # Define the concentrated loads:
    PointLoad = InstantFrame.PointLoad(
        labels     = [""], 
        nodes      = [11], 
        magnitudes = (FX = [-50], FY = [-1], FZ = [0], MX = [0], MY = [0], MZ = [0]))

    # Solve:
    Model = InstantFrame.solve(Node, CrossSection, Material, Connection, Element, Support, UniformLoad, PointLoad, analysis_type = "second order", solution_tolerance = 1E-6)

    # Get the vertical displacement at the free end:
    Δ = -Model.solution.displacements[11][2]

    return Δ
end

# Define the limit state function:
g(x::Vector) = 1 - CantileverBeam(x)

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using the FORM:
Solution = solve(Problem, FORM(), diff = :numeric)
println("FORM:")
println("β: $(Solution.β)")
println("PoF: $(Solution.PoF)")

# Perform the reliability analysis using the SORM:
Solution = solve(Problem, SORM(), diff = :numeric)
println("SORM:")
println("β: $(Solution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β: $(Solution.β₂[2]) (Breitung)")
println("PoF: $(Solution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF: $(Solution.PoF₂[2]) (Breitung)")
```

## Surrogate Models

Most of the time, FE models are very expensive to evaluate, especially when it comes to differentiation of such FE models for the purposes of reliability analysis. Instead, it is possible to build a surrogate of a FE model. For example, [`Surrogates.jl`](https://github.com/SciML/Surrogates.jl) package allows you to build surrogate models. The most important feature of `Surrogates.jl` package is that it permits automatic differentiation of a surrogates model even if the underlying FE model cannot be differentiated automatically, which significantly speeds up the reliability analysis!

#### Example

Let's build a Kriging surrogate model of the FE model of the cantilever beam studied earlier and solve for the reliability indices ``\beta`` and probabilities of failure ``P_{f}`` using First- and Second-Order Reliability Methods again.

```@example 1
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
ĝ(x::Vector) = 1 - CantileverBeamSurrogate(x)

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, ĝ)

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
```