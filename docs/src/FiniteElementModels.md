# Finite Element Models

The real power of `Fortuna.jl` package comes from the ability to define limit state functions using finite element (FE) models. For example, if the roof drift of a frame subjected to wind loads ``\Delta`` must not exceed a certain limit ``\Delta_{0}``, then the limit state function is given by ``g(\vec{X}) = \Delta_{0} - \Delta(\vec{X})``, where ``\Delta`` must be computed from the FE model of this frame with a proper consideration of all random variables involved in the problem formulation ``\vec{X}``. 

`Fortuna.jl` package allows you to solve reliability problems for FE models built using any Julia package with FE modeling capabilities, e.g., [`InstantFrame.jl`](https://github.com/runtosolve/InstantFrame.jl) and [`ONSAS.jl`](https://github.com/ONSAS/ONSAS.jl) packages for structural system problems, [`PowerSystems.jl`](https://github.com/NREL-Sienna/PowerSystems.jl) package for power system problems, [`Gridap.jl`] package for more general FE problems, and many more. Moreover, thanks to the great infrastructure Julia-to-Python developed by the [`JuliaPy`](https://github.com/JuliaPy) organization, you can also build FE models using Python packages with such capabilities. The most prominent example of such Python package is, of course, [`OpenSeesPy`](https://github.com/zhuminjie/OpenSeesPy) which serves as the gold standard to simulate the performance of structural and geotechnical systems subjected to earthquakes.

## `OpenSeesPy`

### Installation

To be able to work with [`OpenSeesPy`](https://openseespydoc.readthedocs.io/en/latest/) package directly from Julia you need to install two packages that allow Julia to talk to Python: 

| Item | Description |
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
Conda.pip_interop(true)
Conda.pip("install", "openseespy")
```

Now you can work with `OpenSeesPy` package directly from Julia!

!!! note
    If you are experiencing any problem please refer to Conda.jl package's [documentation](https://github.com/JuliaPy/Conda.jl).

### Example

Consider the following example of a cantilever beam under axial and transverse loads based on the example provided in [Denavit:2013](@citet) with the only difference that Young's modulus ``E = X_{1}`` and moment of inertia about major axis ``I = X_{2}`` are normally-distributed random variables with the following mean and coefficient of variation values:

| Random Variable | Mean, ``\mu`` | Coefficient of variation, ``V`` |
| :--- | :--- | :--- |
| ``E = X_{1}`` | ``29000 \text{ ksi}`` | ``0.05`` |
| ``I = X_{2}`` | ``110 \text{ in.}^{4}`` | ``0.05`` |

Also, let's define the limit state function as 

```math
g(\vec{X}) = \Delta_{0} - \Delta(\vec{X})
```

where ``\Delta(\vec{X})`` is the downward deflection at the free end, which must not exceed the limit of ``\Delta_{0} = 1 \text{ in.}``. The goal is to find the reliability indices ``\beta`` and probabilities of failure ``P_{f}`` for this problem using first- and second-order reliability method.

```@setup 1
using Conda
Conda.pip_interop(true)
Conda.pip("install", "openseespy")
```

```@example 1
# Preamble:
using Fortuna
using PyCall

# Load OpenSeesPy package:
ops = pyimport("openseespy.opensees")

# Define the random variables:
X₁ = randomvariable("Normal", "M", [29000, 0.05 * 29000]) # Young's modulus
X₂ = randomvariable("Normal", "M", [  110, 0.05 *   110]) # Moment of inertia about major-axis
X  = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0; 0 1]

# Define the limit state function:
function g(x::Vector)
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
    ops.integrator("LoadControl", 0.001)
    ops.analysis("Static")
    ops.analyze(1000)

    # Get the vertical displacement at the free end:
    Δ = -ops.nodeDisp(11, 2)

    return 1 - Δ
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using the FORM:
Solution = solve(Problem, FORM(), Differentiation = :Numeric)
println("FORM:")
println("β:   $(Solution.β)  ")
println("PoF: $(Solution.PoF)")

# Perform the reliability analysis using the SORM:
Solution = solve(Problem, SORM(), Differentiation = :Numeric)
println("SORM:")
println("β:   $(Solution.β₂)  ")
println("PoF: $(Solution.PoF₂)")
```

## `InstantFrame.jl`

### Installation

To install `InstantFrame.jl` package, type `]` in Julia REPL to enter the built-in Julia package manager and execute the following command:

```
pkg> add https://github.com/runtosolve/InstantFrame.jl.git
```

### Example

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

# Define the limit state function:
function g(x::Vector)
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
            ( 0 * 18, 0, 0), 
            ( 1 * 18, 0, 0), 
            ( 2 * 18, 0, 0), 
            ( 3 * 18, 0, 0), 
            ( 4 * 18, 0, 0), 
            ( 5 * 18, 0, 0), 
            ( 6 * 18, 0, 0), 
            ( 7 * 18, 0, 0), 
            ( 8 * 18, 0, 0), 
            ( 9 * 18, 0, 0), 
            (10 * 18, 0, 0)])

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

    return 1 - Δ
end

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using the FORM:
Solution = solve(Problem, FORM(), Differentiation = :Numeric)
println("FORM:")
println("β:   $(Solution.β)  ")
println("PoF: $(Solution.PoF)")

# Perform the reliability analysis using the SORM:
Solution = solve(Problem, SORM(), Differentiation = :Numeric)
println("SORM:")
println("β:   $(Solution.β₂)  ")
println("PoF: $(Solution.PoF₂)")
```