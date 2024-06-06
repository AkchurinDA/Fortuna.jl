using Fortuna
using DelimitedFiles

# Define the random variables:
X₁ = randomvariable("Normal", "M", [29000, 0.05 * 29000]) # Young's modulus
X₂ = randomvariable("Normal", "M", [  110, 0.05 *   110]) # Moment of inertia about major axis
X  = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0; 0 1]

# Define the FE model of the cantilever beam:
OpenSeesPath  = "/Users/.../bin/OpenSees"                 # This must be an absolute path!
WorkDirectory = "/Users/.../Fortuna.jl/examples/OpenSees" # This must be an absolute path!
IFilename     = "CantileverBeam.tcl"
OFilename     = "Output.out"
Placeholders  = [":E", ":I"]
function CantileverBeam(x::Vector)
    # Inject values into the input file:
    IFileString = read(joinpath(WorkDirectory, IFilename), String)
    for (Placeholder, Value) in zip(Placeholders, x)
        IFileString = replace(IFileString, Placeholder => string(Value))
    end

    # Write the modified input file:
    TempIFilename = replace(IFilename, ".tcl" => "Temp.tcl")
    write(joinpath(WorkDirectory, TempIFilename), IFileString)

    # Run the model from the work directory:
    cd(WorkDirectory)
    run(pipeline(`$(OpenSeesPath) $(joinpath(WorkDirectory, TempIFilename))`, 
        stdout = devnull,
        stderr = devnull))
    
    # Extract the output:
    Δ = -readdlm(joinpath(WorkDirectory, OFilename))[end]

    # Delete the created files to prevent cluttering the work directory:
    rm(joinpath(WorkDirectory, TempIFilename))
    rm(joinpath(WorkDirectory,     OFilename))

    # Return the result:
    return Δ
end

# Define the limit state function:
g(x::Vector) = 1 - CantileverBeam(x)

# Define the reliability problem:
Problem = ReliabilityProblem(X, ρˣ, g)

# Perform the reliability analysis using the FORM:
FORMSolution = solve(Problem, FORM(), diff = :numeric)
println("FORM:")
println("β: $(FORMSolution.β)")
println("PoF: $(FORMSolution.PoF)")

# Perform the reliability analysis using the SORM:
SORMSolution = solve(Problem, SORM(), FORMSolution = FORMSolution, diff = :numeric)
println("SORM:")
println("β: $(SORMSolution.β₂[1]) (Hohenbichler and Rackwitz)")
println("β: $(SORMSolution.β₂[2]) (Breitung)")
println("PoF: $(SORMSolution.PoF₂[1]) (Hohenbichler and Rackwitz)")
println("PoF: $(SORMSolution.PoF₂[2]) (Breitung)")