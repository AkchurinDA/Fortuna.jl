using Fortuna
using DelimitedFiles

# Define the random variables:
X₁ = randomvariable("Normal", "M", [29000, 0.10 * 29000]) # Young's modulus
X₂ = randomvariable("Normal", "M", [    2, 0.10 *     2]) # Lateral load
X  = [X₁, X₂]

# Define the correlation matrix:
ρˣ = [1 0; 0 1]

# Define the FE model of the frame under fire:
WorkDirectory = "C:\\Users\\...\\Fortuna.jl\\examples\\SAFIR" # This must be an absolute path!
IFilename     = "FrameUnderFire.IN"
OFilename     = "FrameUnderFireTemp.OUT"
Placeholders  = [":E", ":F"]
function FrameUnderFire(x::Vector)
    # Inject values into the input file:
    IFileString = read(joinpath(WorkDirectory, IFilename), String)
    for (Placeholder, Value) in zip(Placeholders, x)
        IFileString = replace(IFileString, Placeholder => string(Value))
    end

    # Write the modified input file:
    TempIFilename = replace(IFilename, ".IN" => "Temp.IN")
    write(joinpath(WorkDirectory, TempIFilename), IFileString)

    # Run the model from the work directory:
    cd(WorkDirectory)
    run(pipeline(`cmd /C "SAFIR $(replace(TempIFilename, ".IN" => ""))"`, 
        stdout = devnull,
        stderr = devnull))
    
    # Extract the output:
    OFileString = read(OFilename, String)
    SIndex      = findlast("TOTAL DISPLACEMENTS.\r\n --------------------\r\n NODE    DOF 1     DOF 2     DOF 3     DOF 4     DOF 5     DOF 6     DOF 7\r\n", OFileString)
    OFileString = OFileString[(SIndex[end] + 1):end]
    FIndex      = findfirst("\r\n\r\n", OFileString)
    OFileString = OFileString[1:(FIndex[1] - 1)]
    write(joinpath(WorkDirectory, OFilename), OFileString)
    Δ = -readdlm(joinpath(WorkDirectory, OFilename))[38, 3]

    # Delete the created files to prevent cluttering the work directory:
    rm(joinpath(WorkDirectory, TempIFilename))
    rm(joinpath(WorkDirectory,     OFilename))

    # Return the result:
    return Δ
end

# Define the limit state function:
g(x::Vector) = 0.075 - FrameUnderFire(x)

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