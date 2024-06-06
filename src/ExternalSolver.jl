struct ExternalModel
    Solver        ::Symbol
    SolverPath    ::String
    WorkDirectory ::String
    InputFilename ::String
    OutputFilename::String
    Placeholders  ::Vector{String}
    InjectedValues::Vector{Float64}
    Monitor       ::Function
end
