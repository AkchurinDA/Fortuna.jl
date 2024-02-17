function analyze(Problem::SensitivityProblem)
    # Extract the problem data:
    X   = Problem.X
    ρˣ  = Problem.ρˣ
    g   = Problem.g
    θ   = Problem.θ

    # Define a reliability problem for the FORM analysis:
    gFORM(xVariable) = g(xVariable, θ)
    FORMProblem = ReliabilityProblem(X, ρˣ, gFORM)

    # Solve the reliability problem using the FORM:
    FORMSolution    = analyze(FORMProblem, FORM())
    x               = FORMSolution.x[:, end]
    u               = FORMSolution.u[:, end]
    β               = FORMSolution.β

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Define gradient functions of the limit state function in X- and U-spaces:
    ∇gVariable(x, θ) = gradient(Unknown -> g(x, Unknown), θ)
    ∇GVariable(u, θ) = gradient(Unknown -> G(g, θ, NatafObject, Unknown), u)
    
    # Compute the sensitivities w.r.t. the reliability index:
    ∇β = ∇gVariable(x, θ) / norm(∇GVariable(u, θ))

    # Compute the sensitivities w.r.t. the probability of failure:
    ∇PoF = -pdf(Normal(0, 1), β) * ∇β

    return SensitivityProblemCache(FORMSolution, ∇β, ∇PoF)
end

function G(g::Function, θ::Vector{Float64}, NatafObject::NatafTransformation, USample::AbstractVector)
    # Transform samples:
    XSample = transformsamples(NatafObject, USample, "U2X")

    # Evaluate the limit state function at the transform samples:
    GSample = g(XSample, θ)

    # Return the result:
    return GSample
end