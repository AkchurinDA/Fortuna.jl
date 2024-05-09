"""
    randomvariable(Distribution::AbstractString, DefineBy::AbstractString, Values::Union{Real, AbstractVector{<:Real}})

Function used to define random variables.
"""
function randomvariable(Distribution::AbstractString, DefineBy::AbstractString, Values::Union{Real, AbstractVector{<:Real}})
    # Supported distributions:
    SupportedDistributions = ["Exponential", "Frechet", "Gamma", "Gumbel", "LogNormal", "Normal", "Poisson", "Uniform", "Weibull"]

    # Error-catching:
    Distribution in SupportedDistributions || throw(ArgumentError("Provided distribution is not supported!"))
    (DefineBy == "M" || DefineBy == "P")   || throw(ArgumentError("""Random variables can only be defined by "Moments" ("M") and "Parameters ("P")"!"""))

    # Convert moments into parameters if needed:
    DefineBy == "M" && (Values = convertmoments(Distribution, Values))

    # Define a random variable:
    RandomVariable = getfield((@__MODULE__).Distributions, Symbol(Distribution))(Values...)

    # Return the result:
    return RandomVariable
end

function convertmoments(Distribution::AbstractString, Moments::Union{Real, AbstractVector{<:Real}})
    # Error-catching:
    length(Moments) == 2 || throw(ArgumentError("Too few or many moments are provided! Provide only the mean (μ) and standard deviation (σ) in a vector format (i.e., [μ, σ])!"))

    # Convert moments to parameters:
    if Distribution == "Exponential"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Error catching:
        if Mean != STD
            throw(DomainError(Moments, "Mean and standard deviation values of must be the same!"))
        end

        # Convert moments to parameters:
        θ          = Mean
        Parameters = θ
    end

    if Distribution == "Frechet"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        FFrechet(u, p) = sqrt(SpecialFunctions.gamma(1 - 2 / u) - SpecialFunctions.gamma(1 - 1 / u) ^ 2) / SpecialFunctions.gamma(1 - 1 / u) - p
        u₀             = (2 + 1E-1, 1E+6)
        p₀             = STD / Mean
        Problem        = NonlinearSolve.IntervalNonlinearProblem(FFrechet, u₀, p₀)
        Solution       = NonlinearSolve.solve(Problem, nothing, abstol = 1E-9, reltol = 1E-9)
        α              = Solution.u
        if !isapprox(FFrechet(α, p₀), 0, atol = 1E-9)
            throw(DomainError(Moments, "Conversion of the provided moments to parameters has failed!"))
        end
        θ              = Mean / SpecialFunctions.gamma(1 - 1 / α)
        Parameters     = [α, θ]
    end

    if Distribution == "Gamma"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        α          = Mean ^ 2 / STD ^ 2
        θ          = STD ^ 2 / Mean
        Parameters = [α, θ]
    end
    
    if Distribution == "Gumbel"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        γ          = Base.MathConstants.eulergamma
        μ          = Mean - (STD * γ * sqrt(6)) / π
        θ          = (STD * sqrt(6)) / π
        Parameters = [μ, θ]
    end
    
    if Distribution == "LogNormal"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        μ          = log(Mean) - log(sqrt(1 + (STD / Mean) ^ 2))
        σ          = sqrt(log(1 + (STD / Mean) ^ 2))
        Parameters = [μ, σ]
    end
    
    if Distribution == "Normal"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        μ          = Mean
        σ          = STD
        Parameters = [μ, σ]
    end
    
    if Distribution == "Poisson"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Error catching:
        if !(Mean ≈ (STD ^ 2))
            throw(DomainError(Moments, "Standard deviation must be equal to square root of mean!"))
        end

        # Convert moments to parameters:
        λ          = Mean
        Parameters = λ
    end
    
    if Distribution == "Uniform"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        a          = Mean - STD * sqrt(3)
        b          = Mean + STD * sqrt(3)
        Parameters = [a, b]
    end
    
    if Distribution == "Weibull"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        FWeibull(u, p) = sqrt(SpecialFunctions.gamma(1 + 2 / u) - SpecialFunctions.gamma(1 + 1 / u) ^ 2) / SpecialFunctions.gamma(1 + 1 / u) - p
        u₀             = (1E-1, 1E+6)
        p₀             = STD / Mean
        Problem        = NonlinearSolve.IntervalNonlinearProblem(FWeibull, u₀, p₀)
        Solution       = NonlinearSolve.solve(Problem, nothing, abstol = 1E-9, reltol = 1E-9)
        α              = Solution.u
        if !isapprox(FWeibull(α, p₀), 0, atol = 1E-9)
            throw(DomainError(Moments, "Conversion of the provided moments to parameters has failed!"))
        end
        θ          = Mean / SpecialFunctions.gamma(1 + 1 / α)
        Parameters = [α, θ]
    end

    # Return the result:
    return Parameters
end