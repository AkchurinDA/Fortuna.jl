"""
    randomvariable(Distribution::AbstractString, DefineBy::AbstractString, Values::Union{Real, AbstractVector{<:Real}})

Function used to define random variables.
"""
function randomvariable(Distribution::AbstractString, DefineBy::AbstractString, Values::Union{Real, AbstractVector{<:Real}})
    # Convert strings to lowercase:
    Distribution = lowercase(Distribution)
    DefineBy     = lowercase(DefineBy)

    # Error-catching:
    if DefineBy != "m" && DefineBy != "p"
        throw(ArgumentError("""Random variables can only be defined by "Moments" ("M") and "Parameters ("P")"!"""))
    end

    # Convert moments of a random variable into its parameters if the random variable is defined by its moments:
    if DefineBy == "m"
        Values = convert(Distribution, Values)
    end

    # Create a random variable:
    if Distribution == "exponential"
        RandomVariable = Distributions.Exponential(Values...)
    elseif Distribution == "frechet"
        RandomVariable = Distributions.Frechet(Values...)
    elseif Distribution == "gamma"
        RandomVariable = Distributions.Gamma(Values...)
    elseif Distribution == "gumbel"
        RandomVariable = Distributions.Gumbel(Values...)
    elseif Distribution == "lognormal"
        RandomVariable = Distributions.LogNormal(Values...)
    elseif Distribution == "normal"
        RandomVariable = Distributions.Normal(Values...)
    elseif Distribution == "poisson"
        RandomVariable = Distributions.Poisson(Values...)
    elseif Distribution == "uniform"
        RandomVariable = Distributions.Uniform(Values...)
    elseif Distribution == "weibull"
        RandomVariable = Distributions.Weibull(Values...)
    else
        error("Provided distribution is not supported!")
    end

    # Return the result:
    return RandomVariable
end

function convert(Distribution::AbstractString, Moments::Union{Real, AbstractVector{<:Real}})
    # Convert strings to lowercase:
    Distribution = lowercase(Distribution)

    # Error-catching:
    if length(Moments) > 2
        throw(ArgumentError("Too many moments are provided!"))
    end

    # Convert moments to parameters:
    if Distribution == "exponential"
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
    elseif Distribution == "frechet"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        FFrechet(u, p) = sqrt(SpecialFunctions.gamma(1 - 2 / u) - SpecialFunctions.gamma(1 - 1 / u) ^ 2) / SpecialFunctions.gamma(1 - 1 / u) - STD / Mean
        u₀             = (2 + 0.1, 10000)
        Problem        = NonlinearSolve.IntervalNonlinearProblem(FFrechet, u₀)
        Solution       = NonlinearSolve.solve(Problem, NonlinearSolve.Bisection(), abstol = 10 ^ (-9))
        if !isapprox(FFrechet(Solution.u, 0), 0, atol = 10 ^ (-9))
            throw(DomainError(Moments, "Conversion of the provided moments to parameters has failed!"))
        end
        α              = Solution.u
        θ              = Mean / SpecialFunctions.gamma(1 - 1 / α)
        Parameters     = [α, θ]
    elseif Distribution == "gamma"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        α          = Mean ^ 2 / STD ^ 2
        θ          = STD ^ 2 / Mean
        Parameters = [α, θ]
    elseif Distribution == "gumbel"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        γ          = Base.MathConstants.eulergamma
        μ          = Mean - (STD * γ * sqrt(6)) / π
        θ          = (STD * sqrt(6)) / π
        Parameters = [μ, θ]
    elseif Distribution == "lognormal"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        μ          = log(Mean) - log(sqrt(1 + (STD / Mean) ^ 2))
        σ          = sqrt(log(1 + (STD / Mean) ^ 2))
        Parameters = [μ, σ]
    elseif Distribution == "normal"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        μ          = Mean
        σ          = STD
        Parameters = [μ, σ]
    elseif Distribution == "poisson"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Error catching:
        if Mean != STD
            throw(DomainError(Moments, "Mean and standard deviation values of must be the same!"))
        end

        # Convert moments to parameters:
        λ          = Mean
        Parameters = λ
    elseif Distribution == "uniform"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        a          = Mean - STD * sqrt(3)
        b          = Mean + STD * sqrt(3)
        Parameters = [a, b]
    elseif Distribution == "weibull"
        # Extract moments:
        Mean = Moments[1]
        STD  = Moments[2]

        # Convert moments to parameters:
        FWeibull(u, p) = sqrt(SpecialFunctions.gamma(1 + 2 / u) - SpecialFunctions.gamma(1 + 1 / u) ^ 2) / SpecialFunctions.gamma(1 + 1 / u) - STD / Mean
        u₀             = (0.1, 10000)
        Problem        = NonlinearSolve.IntervalNonlinearProblem(FWeibull, u₀)
        Solution       = NonlinearSolve.solve(Problem, NonlinearSolve.Bisection(), abstol = 10 ^ (-9))
        if !isapprox(FWeibull(Solution.u, 0), 0, atol = 10 ^ (-9))
            throw(DomainError(Moments, "Conversion of the provided moments to parameters has failed!"))
        end
        α              = Solution.u
        θ              = Mean / SpecialFunctions.gamma(1 + 1 / α)
        Parameters     = [α, θ]
    else
        error("Provided distribution is not supported!")
    end

    # Return the result:
    return Parameters
end