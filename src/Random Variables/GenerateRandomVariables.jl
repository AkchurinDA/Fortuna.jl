"""
    randomvariable(DistributionName::AbstractString, DefineBy::AbstractString, Values::Union{Real, AbstractVector{<:Real}})

Function used to define random variables.
"""
function randomvariable(DistributionName::AbstractString, DefineBy::AbstractString, Values::Union{Real, AbstractVector{<:Real}})
    # Convert strings to lowercase:
    DistributionName = lowercase(DistributionName)
    DefineBy = lowercase(DefineBy)

    # Error-catching:
    if DefineBy != "m" && DefineBy != "p"
        error("""Random variables can only be defined by "Moments" ("M") and "Parameters ("P")".""")
    end

    # Convert moments of a random variable into its parameters if the random variable is defined by its moments:
    if DefineBy == "m"
        Values = convert(DistributionName, Values)
    end

    # Create a random variable:
    if DistributionName == "exponential"
        RandomVariable = Distributions.Exponential(Values)
    elseif DistributionName == "gamma"
        RandomVariable = Distributions.Gamma(Values[1], Values[2])
    elseif DistributionName == "gumbel"
        RandomVariable = Distributions.Gumbel(Values[1], Values[2])
    elseif DistributionName == "lognormal"
        RandomVariable = Distributions.LogNormal(Values[1], Values[2])
    elseif DistributionName == "normal"
        RandomVariable = Distributions.Normal(Values[1], Values[2])
    elseif DistributionName == "poisson"
        RandomVariable = Distributions.Poisson(Values[1])
    elseif DistributionName == "uniform"
        RandomVariable = Distributions.Uniform(Values[1], Values[2])
    elseif DistributionName == "weibull"
        RandomVariable = Distributions.Weibull(Values[1], Values[2])
    else
        error("Provided distribution is not supported.")
    end

    # Return the result:
    return RandomVariable
end

function convert(DistributionName::AbstractString, Moments::Union{Real, AbstractVector{<:Real}})
    # Convert strings to lowercase:
    DistributionName = lowercase(DistributionName)

    # Error-catching:
    if length(Moments) > 2
        error("Too many moments are provided.")
    end

    # Convert moments to parameters:
    if DistributionName == "exponential"
        # Extract moments:
        Mean        = Moments[1]
        STD         = Moments[2]

        # Error catching:
        if Mean != STD
            error("Mean and standard deviation values of an exponential random variables must be the same.")
        end

        # Convert moments to parameters:
        θ           = Mean
        Parameters  = θ
    elseif DistributionName == "gamma"
        # Extract moments:
        Mean        = Moments[1]
        STD         = Moments[2]

        # Convert moments to parameters:
        α           = Mean^2 / STD^2
        θ           = STD^2 / Mean
        Parameters  = [α, θ]
    elseif DistributionName == "gumbel"
        # Extract moments:
        Mean        = Moments[1]
        STD         = Moments[2]

        # Convert moments to parameters:
        γ           = Base.MathConstants.eulergamma
        μ           = Mean - (STD * γ * sqrt(6)) / π
        θ           = (STD * sqrt(6)) / π
        Parameters  = [μ, θ]
    elseif DistributionName == "lognormal"
        # Extract moments:
        Mean        = Moments[1]
        STD         = Moments[2]

        # Convert moments to parameters:
        μ           = log(Mean) - log(sqrt(1 + (STD / Mean)^2))
        σ           = sqrt(log(1 + (STD / Mean)^2))
        Parameters  = [μ, σ]
    elseif DistributionName == "normal"
        # Extract moments:
        Mean        = Moments[1]
        STD         = Moments[2]

        # Convert moments to parameters:
        μ           = Mean
        σ           = STD
        Parameters  = [μ, σ]
    elseif DistributionName == "poisson"
        # Extract moments:
        Mean        = Moments[1]

        # Convert moments to parameters:
        λ           = Mean
        Parameters  = λ
    elseif DistributionName == "uniform"
        # Extract moments:
        Mean        = Moments[1]
        STD         = Moments[2]

        # Convert moments to parameters:
        a           = Mean - STD * sqrt(3)
        b           = Mean + STD * sqrt(3)
        Parameters  = [a, b]
    elseif DistributionName == "weibull"
        # Extract moments:
        Mean        = Moments[1]
        STD         = Moments[2]

        # Convert moments to parameters:
        F(u, p)     = sqrt(SpecialFunctions.gamma(1 + 2 / u) - SpecialFunctions.gamma(1 + 1 / u)^2) / SpecialFunctions.gamma(1 + 1 / u) - STD / Mean
        u₀          = (0.1, 1000)
        Problem     = NonlinearSolve.IntervalNonlinearProblem(F, u₀)
        Solution    = NonlinearSolve.solve(Problem, NonlinearSolve.Bisection(), abstol=10^(-9), reltol=10^(-9))
        α           = Solution.u
        θ           = Mean / SpecialFunctions.gamma(1 + 1 / α)
        Parameters  = [α, θ]
    else
        error("Provided distribution is not supported.")
    end

    # Return the result:
    return Parameters
end