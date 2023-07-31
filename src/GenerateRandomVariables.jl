function generaterv(DistributionName::String, DefineBy::String, Values::Union{Real,Vector{<:Real}})
    # Convert strings to lowercase:
    DistributionName = lowercase(DistributionName)
    DefineBy = lowercase(DefineBy)

    # Error-catching:
    if DefineBy != "moments" && DefineBy != "m" && DefineBy != "parameters" && DefineBy != "p"
        error("""Random variables can only be defined by "Moments" ("M") and "Parameters ("P")".""")
    end

    # Convert moments of a random variable into its parameters if the random variable is defined by its moments:
    if DefineBy == "moments" || DefineBy == "m"
        Values = convertmoments(DistributionName, Values)
    end

    # Create a random variable:
    if DistributionName == "exponential"
        RV = Exponential(Values)
    elseif DistributionName == "gamma"
        RV = Gamma(Values[1], Values[2])
    elseif DistributionName == "gumbel"
        RV = Gumbel(Values[1], Values[2])
    elseif DistributionName == "lognormal"
        RV = LogNormal(Values[1], Values[2])
    elseif DistributionName == "normal"
        RV = Normal(Values[1], Values[2])
    elseif DistributionName == "poisson"
        RV = Poisson(Values[1])
    elseif DistributionName == "uniform"
        RV = Uniform(Values[1], Values[2])
    else
        error("Provided distribution is not supported.")
    end

    return RV
end

function convertmoments(DistributionName::String, Moments)
    # Convert strings to lowercase:
    DistributionName = lowercase(DistributionName)

    # Error catching:
    if length(Moments) > 2
        error("Too many moments are provided.")
    end

    if DistributionName == "exponential"
        # Extract moments:
        Mean = Moments[1]
        STD = Moments[2]

        # Error catching:
        if Mean != STD
            error("Mean and standard deviation values of an exponential random variables must be the same.")
        end

        # Convert moments to parameters:
        θ = Mean
        Parameters = θ
    elseif DistributionName == "gamma"
        # Extract moments:
        Mean = Moments[1]
        STD = Moments[2]

        # Convert moments to parameters:
        α = Mean^2 / STD^2
        θ = STD^2 / Mean
        Parameters = [α, θ]
    elseif DistributionName == "gumbel"
        # Extract moments:
        Mean = Moments[1]
        STD = Moments[2]

        # Convert moments to parameters:
        γ = MathConstants.eulergamma
        μ = Mean - STD * γ * sqrt(6) / π
        θ = STD * sqrt(6) / π
        Parameters = [μ, θ]
    elseif DistributionName == "lognormal"
        # Extract moments:
        Mean = Moments[1]
        STD = Moments[2]

        # Convert moments to parameters:
        μ = log(Mean) - log(sqrt(1 + (STD / Mean)^2))
        σ = sqrt(log(1 + (STD / Mean)^2))
        Parameters = [μ, σ]
    elseif DistributionName == "normal"
        # Extract moments:
        Mean = Moments[1]
        STD = Moments[2]

        # Convert moments to parameters:
        μ = Mean
        σ = STD
        Parameters = [μ, σ]
    elseif DistributionName == "poisson"
        # Extract moments:
        Mean = Moments[1]

        # Convert moments to parameters:
        λ = Mean
        Parameters = λ
    elseif DistributionName == "uniform"
        # Extract moments:
        Mean = Moments[1]
        STD = Moments[2]

        # Convert moments to parameters:
        a = Mean - STD * sqrt(3)
        b = Mean + STD * sqrt(3)
        Parameters = [a, b]
    else
        error("Provided distribution is not supported.")
    end

    return Parameters
end