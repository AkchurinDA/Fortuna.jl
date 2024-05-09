@testset "Sampling Techniques #1" begin
    # Set an RNG seed:
    Random.seed!(123)

    # Generate a random vector X with uncorrelated marginal random variables X₁ and X₂:
    X₁ = randomvariable("Gamma", "M", [10, 1.5])
    X₂ = randomvariable("Gamma", "M", [15, 2.5])
    X  = [X₁, X₂]

    # Generate samples:
    NumSamples  = 10 ^ 6
    XSamplesITS = rand(X, NumSamples, ITS())
    XSamplesLHS = rand(X, NumSamples, LHS())

    # Test the results:
    @test isapprox(mean(XSamplesITS[1, :]),    10,         rtol = 1E-2) # Inverse Transform Sampling
    @test isapprox(mean(XSamplesITS[2, :]),    15,         rtol = 1E-2)
    @test isapprox(std(XSamplesITS[1, :]),     1.5,        rtol = 1E-2)
    @test isapprox(std(XSamplesITS[2, :]),     2.5,        rtol = 1E-2)
    @test isapprox(cor(XSamplesITS, dims = 2), [1 0; 0 1], rtol = 1E-2)
    @test isapprox(mean(XSamplesLHS[1, :]),    10,         rtol = 1E-2) # Latin Hypercube Sampling
    @test isapprox(mean(XSamplesLHS[2, :]),    15,         rtol = 1E-2)
    @test isapprox(std(XSamplesLHS[1, :]),     1.5,        rtol = 1E-2)
    @test isapprox(std(XSamplesLHS[2, :]),     2.5,        rtol = 1E-2)
    @test isapprox(cor(XSamplesLHS, dims = 2), [1 0; 0 1], rtol = 1E-2)
end

@testset "Sampling Techniques #2" begin
    # Set an RNG seed:
    Random.seed!(123)

    # Define a list of reliability indices of interest:
    ρList = (-0.75):(0.25):(+0.75)

    for i in eachindex(ρList)
        # Generate a random vector X of correlated marginal distributions:
        X₁ = randomvariable("Gamma", "M", [10, 1.5])
        X₂ = randomvariable("Gamma", "M", [15, 2.5])
        X  = [X₁, X₂]
        ρˣ = [1 ρList[i]; ρList[i] 1]

        # Perform Nataf transformation of the correlated marginal random variables:
        NatafObject = NatafTransformation(X, ρˣ)

        # Generate samples:
        NumSamples        = 10 ^ 6
        XSamplesITS, _, _ = rand(NatafObject, NumSamples, ITS())
        XSamplesLHS, _, _ = rand(NatafObject, NumSamples, LHS())

        # Test the results:
        @test isapprox(mean(XSamplesITS[1, :]),    10,  rtol = 1E-2) # Inverse Transform Sampling
        @test isapprox(mean(XSamplesITS[2, :]),    15,  rtol = 1E-2)
        @test isapprox(std(XSamplesITS[1, :]),     1.5, rtol = 1E-2)
        @test isapprox(std(XSamplesITS[2, :]),     2.5, rtol = 1E-2)
        @test isapprox(cor(XSamplesITS, dims = 2), ρˣ,  rtol = 1E-2)
        @test isapprox(mean(XSamplesLHS[1, :]),    10,  rtol = 1E-2) # Latin Hypercube Sampling
        @test isapprox(mean(XSamplesLHS[2, :]),    15,  rtol = 1E-2)
        @test isapprox(std(XSamplesLHS[1, :]),     1.5, rtol = 1E-2)
        @test isapprox(std(XSamplesLHS[2, :]),     2.5, rtol = 1E-2)
        @test isapprox(cor(XSamplesLHS, dims = 2), ρˣ,  rtol = 1E-2)
    end
end

@testset "Sampling Techniques #3" begin
    # Generate a random vector X of correlated marginal distributions:
    X₁ = randomvariable("Gamma", "M", [10, 1.5])
    X₂ = randomvariable("Gamma", "M", [15, 2.5])
    X  = [X₁, X₂]
    ρˣ = [1 0.75; 0.75 1]

    # Perform Nataf transformation of the correlated marginal random variables:
    NatafObject = NatafTransformation(X, ρˣ)

    # Define number of samples:
    NumSamples = 10 ^ 6

    # Generate samples from a random variable:
    Random.seed!(123)
    XSamplesITS₁ = rand(X₁, NumSamples, ITS())
    XSamplesLHS₁ = rand(X₁, NumSamples, LHS())

    Random.seed!(123)
    XSamplesITS₂ = rand(X₁, NumSamples, ITS())
    XSamplesLHS₂ = rand(X₁, NumSamples, LHS())

    # Test the results:
    @test XSamplesITS₁ == XSamplesITS₂
    @test XSamplesLHS₁ == XSamplesLHS₂

    # Generate samples from a random variable:
    Random.seed!(123)
    XSamplesITS₁ = rand(X₂, NumSamples, ITS())
    XSamplesLHS₁ = rand(X₂, NumSamples, LHS())

    Random.seed!(123)
    XSamplesITS₂ = rand(X₂, NumSamples, ITS())
    XSamplesLHS₂ = rand(X₂, NumSamples, LHS())

    # Test the results:
    @test XSamplesITS₁ == XSamplesITS₂
    @test XSamplesLHS₁ == XSamplesLHS₂

    # Generate samples from a random vector:
    Random.seed!(123)
    XSamplesITS₁ = rand(X, NumSamples, ITS())
    XSamplesLHS₁ = rand(X, NumSamples, LHS())

    Random.seed!(123)
    XSamplesITS₂ = rand(X, NumSamples, ITS())
    XSamplesLHS₂ = rand(X, NumSamples, LHS())

    # Test the results:
    @test XSamplesITS₁ == XSamplesITS₂
    @test XSamplesLHS₁ == XSamplesLHS₂

    # Generate samples from a transformation object:
    Random.seed!(123)
    XSamplesITS₁, ZSamplesITS₁, USamplesITS₁ = rand(NatafObject, NumSamples, ITS())
    XSamplesLHS₁, ZSamplesLHS₁, USamplesLHS₁ = rand(NatafObject, NumSamples, LHS())

    Random.seed!(123)
    XSamplesITS₂, ZSamplesITS₂, USamplesITS₂ = rand(NatafObject, NumSamples, ITS())
    XSamplesLHS₂, ZSamplesLHS₂, USamplesLHS₂ = rand(NatafObject, NumSamples, LHS())

    # Test the results:
    @test XSamplesITS₁ == XSamplesITS₂
    @test ZSamplesITS₁ == ZSamplesITS₂
    @test USamplesITS₁ == USamplesITS₂
    @test XSamplesLHS₁ == XSamplesLHS₂
    @test ZSamplesLHS₁ == ZSamplesLHS₂
    @test USamplesLHS₁ == USamplesLHS₂
end