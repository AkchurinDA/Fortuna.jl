@testset "Nataf Transformation: Distorted Correlation Matrix #1 (Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector:
    X₁  = randomvariable("Normal", "M", [0, 1])
    X₂  = randomvariable("Normal", "M", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Test the results:
    ρᶻ = NatafObject.ρᶻ
    @test ρᶻ == ρˣ
end

@testset "Nataf Transformation: Distorted Correlation Matrix #2 (Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector:
    X₁  = randomvariable("Uniform", "P", [0, 1])
    X₂  = randomvariable("Uniform", "P", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Test the results:
    ρᶻ = NatafObject.ρᶻ
    @test ρᶻ == ρˣ
end

@testset "Nataf Transformation: Distorted Correlation Matrix #3 (Non-Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector:
    X₁  = randomvariable("Normal", "M", [0, 1])
    X₂  = randomvariable("Normal", "M", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0.8; 0.8 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Test the results:
    ρᶻ = NatafObject.ρᶻ
    @test isapprox(ρᶻ, [1 0.8; 0.8 1], rtol = 10^(-6))
end

@testset "Nataf Transformation: Distorted Correlation Matrix #4 (Non-Identity)" begin
    # Test from UQPy package (https://github.com/SURGroup/UQpy/tree/master)

    # Define a random vector:
    X₁  = randomvariable("Uniform", "P", [0, 1])
    X₂  = randomvariable("Uniform", "P", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0.8; 0.8 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Test the results:
    ρᶻ = NatafObject.ρᶻ
    @test isapprox(ρᶻ, [1 0.8134732861515996; 0.8134732861515996 1], rtol = 10 ^ (-6))
end

@testset "Nataf Transformation: Transformed Samples #1 (Identity)" begin
    # Define a random vector:
    X₁  = randomvariable("Normal", "M", [0, 1])
    X₂  = randomvariable("Normal", "M", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Define samples:
    XRange₁ = range(-6, 6, length = 100)
    XRange₂ = range(-6, 6, length = 100)
    Samples = Matrix{Float64}(undef, 2, 100 * 100)
    for i in 1:100
        for j in 1:100
            Samples[1, (i - 1) * 100 + j] = XRange₁[i]
            Samples[2, (i - 1) * 100 + j] = XRange₂[j]
        end
    end

    # Perform transformation:
    TransformedSamplesX2U = transformsamples(NatafObject, Samples, "X2U")
    TransformedSamplesU2X = transformsamples(NatafObject, Samples, "U2X")

    # Test the results:
    @test TransformedSamplesX2U == TransformedSamplesU2X
end

@testset "Nataf Transformation: Jacobians #1 (Identity)" begin
    # Define a random vector:
    X₁  = randomvariable("Normal", "M", [0, 1])
    X₂  = randomvariable("Normal", "M", [0, 1])
    X   = [X₁, X₂]
    ρˣ  = [1 0; 0 1]

    # Perform Nataf transformation:
    NatafObject = NatafTransformation(X, ρˣ)

    # Define samples:
    XRange₁ = range(-6, 6, length = 100)
    XRange₂ = range(-6, 6, length = 100)
    Samples = Matrix{Float64}(undef, 2, 100 * 100)
    for i in 1:100
        for j in 1:100
            Samples[1, (i - 1) * 100 + j] = XRange₁[i]
            Samples[2, (i - 1) * 100 + j] = XRange₂[j]
        end
    end

    # Perform transformation:
    JX2U = getjacobian(NatafObject, Samples, "X2U")
    JU2X = getjacobian(NatafObject, Samples, "U2X")

    # Test the results:
    @test JX2U == JU2X
end