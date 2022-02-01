using Compat: cispi
# using DoubleFloats
using FFT
using Random
using Test

# `DoubleFloat` does not work: <https://github.com/JuliaMath/DoubleFloats.jl/issues/139>
const Types = [Complex{Float32}, Complex{Float64}, Complex{BigFloat}]

# # Random.rand(rng::AbstractRNG, ::Random.Random.Random.CloseOpen01{DoubleFloat}) = ???
# Random.randn(::Type{DoubleFloat}, dims::Integer...) = DoubleFloat.(randn(Float64, dims...))
# Random.randn(::Type{Complex{DoubleFloat}}, dims::Integer...) = Complex{DoubleFloat}.(randn(Complex{Float64}, dims...))
# Base.eps(::Type{DoubleFloat}) = eps(Float64)^2

# Random.rand(rng::AbstractRNG, ::Random.Random.Random.CloseOpen01{BigFloat}) = ???
Random.randn(::Type{BigFloat}, dims::Integer...) = BigFloat.(randn(Float64, dims...))
Random.randn(::Type{Complex{BigFloat}}, dims::Integer...) = Complex{BigFloat}.(randn(Complex{Float64}, dims...))

Random.seed!(100)
@testset "Direct Fourier transforms (T=$T)" for T in Types
    maxbits = T === Complex{BigFloat} ? 8 : 10
    for bits in 1:maxbits
        N = 2^bits
        x = randn(T, N)
        X = FFT.direct_ft(x)
        x′ = conj(FFT.direct_ft(conj(X)) / N)
        @test x′ ≈ x
    end
end

Random.seed!(100)
@testset "Radix-2 Cooley–Tukey Fast Fourier transforms (T=$T)" for T in Types
    maxbits = T === Complex{BigFloat} ? 15 : 20
    for bits in 1:maxbits
        N = 2^bits
        x = randn(T, N)
        X = FFT.ditfft2(x)
        x′ = conj(FFT.ditfft2(conj(X)) / N)
        @test x′ ≈ x
    end
end

Random.seed!(100)
@testset "General radix Cooley–Tukey Fast Fourier transforms (T=$T)" for T in Types
    maxbits = T === Complex{BigFloat} ? 15 : 20
    for bits in 1:maxbits
        N = 2^bits
        x = randn(T, N)
        X = FFT.radix_fft(x; radix=2)
        x′ = conj(FFT.radix_fft(conj(X); radix=2) / N)
        @test x′ ≈ x
    end
end

Random.seed!(100)
@testset "Comparing algorithms for Fourier transforms (T=$T)" for T in Types
    maxbits = T === Complex{BigFloat} ? 8 : 10
    for bits in 1:maxbits
        N = 2^bits
        x = randn(T, N)
        X = FFT.direct_ft(x)
        if N == 1
            X1 = FFT.direct_ft_1(x)
            @test X1 ≈ X
        end
        if N == 2
            X2 = FFT.direct_ft_2(x)
            @test X2 ≈ X
        end
        if N == 4
            X4 = FFT.direct_ft_4(x)
            @test X4 ≈ X
        end
        Xct2 = FFT.ditfft2(x)
        @test Xct2 ≈ X
        for radix in 2:8
            if N % radix == 0
                Xr = FFT.radix_fft(x; radix=radix)
                @test Xr ≈ X
            end
        end
    end
end

Random.seed!(100)
@testset "Simple Fourier transforms (T=$T)" for T in Types
    RT = typeof(real(zero(T)))
    maxbits = T === Complex{BigFloat} ? 15 : 20
    for bits in 1:maxbits
        for iter in 1:10
            N = 2^bits
            k = rand(1:N)
            x = [cispi(2 * ((k - 1) * (n - 1) % N) / RT(N)) for n in 1:N]
            X₀ = [(k′ == k) * N for k′ in 1:N]
            X = fft(x)
            @test isapprox(X, X₀; atol=N^2 * sqrt(eps()))
        end
    end
end

Random.seed!(100)
@testset "Inverse Fourier transforms (T=$T)" for T in Types
    maxbits = T === Complex{BigFloat} ? 15 : 20
    for bits in 1:maxbits
        for iter in 1:10
            N = 2^bits
            x = randn(T, N)
            X = fft(x)
            x′ = inv_fft(X)
            @test x′ ≈ x
        end
    end
end

Random.seed!(100)
@testset "Linearity of Fourier transforms (T=$T)" for T in Types
    maxbits = T === Complex{BigFloat} ? 10 : 15
    for bits in 1:maxbits
        for iter in 1:10
            N = 2^bits
            x = randn(T, N)
            y = randn(T, N)
            α = randn(T)
            z = α * x + y
            X = fft(x)
            Y = fft(y)
            Z = fft(z)
            @test α * X + Y ≈ Z
        end
    end
end
