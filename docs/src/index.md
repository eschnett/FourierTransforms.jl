# FFT.jl

Fast Fourier Transforms implemented in pure Julia.

The main goal of this package is to calculate Fourier Transforms for
all Julia array types and size and element types with reasonable
efficiency (i.e. with `O(n log n)` operations and `O(n)` storage).
This should also make this package differentiable. You can expect
special-purpose packages that specialize on array or element types
such as [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) to achieve a
higher performance.

Currently, not all input lengths are supported. Powers of 2 and
products of powers small primes work fine. However, each prime factor
of the input length is handled via a direct Fourier transform, and
this becomes inefficient if there are large prime factors.

Most functions come in two versions, one that mutates its arguments
and one that allocates its output.

## Example

The Fourier Transform is a linear operator:
```Julia
using FFT
Random.seed!(100)

T = Complex{BigFloat}
N = 768

# Create some data
x = T.(randn(Complex{Float64}, N));
y = T.(randn(Complex{Float64}, N));
α = T.(randn(Complex{Float64}));

# Calculate a linear combination
z = α * x + y;

# Evaluate the Fourier Transform
X = fft(x);
Y = fft(y);
Z = fft(z);

# Check linearity
println("Error: ", maximum(abs.(α * X + Y - Z)))
```

## Fourier Transforms

The functions with `dh` in their names use the Driscoll & Healy
quadrature points. With `L = lmax+1` modes there are `nphi * ntheta`
quadrature points on the sphere, with `nphi = 2L-1` and `ntheta = 2L`.
The points are equispaced in the angles `theta` (latitude) and `phi`
(longitude), and they straddle (avoid) the poles.

```@docs
fft!
fft
inv_fft!
inv_fft
```
