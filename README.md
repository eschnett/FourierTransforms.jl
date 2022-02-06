# FourierTransforms.jl

Fast Fourier Transforms implemented in pure Julia.

* [![Documenter](https://img.shields.io/badge/docs-dev-blue.svg)](https://eschnett.github.io/FourierTransforms.jl/dev)
* [![GitHub
  CI](https://github.com/eschnett/FourierTransforms.jl/workflows/CI/badge.svg)](https://github.com/eschnett/FourierTransforms.jl/actions)
* [![Codecov](https://codecov.io/gh/eschnett/FourierTransforms.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eschnett/FourierTransforms.jl)

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

## Example

The Fourier Transform is a linear operator:
```Julia
using FourierTransforms
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

## Related Packages

- [FFTW.jl](https://github.com/JuliaMath/FFTW.jl)
- [FourierTransforms.jl](https://github.com/JuliaComputing/FourierTransforms.jl)
  (unregistered, apparently unmaintained)
