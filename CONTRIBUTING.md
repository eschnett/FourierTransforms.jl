# Contributing

Contributions are welcome!

Please open an issue if you have a suggestion or found an error, or
want to discuss a larger change before opening a pull request. You are
also welcome to open a pull request if you have a concrete
modification.

## Planned improvements (looking for help)

This package is still young, and is in need of several particular
improvements:

- Benchmark results
- Improve efficiency for small transforms (loop unrolling, SIMD
  operations)
- Investigate transposing arrays to improve memory access efficiency
- Automatically fall back to efficient packages (e.g. FFTW.jl) if
  possible
- Investigate various radix choices
- Implement a true in-place transform, or reduce memory requirements
- Handle large prime factors more efficiently
- Implement an operator interface, i.e. make a Fourier transform look
  like an abstract matrix
- Support 2d and higher-dimensional transforms
- Support distributed computing, e.g. DistributedArrays.jl
- Support GPUs

## Details

The main goal of this package is to calculate Fourier Transforms for
all Julia array types and size and element types with reasonable
efficiency (i.e. with `O(n log n)` operations and `O(n)` storage).

Code additions or changes should be accompanied by test cases.
Improvements that work only for special cases should be accompanied by
generic code that works on all architectures and for all types.
