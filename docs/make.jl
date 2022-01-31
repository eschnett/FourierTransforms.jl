# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using FFT

makedocs(; sitename="FFT", format=Documenter.HTML(), modules=[FFT])

deploydocs(; repo="github.com/eschnett/FFT.jl.git", devbranch="main", push_preview=true)
