# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using FourierTransforms

makedocs(; sitename="FourierTransforms", format=Documenter.HTML(), modules=[FourierTransforms])

deploydocs(; repo="github.com/eschnett/FourierTransforms.jl.git", devbranch="main", push_preview=true)
