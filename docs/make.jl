using Documenter
using ControlTimingSafety

makedocs(
    sitename = "ControlTimingSafety",
    format = Documenter.HTML(),
    modules = [ControlTimingSafety]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/Ratfink/ControlTimingSafety.jl.git"
)
